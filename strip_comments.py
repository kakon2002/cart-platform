"""Remove every comment and reduce each docstring to one line."""

from __future__ import annotations

import argparse
import ast
import io
import sys
import tokenize
from pathlib import Path

ROOT = Path(__file__).resolve().parent

SKIP_ANYWHERE = {".venv", ".git", "__pycache__"}

SKIP_PATHS = ("data", "reports")


def targets() -> list[Path]:
    """Every Python file in the project except vendored and generated trees."""
    skip_roots = [(ROOT / name).resolve() for name in SKIP_PATHS]
    out = []
    for path in ROOT.rglob("*.py"):
        if any(part in SKIP_ANYWHERE for part in path.parts):
            continue

        resolved = path.resolve()
        if any(resolved.is_relative_to(root) for root in skip_roots):
            continue
        out.append(path)
    return sorted(out)


def strip_comments(source: str) -> str:
    """Return the source with comment tokens removed."""
    result: list[str] = []
    last_row, last_col = 1, 0
    for tok in tokenize.generate_tokens(io.StringIO(source).readline):
        kind, text, (srow, scol), (erow, ecol), _line = tok
        if srow > last_row:
            result.append("\n" * (srow - last_row))
            last_col = 0
        if scol > last_col:
            result.append(" " * (scol - last_col))
        if kind == tokenize.COMMENT:
            last_row, last_col = erow, ecol
            continue

        if kind in (tokenize.NEWLINE, tokenize.NL):
            last_row, last_col = erow, ecol
            continue

        if kind == getattr(tokenize, "FSTRING_MIDDLE", None):
            text = text.replace("{", "{{").replace("}", "}}")
            result.append(text)

            last_row = erow
            last_col = (scol if erow == srow else 0) + len(text.split(chr(10))[-1])
            continue
        result.append(text)
        last_row, last_col = erow, ecol
    return "".join(result)


def summary_line(text: str) -> str:
    """The docstring's opening paragraph, joined into a single line."""
    out: list[str] = []
    for line in text.strip().splitlines():
        if not line.strip():
            break
        out.append(line.strip())
    return " ".join(out)


def shorten_docstrings(source: str) -> tuple[str, list[str]]:
    """Reduce every docstring to its opening paragraph, joined onto one line."""
    tree = ast.parse(source)
    lines = source.splitlines(keepends=True)
    edits: list[tuple[int, int, str]] = []
    joined: dict[int, str] = {}

    def visit(node):
        """Record the docstring edit for this node and its nested definitions."""
        body = getattr(node, "body", None)
        if not isinstance(body, list) or not body:
            return
        if isinstance(body[0], ast.Expr) and isinstance(
                body[0].value, ast.Constant) and isinstance(body[0].value.value, str):
            expr = body[0]
            raw = expr.value.value.strip()
            summary = summary_line(raw)
            paragraph = []
            for line in raw.splitlines():
                if not line.strip():
                    break
                paragraph.append(line)
            if len(paragraph) > 1:
                joined[expr.lineno] = f"line {expr.lineno}: {summary}"
            indent = " " * expr.col_offset
            edits.append((expr.lineno, expr.end_lineno,
                          f'{indent}"""{summary}"""\n'))
        for child in body:
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef,
                                  ast.ClassDef)):
                visit(child)

    visit(tree)
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            visit(node)

    seen: set[int] = set()
    for start, end, replacement in sorted(edits, key=lambda e: -e[0]):
        if start in seen:
            continue
        seen.add(start)
        lines[start - 1:end] = [replacement]
    return "".join(lines), [joined[k] for k in sorted(joined)]


def process(path: Path, write: bool) -> tuple[int, int, list[str]]:
    """Strip one file; return (comments removed, lines removed, summaries joined)."""
    original = path.read_text(encoding="utf-8")
    before_comments = sum(
        1 for tok in tokenize.generate_tokens(io.StringIO(original).readline)
        if tok.type == tokenize.COMMENT)

    stripped, joined = shorten_docstrings(original)
    stripped = strip_comments(stripped)

    out_lines: list[str] = []
    pending = 0
    for line in stripped.splitlines():
        if not line.strip():
            pending += 1
            continue

        allowed = 2 if line[:1] not in (" ", chr(9)) else 1

        prev = out_lines[-1].rstrip() if out_lines else ""
        if prev.endswith((":", "(", "[", "{", ",")) or                 line.lstrip()[:1] in (")", "]", "}"):
            allowed = 0
        if out_lines:
            out_lines.extend([""] * min(pending, allowed))
        pending = 0
        out_lines.append(line.rstrip())
    text = "\n".join(out_lines).rstrip() + "\n"

    a = ast.dump(ast.parse(original))
    b = ast.dump(ast.parse(text))
    if a != b:
        a = ast.dump(_drop_docstrings(ast.parse(original)))
        b = ast.dump(_drop_docstrings(ast.parse(text)))
        if a != b:
            raise SystemExit(
                f"{path}: the stripped file does not parse to the same tree. "
                "Refusing to write - this would be a behaviour change.")

    delta = len(original.splitlines()) - len(text.splitlines())
    if write:
        path.write_text(text, encoding="utf-8")
    return before_comments, delta, joined


def _drop_docstrings(tree: ast.AST) -> ast.AST:
    """Blank every docstring so two trees compare on code alone."""
    for node in ast.walk(tree):
        body = getattr(node, "body", None)

        if not isinstance(body, list) or not body:
            continue
        if isinstance(body[0], ast.Expr) and isinstance(
                body[0].value, ast.Constant) and isinstance(body[0].value.value, str):
            body[0].value.value = ""
    return tree


def main() -> int:
    """Strip every file, or report what would change."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true",
                        help="report what would change and write nothing")
    args = parser.parse_args()

    total_comments = total_lines = 0
    rewrapped: list[tuple[str, str]] = []
    for path in targets():
        rel = path.relative_to(ROOT).as_posix()
        comments, delta, joined = process(path, write=not args.check)
        total_comments += comments
        total_lines += delta
        rewrapped.extend((rel, note) for note in joined)
        if comments or delta:
            print(f"  {rel:44s} {comments:4d} comments, {delta:4d} lines")
    verb = "would remove" if args.check else "removed"
    print(f"\n{verb} {total_comments} comments and {total_lines} lines "
          f"across {len(targets())} files")

    print(f"{len(rewrapped)} wrapped summaries joined onto one line; the "
          "AST guard cannot see docstrings, so this list is the only report "
          "of what was rewritten")
    for rel, note in rewrapped:
        print(f"  {rel} {note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
