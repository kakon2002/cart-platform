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


def first_line(text: str) -> str:
    """The first non-empty line of a docstring, trimmed."""
    for line in text.strip().splitlines():
        if line.strip():
            return line.strip()
    return ""


def shorten_docstrings(source: str) -> str:
    """Reduce every module, class and function docstring to its first line."""
    tree = ast.parse(source)
    lines = source.splitlines(keepends=True)
    edits: list[tuple[int, int, str]] = []

    def visit(node):
        """Record the docstring edit for this node and its nested definitions."""
        body = getattr(node, "body", None)
        if not isinstance(body, list) or not body:
            return
        if isinstance(body[0], ast.Expr) and isinstance(
                body[0].value, ast.Constant) and isinstance(body[0].value.value, str):
            expr = body[0]
            summary = first_line(expr.value.value)
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
    return "".join(lines)


def process(path: Path, write: bool) -> tuple[int, int]:
    """Strip one file; return (comment lines removed, docstring lines removed)."""
    original = path.read_text(encoding="utf-8")
    before_comments = sum(
        1 for tok in tokenize.generate_tokens(io.StringIO(original).readline)
        if tok.type == tokenize.COMMENT)

    stripped = shorten_docstrings(original)
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
    return before_comments, delta


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
    for path in targets():
        comments, delta = process(path, write=not args.check)
        total_comments += comments
        total_lines += delta
        if comments or delta:
            print(f"  {path.relative_to(ROOT).as_posix():44s} "
                  f"{comments:4d} comments, {delta:4d} lines")
    verb = "would remove" if args.check else "removed"
    print(f"\n{verb} {total_comments} comments and {total_lines} lines "
          f"across {len(targets())} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
