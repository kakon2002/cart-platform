"""Render the full-run report as a single page.

Reads ``reports/full-run.json`` — the structured result ``run_all.py`` writes
straight from the objects it already holds. Deliberately not the markdown: a
regex pinned to another file's f-strings breaks when that file is reformatted,
and worse, can half-match and render confident headings above empty lists.

    .venv\\Scripts\\python.exe run_all.py --fresh
    .venv\\Scripts\\python.exe make_artifact.py [reports/full-run.json]
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

#: Stages 7 and 8 are absent from the pipeline, not from this page. Showing the
#: gap is information: a reader counting 1..11 and finding nine entries should
#: be told why rather than left to wonder whether two were dropped from the
#: report.
ABSENT = {
    "7": ("Manufacturing", "schema only; not on the path to the end state"),
    "8": ("Trial design", "schema only; not on the path to the end state"),
}


def attrition(api_lines: list[str]) -> list[dict]:
    """The attrition chain, read back out of the API's own closing block."""
    steps = []
    for line in api_lines:
        hit = re.match(r"^\s+(.+?)\s+-\s*(\d+)\s+(\d+) remain$", line)
        if hit:
            steps.append({"gate": hit.group(1).strip(),
                          "dropped": int(hit.group(2)),
                          "remaining": int(hit.group(3))})
    return steps


def statuses(api_lines: list[str]) -> list[tuple[str, str, list[str]]]:
    """Each endpoint's named status and its reasons."""
    out: list[tuple[str, str, list[str]]] = []
    for line in api_lines:
        hit = re.match(r"^\s*(GET \S+)\s*->\s*(\S+)$", line)
        if hit:
            out.append((hit.group(1), hit.group(2), []))
        elif out and line.strip().startswith("- "):
            out[-1][2].append(line.strip()[2:])
    return out


E = html.escape


def render(run: dict) -> str:
    steps = attrition(run["api"])
    ends = statuses(run["api"])
    pool = steps[0]["dropped"] + steps[0]["remaining"] if steps else 0

    ledger = []
    for stage in run["stages"]:
        state = "clear" if stage["ok"] else "tripped"
        mark = "clear" if stage["ok"] else f"{stage['total'] - stage['clear']} tripped"
        ledger.append(f"""
        <tr class="row row--{state}">
          <td class="num">{E(stage['number'])}</td>
          <td class="nm">{E(stage['name'])}</td>
          <td class="ct"><span class="frac">{stage['clear']}/{stage['total']}</span></td>
          <td class="st"><span class="pill pill--{state}">{E(mark)}</span></td>
          <td class="tm">{stage['seconds']}s</td>
        </tr>""")
        if stage["number"] == "6":
            for number, (name, why) in ABSENT.items():
                ledger.append(f"""
        <tr class="row row--absent">
          <td class="num">{number}</td>
          <td class="nm">{E(name)}</td>
          <td class="ct">—</td>
          <td class="st"><span class="pill pill--absent">not implemented</span></td>
          <td class="tm">{E(why)}</td>
        </tr>""")

    chain = []
    for index, step in enumerate(steps):
        width = (step["dropped"] / pool * 100) if pool else 0
        chain.append(f"""
        <li class="gate">
          <div class="gate__head">
            <span class="gate__name">{E(step['gate'])}</span>
            <span class="gate__drop">&minus;{step['dropped']}</span>
          </div>
          <div class="gate__bar"><span style="width:{width:.4f}%"></span></div>
          <div class="gate__rem">{step['remaining']} remain</div>
        </li>""")

    endings = []
    for endpoint, status, reasons in ends:
        items = "".join(f"<li>{E(r)}</li>" for r in reasons)
        endings.append(f"""
        <div class="ending">
          <p class="ending__ep">{E(endpoint)}</p>
          <p class="ending__st">{E(status)}</p>
          {f'<ul class="ending__why">{items}</ul>' if items else ''}
        </div>""")

    # A tripped criterion with a decision on record and one without are not the
    # same thing, and rendering them identically is exactly the silencing this
    # page claims not to do.
    accepted = run.get("accepted", {})
    sections = []
    for stage in run["stages"]:
        rows = []
        for c in stage["criteria"]:
            why = accepted.get(f"{stage['number']}/{c['id']}")
            if not c["tripped"]:
                kind, note = "clear", ""
            elif why:
                kind = "accepted"
                note = f'<span class="crit__note">Accepted — {E(why)}</span>'
            else:
                kind = "regression"
                note = ('<span class="crit__note crit__note--new">Not on the '
                        'accepted list. This is new.</span>')
            rows.append(f"""
          <li class="crit crit--{kind}">
            <span class="crit__id">{E(c['id'])}</span>
            <span class="crit__detail">{E(c['detail'])}{note}</span>
          </li>""")
        rows = "".join(rows)
        label = "API" if stage["number"] == "API" else f"Stage {stage['number']}"
        sections.append(f"""
      <section class="stage" id="s{E(stage['number'])}">
        <h3 class="stage__h">
          <span class="stage__label">{E(label)}</span>
          {E(stage['name'])}
          <span class="stage__score">{stage['clear']}/{stage['total']}</span>
        </h3>
        <ul class="crits">{rows}</ul>
      </section>""")

    tripped_total = run.get("total", 0) - run.get("clear", 0)
    new = run.get("unexpected", [])
    stale = run.get("stale", [])
    # The lede must not claim every trip is a recorded limitation when one of
    # them is not. This is the sentence the accepted list exists to make true.
    if new:
        verdict = (f"<strong>{len(new)} of them are not on the accepted list "
                   "and are new.</strong> They are marked below.")
    else:
        verdict = ("Every one has a decision on record, shown against it — a "
                   "tripped criterion here is a recorded limitation, not a "
                   "silenced one.")
    if stale:
        verdict += (f" {len(stale)} accepted exemption(s) did not trip on this "
                    "run and should be removed from the list.")

    return f"""<title>Full Run Ledger</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Instrument+Serif:ital@0;1&amp;family=IBM+Plex+Mono:wght@400;500&amp;family=IBM+Plex+Sans:wght@400;500;600&amp;display=swap">
<style>
:root {{
  --paper:#EDEFF1; --surface:#F8F9FA; --sunk:#E3E7EA;
  --ink:#131C2B; --ink-soft:#515E70; --ink-faint:#7A879A;
  --rule:#D2D8DE; --rule-soft:#E1E6EA;
  --madder:#94301D; --madder-soft:#B8593F;
  --clear:#2C6046; --tripped:#8C5A0C;
  --clear-bg:#DFEAE3; --tripped-bg:#F2E7CF; --madder-bg:#F3DFDA;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper:#0D131B; --surface:#161E29; --sunk:#111925;
    --ink:#E3E8EE; --ink-soft:#9AA6B6; --ink-faint:#6E7C8D;
    --rule:#28323F; --rule-soft:#1D2631;
    --madder:#DE7A63; --madder-soft:#C4614A;
    --clear:#71C398; --tripped:#D6A445;
    --clear-bg:#152A21; --tripped-bg:#2A2213; --madder-bg:#331A15;
  }}
}}
:root[data-theme="dark"] {{
  --paper:#0D131B; --surface:#161E29; --sunk:#111925;
  --ink:#E3E8EE; --ink-soft:#9AA6B6; --ink-faint:#6E7C8D;
  --rule:#28323F; --rule-soft:#1D2631;
  --madder:#DE7A63; --madder-soft:#C4614A;
  --clear:#71C398; --tripped:#D6A445;
  --clear-bg:#152A21; --tripped-bg:#2A2213; --madder-bg:#331A15;
}}

* {{ box-sizing:border-box; }}
body {{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  font-size:16px; line-height:1.6;
  -webkit-font-smoothing:antialiased;
}}
.wrap {{ max-width:1080px; margin:0 auto; padding:0 28px 96px; }}
.prose {{ max-width:68ch; }}

/* ---- masthead ---------------------------------------------------------- */
.mast {{
  border-bottom:2px solid var(--ink);
  padding:56px 0 22px; margin-bottom:36px;
  display:flex; flex-wrap:wrap; gap:28px; align-items:flex-end;
  justify-content:space-between;
}}
.mast__t {{
  font-family:"Instrument Serif",Georgia,"Times New Roman",serif;
  font-weight:400; font-size:clamp(2.6rem,6vw,4.1rem); line-height:1.02;
  margin:0; letter-spacing:-0.012em; text-wrap:balance;
}}
.mast__sub {{
  margin:12px 0 0; color:var(--ink-soft); max-width:52ch; font-size:1.02rem;
}}
.tally {{ text-align:right; flex-shrink:0; }}
.tally__n {{
  font-family:"Instrument Serif",Georgia,serif;
  font-size:clamp(3rem,8vw,4.6rem); line-height:0.9; display:block;
  font-variant-numeric:tabular-nums; color:var(--madder);
}}
.tally__l {{
  font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:0.7rem; letter-spacing:0.13em; text-transform:uppercase;
  color:var(--ink-faint); display:block; margin-top:10px;
}}

/* ---- shared ------------------------------------------------------------ */
.eyebrow {{
  font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:0.7rem; letter-spacing:0.15em; text-transform:uppercase;
  color:var(--madder); margin:0 0 6px;
}}
h2 {{
  font-family:"Instrument Serif",Georgia,serif; font-weight:400;
  font-size:2rem; line-height:1.12; margin:0 0 10px; text-wrap:balance;
}}
.block {{ margin-top:64px; }}
.lede {{ color:var(--ink-soft); margin:0 0 24px; max-width:64ch; }}

/* ---- ledger ------------------------------------------------------------ */
.scroll {{ overflow-x:auto; }}
table {{ border-collapse:collapse; width:100%; min-width:560px; }}
thead th {{
  font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:0.68rem; letter-spacing:0.12em; text-transform:uppercase;
  color:var(--ink-faint); font-weight:400; text-align:left;
  padding:0 12px 9px; border-bottom:1px solid var(--rule);
}}
tbody td {{ padding:11px 12px; border-bottom:1px solid var(--rule-soft); }}
.num {{
  font-family:"IBM Plex Mono",ui-monospace,monospace; color:var(--ink-faint);
  width:1%; white-space:nowrap; font-size:0.85rem;
}}
.nm {{ font-weight:500; }}
.frac {{
  font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-variant-numeric:tabular-nums; font-size:0.92rem;
}}
.tm {{
  font-family:"IBM Plex Mono",ui-monospace,monospace; color:var(--ink-faint);
  font-size:0.82rem; font-variant-numeric:tabular-nums; text-align:right;
}}
.row--absent td {{ color:var(--ink-faint); }}
.row--absent .nm {{ font-weight:400; font-style:italic; }}
.row--absent .tm {{ text-align:right; font-style:italic; }}
.pill {{
  display:inline-block; padding:2px 9px; border-radius:2px;
  font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:0.7rem; letter-spacing:0.05em; white-space:nowrap;
}}
.pill--clear {{ background:var(--clear-bg); color:var(--clear); }}
.pill--tripped {{ background:var(--tripped-bg); color:var(--tripped); }}
.pill--absent {{ border:1px solid var(--rule); color:var(--ink-faint); }}

/* ---- attrition --------------------------------------------------------- */
.gates {{ list-style:none; margin:0; padding:0; display:flex;
  flex-direction:column; gap:20px; max-width:640px; }}
.gate__head {{ display:flex; justify-content:space-between; align-items:baseline;
  gap:16px; margin-bottom:6px; }}
.gate__name {{ font-weight:500; }}
.gate__drop {{
  font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-variant-numeric:tabular-nums; color:var(--madder); font-size:0.95rem;
}}
.gate__bar {{ height:9px; background:var(--sunk); border-radius:1px; overflow:hidden; }}
.gate__bar span {{ display:block; height:100%; background:var(--madder);
  min-width:2px; }}
.gate__rem {{
  font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.76rem;
  color:var(--ink-faint); margin-top:5px; font-variant-numeric:tabular-nums;
}}

/* ---- endings ----------------------------------------------------------- */
.endings {{ display:grid; gap:20px; grid-template-columns:repeat(auto-fit,minmax(300px,1fr)); }}
.ending {{
  background:var(--surface); border:1px solid var(--rule);
  border-left:3px solid var(--madder); padding:20px 22px; border-radius:2px;
}}
.ending__ep {{
  font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.76rem;
  color:var(--ink-faint); margin:0 0 8px; letter-spacing:0.04em;
}}
.ending__st {{
  font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.95rem;
  font-weight:500; color:var(--madder); margin:0; word-break:break-word;
}}
.ending__why {{ margin:14px 0 0; padding-left:18px; color:var(--ink-soft);
  font-size:0.9rem; display:flex; flex-direction:column; gap:7px; }}

/* ---- criteria ---------------------------------------------------------- */
.stage {{ margin-top:40px; }}
.stage__h {{
  display:flex; align-items:baseline; gap:12px; flex-wrap:wrap;
  font-family:"Instrument Serif",Georgia,serif; font-weight:400;
  font-size:1.5rem; margin:0 0 14px; padding-bottom:9px;
  border-bottom:1px solid var(--rule);
}}
.stage__label {{
  font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.7rem;
  letter-spacing:0.12em; text-transform:uppercase; color:var(--madder);
}}
.stage__score {{
  margin-left:auto; font-family:"IBM Plex Mono",ui-monospace,monospace;
  font-size:0.85rem; color:var(--ink-faint); font-variant-numeric:tabular-nums;
}}
.crits {{ list-style:none; margin:0; padding:0;
  display:flex; flex-direction:column; gap:2px; }}
.crit {{
  display:grid; grid-template-columns:82px 1fr; gap:14px;
  padding:7px 10px 7px 12px; border-left:2px solid transparent;
  align-items:baseline;
}}
.crit--clear {{ border-left-color:var(--rule); }}
.crit--accepted {{ border-left-color:var(--tripped); background:var(--tripped-bg); }}
.crit--regression {{ border-left-color:var(--madder); background:var(--madder-bg); }}
.crit__id {{
  font-family:"IBM Plex Mono",ui-monospace,monospace; font-size:0.78rem;
  color:var(--ink-faint);
}}
.crit--accepted .crit__id {{ color:var(--tripped); font-weight:500; }}
.crit--regression .crit__id {{ color:var(--madder); font-weight:500; }}
.crit__detail {{ font-size:0.9rem; color:var(--ink-soft); }}
.crit--accepted .crit__detail, .crit--regression .crit__detail {{ color:var(--ink); }}
.crit__note {{
  display:block; margin-top:5px; font-size:0.82rem; color:var(--ink-faint);
  border-left:2px solid var(--rule); padding-left:10px;
}}
.crit__note--new {{ color:var(--madder); border-left-color:var(--madder);
  font-weight:500; }}

.foot {{
  margin-top:80px; padding-top:22px; border-top:1px solid var(--rule);
  color:var(--ink-faint); font-size:0.83rem;
  font-family:"IBM Plex Mono",ui-monospace,monospace;
}}
@media (max-width:620px) {{
  .crit {{ grid-template-columns:1fr; gap:3px; }}
  .mast {{ padding-top:36px; }}
  .tally {{ text-align:left; }}
}}
</style>

<div class="wrap">
  <header class="mast">
    <div>
      <h1 class="mast__t">Full run ledger</h1>
      <p class="mast__sub">Every stage of the design platform, run end to end in
      one command against the cached sources. {E(run.get('provenance',''))}</p>
    </div>
    <div class="tally">
      <span class="tally__n">{run.get('clear', 0)}<span style="color:var(--ink-faint)">/{run.get('total', 0)}</span></span>
      <span class="tally__l">criteria clear &middot; {run.get('minutes','?')} min</span>
    </div>
  </header>

  <section class="block" style="margin-top:0">
    <p class="eyebrow">Stage ledger</p>
    <h2>What ran, and what it cost</h2>
    <p class="lede">Each stage is verified by its own process against criteria
    written before the run. {tripped_total} criteria are tripped. {verdict}</p>
    <div class="scroll">
      <table>
        <thead><tr>
          <th>#</th><th>Stage</th><th>Criteria</th><th>State</th><th>Time</th>
        </tr></thead>
        <tbody>{''.join(ledger)}</tbody>
      </table>
    </div>
  </section>

  <section class="block">
    <p class="eyebrow">The result</p>
    <h2>Where {pool} candidate targets stop</h2>
    <p class="lede">Every candidate is attributed to the first gate it fails, so
    the drops sum to the pool rather than overlapping. Nothing reaches the end,
    and that is a measurement of the constraints — not a failure of the run.</p>
    <ol class="gates">{''.join(chain)}</ol>
  </section>

  <section class="block">
    <p class="eyebrow">Served over HTTP</p>
    <h2>An empty pipeline answers 200</h2>
    <p class="lede">The endpoints return a named status with the reasons
    computed from this run — never a 404, a 500, or a bare empty list. An empty
    list would read as having looked and found nothing to say.</p>
    <div class="endings">{''.join(endings)}</div>
  </section>

  <section class="block">
    <p class="eyebrow">Every criterion</p>
    <h2>The full ledger</h2>
    <p class="lede">All {run.get('total', 0)} criteria, in pipeline order, exactly as
    each verifier reported them.</p>
    {''.join(sections)}
  </section>

  <p class="foot">Generated from {E(run.get('report', 'the run'))} &middot;
  {'derived artifacts rebuilt from the cache' if run.get('fresh') else 'derived artifacts reused'}
  &middot; no value on this page is transcribed by hand</p>
</div>
"""


def main() -> int:
    # Taken as an argument rather than hardcoded: run_all.py accepts --report,
    # and a fixed path here would happily render a *previous* run's page while
    # printing a confident success line.
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else (
        ROOT / "reports" / "full-run.json")
    if not source.exists():
        print(f"No run at {source}. Run run_all.py first.")
        return 1

    run = json.loads(source.read_text(encoding="utf-8"))
    if not run.get("stages"):
        print("The run records no stages; refusing to write a blank page.")
        return 2
    missing = [s["number"] for s in run["stages"] if not s.get("criteria")]
    if missing:
        print(f"Stages {', '.join(missing)} record no criteria; refusing to "
              "render a page that would look complete.")
        return 2

    out = source.with_suffix(".html")
    out.write_text(render(run), encoding="utf-8")
    criteria = sum(len(s["criteria"]) for s in run["stages"])
    new = len(run.get("unexpected", []))
    stale = len(run.get("stale", []))
    print(f"{len(run['stages'])} stages, {criteria} criteria -> {out}")
    if new or stale:
        print(f"  {new} regression(s) and {stale} stale exemption(s) marked "
              "on the page")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
