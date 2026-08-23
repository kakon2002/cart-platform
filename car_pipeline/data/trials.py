"""Interventional trial counts per antigen, from the public registry.

**Matching is inexact by construction and the output says so.** The query is free
text over the registry, so a returned study is one that *mentions* the antigen,
not one testing a binder against it. A count reported as the second would be a
claim the registry cannot support — this project has already been bitten once by
a free-text search returning a plausible integer that was entirely spurious.

What the count is good for is the terminated, withdrawn and suspended tally. That
is not proof of a safety problem and is not reported as one; it is a signal that
something happened to people already and is worth reading before dosing more.

**Two limits, both measured, both carried into every row rather than assumed
away.**

*The tallies cover one page.* `total` is the registry-wide count; the stopped,
phase and CAR tallies are computed over the studies actually returned. For an
antigen with more studies than the page holds those tallies are a floor, and the
row says so with `truncated`. Paginating an antigen with 26,605 studies to count
its terminations is not proportionate to what the number is used for.

*The query is the gene symbol only.* Synonyms are not searched, and the
undercount that causes is large rather than marginal: measured live, `CLDN18`
returns 3 studies while `Claudin 18.2` returns 156. A zero here means "no study
mentions this symbol", never "this antigen is untried".
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import CacheEntry, DataSource, _write_json_atomic

BASE = "https://clinicaltrials.gov/api/v2/studies"
USER_AGENT = "car-platform/stage9"
RELEASE_PIN = "v2"
PAGE_SIZE = 200

#: Statuses that mean a trial stopped before its planned end.
STOPPED = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}


@dataclass
class TrialSummary:
    antigen: str
    total: int = 0
    #: Studies actually inspected. The tallies below cover these, not `total`.
    returned: int = 0
    truncated: bool = False
    phases: dict[str, int] = field(default_factory=dict)
    stopped: int = 0
    stopped_ids: list[str] = field(default_factory=list)
    car_mentioning: int = 0

    @property
    def has_stopped(self) -> bool:
        return self.stopped > 0


class TrialSource(DataSource):
    name = "ClinicalTrials.gov"
    namespace = "trials"

    def __init__(self, antigens: list[str] | None = None, **kwargs):
        super().__init__(**kwargs)
        self.antigens = sorted(antigens or [])

    def cache_entries(self) -> Iterable[CacheEntry]:
        return [
            CacheEntry(
                key="counts",
                filename="trial_counts.json",
                fingerprint={
                    "release": RELEASE_PIN,
                    "antigens": self.antigens,
                    "measure": "interventional_counts",
                },
            )
        ]

    def _query(self, antigen: str) -> dict:
        params = {
            "query.term": antigen,
            "filter.overallStatus": "",
            "pageSize": PAGE_SIZE,
            "countTotal": "true",
            "fields": "NCTId|OverallStatus|Phase|BriefTitle",
        }
        params = {k: v for k, v in params.items() if v != ""}
        url = f"{BASE}?{urllib.parse.urlencode(params)}"
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=90) as response:
            return json.loads(response.read())

    def fetch(self) -> Path:
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            print(f"  querying the trial registry for {len(self.antigens)} antigens",
                  flush=True)
            payload = {}
            for n, antigen in enumerate(self.antigens, 1):
                data = self._query(antigen)
                studies = data.get("studies", [])
                phases: dict[str, int] = {}
                stopped, stopped_ids, car = 0, [], 0
                for study in studies:
                    protocol = study.get("protocolSection", {})
                    ident = protocol.get("identificationModule", {})
                    status = protocol.get("statusModule", {}).get("overallStatus", "")
                    for phase in protocol.get("designModule", {}).get("phases", []):
                        phases[phase] = phases.get(phase, 0) + 1
                    if status in STOPPED:
                        stopped += 1
                        stopped_ids.append(ident.get("nctId", ""))
                    title = (ident.get("briefTitle", "") or "").upper()
                    if "CAR" in title.split() or "CAR-T" in title:
                        car += 1
                total = data.get("totalCount", len(studies))
                payload[antigen] = {
                    "total": total,
                    "returned": len(studies),
                    # The tallies below cover `returned`, not `total`.
                    "truncated": total > len(studies),
                    "phases": phases,
                    "stopped": stopped,
                    "stopped_ids": stopped_ids[:10],
                    "car_mentioning": car,
                }
                if n % 25 == 0:
                    print(f"    {n}/{len(self.antigens)}", flush=True)
            _write_json_atomic(tmp, payload)
            import hashlib

            blob = json.dumps(payload, sort_keys=True).encode("utf-8")
            return {
                "digest": hashlib.sha256(blob).hexdigest(),
                "declared_rows": len(self.antigens),
                "observed_rows": len(payload),
                "extra": {"source": BASE},
            }

        return self.cache.ensure(entry, fetcher)

    def load(self) -> dict[str, TrialSummary]:
        entry = next(iter(self.cache_entries()))
        if not self.cache.is_valid(entry):
            self.fetch()
        raw = json.loads(self.cache.path(entry).read_text(encoding="utf-8"))
        return {
            antigen: TrialSummary(
                antigen=antigen,
                total=row["total"],
                returned=row.get("returned", 0),
                truncated=bool(row.get("truncated")),
                phases=row["phases"],
                stopped=row["stopped"],
                stopped_ids=row["stopped_ids"],
                car_mentioning=row["car_mentioning"],
            )
            for antigen, row in raw.items()
        }
