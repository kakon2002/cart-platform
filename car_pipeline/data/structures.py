"""Deposited structures, retrieved by accession rather than by name.

Retrieval is anchored on the UniProt accession, never on a full-text search of
the symbol. Measured while writing the specification: a full-text query for one
pool target returned 369 entries of which the top hits were a bacterial RNA
chaperone, a sulfur transferase and a photosystem supercomplex — a 200 response
carrying valid JSON and a plausible integer that was entirely spurious. Name
matching is how this project has produced its worst answers and it is not used
here.

**The zero contract.** A query with no hits answers **HTTP 204 with an empty
body**, not 200 with a count of zero. Parsing that body raises, and the obvious
repair — catching the exception and returning "no structures" — makes a broken
query indistinguishable from an honest absence. 204 with zero bytes is the only
accepted representation of zero; anything else is an error and is raised.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from typing import Iterable

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/rest/v1/core/entry"
ACCESSION_ATTRIBUTE = (
    "rcsb_polymer_entity_container_identifiers"
    ".reference_sequence_identifiers.database_accession"
)
DATABASE_ATTRIBUTE = (
    "rcsb_polymer_entity_container_identifiers"
    ".reference_sequence_identifiers.database_name"
)
USER_AGENT = "car-platform/stage5"
#: Page size. Exceeding it raises rather than truncating.
PAGE_ROWS = 500
TIMEOUT = 60
#: Transport failures are retried; an answer is never invented from one.
RETRIES = 3
RETRY_BACKOFF = 2.0


class RetrievalError(RuntimeError):
    """A response that is neither a parsed result nor an honest zero."""


def _query_body(accession: str) -> dict:
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": ACCESSION_ATTRIBUTE,
                        "operator": "exact_match",
                        "value": accession,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": DATABASE_ATTRIBUTE,
                        "operator": "exact_match",
                        "value": "UniProt",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": PAGE_ROWS}},
    }


def entries_for(accession: str) -> list[str]:
    """Entry identifiers whose polymer entities cross-reference this accession.

    Returns an empty list only for a 204 with a zero-length body. Every other
    shape raises: a spurious empty list here would be recorded downstream as
    "the literature holds nothing for this target", which is a claim about
    biology made out of a transport failure.
    """
    body = json.dumps(_query_body(accession)).encode("utf-8")
    request = urllib.request.Request(
        SEARCH_URL,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": USER_AGENT},
    )
    # A dropped connection is a transport failure, not an answer. Retried a
    # fixed number of times, then raised — never converted into an empty list,
    # which downstream would record as "the literature holds nothing".
    last: Exception | None = None
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                status = response.status
                payload = response.read()
            break
        except urllib.error.HTTPError as exc:
            if exc.code == 204:
                return []
            raise RetrievalError(
                f"{accession}: search returned HTTP {exc.code}"
            ) from exc
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError) as exc:
            last = exc
            if attempt + 1 < RETRIES:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
    else:
        raise RetrievalError(
            f"{accession}: {RETRIES} attempts failed, last {type(last).__name__}"
        ) from last

    if status == 204:
        if payload:
            raise RetrievalError(
                f"{accession}: HTTP 204 carried {len(payload)} bytes; a no-hit "
                "response must be empty"
            )
        return []
    if not payload:
        raise RetrievalError(
            f"{accession}: HTTP {status} with an empty body; only 204 may be empty"
        )
    parsed = json.loads(payload)
    rows = [row["identifier"] for row in parsed.get("result_set", [])]
    total = parsed.get("total_count")
    # A page cap that quietly drops results would undercount both the entries and
    # the candidates drawn from them, and would look like a protein with fewer
    # structures rather than like a truncated read.
    if total is not None and total > len(rows):
        raise RetrievalError(
            f"{accession}: {total} entries but only {len(rows)} returned; the "
            f"page size of {PAGE_ROWS} truncated the result"
        )
    return rows


def entry_summary(entry_id: str) -> dict:
    """Title and experimental method for one entry."""
    request = urllib.request.Request(
        f"{DATA_URL}/{entry_id}", headers={"User-Agent": USER_AGENT}
    )
    with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
        parsed = json.loads(response.read())
    methods = [m.get("method", "") for m in parsed.get("exptl", [])]
    return {
        "id": entry_id,
        "title": parsed.get("struct", {}).get("title", ""),
        "methods": methods,
        # Recorded so a computed model can never be reported as retrieved
        # evidence. §1 forbids presenting a prediction beside a measurement.
        "is_model": any("THEORETICAL" in m.upper() for m in methods),
    }


def entries_for_all(accessions: Iterable[str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for accession in accessions:
        out[accession] = entries_for(accession)
    return out
