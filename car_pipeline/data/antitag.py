"""The adaptor receptor's anti-tag binder, retrieved from a deposited structure."""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterable

from car_pipeline.data.source import CacheEntry, DataSource, _write_json_atomic

ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry"
ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity"
USER_AGENT = "car-platform/antitag"
TIMEOUT = 60
RETRIES = 3
RETRY_BACKOFF = 2.0

ENTRY_ID = "1P4B"
REVISION_PIN = "1.4"
BINDER_ENTITIES = ("1", "2")
ANTIGEN_ENTITY = "3"

TAG_SYSTEM = "peptide neo-epitope, GCN4(7P-14P)"


class AntiTagError(RuntimeError):
    """The binder could not be retrieved under its pinned terms."""


def _get(url: str) -> dict:
    """One JSON GET, retried on transport failure and never invented."""
    last = None
    for attempt in range(RETRIES):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
                return json.loads(response.read())
        except urllib.error.HTTPError as exc:
            raise AntiTagError(f"{url} answered {exc.code}") from exc
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError) as exc:
            last = exc
            if attempt < RETRIES - 1:
                time.sleep(RETRY_BACKOFF * (attempt + 1))
    raise AntiTagError(f"{url} unreachable after {RETRIES} attempts: {last}")


class AntiTagSource(DataSource):
    name = "RCSB"
    namespace = "antitag"

    def cache_entries(self) -> Iterable[CacheEntry]:
        return [
            CacheEntry(
                key=f"binder__{ENTRY_ID}",
                filename=f"binder__{ENTRY_ID}.json",
                fingerprint={
                    "entry": ENTRY_ID,
                    "entities": list(BINDER_ENTITIES),
                    "revision": REVISION_PIN,
                    "measure": "anti_tag_binder",
                },
            )
        ]

    def fetch(self) -> Path:
        """Retrieve the pinned entity sequences, refusing a moved revision."""
        entry = next(iter(self.cache_entries()))

        def fetcher(tmp: Path) -> dict:
            print(f"  retrieving anti-tag binder from {ENTRY_ID}", flush=True)
            header = _get(f"{ENTRY_URL}/{ENTRY_ID}")
            accession = header.get("rcsb_accession_info") or {}
            revision = "%s.%s" % (accession.get("major_revision"),
                                  accession.get("minor_revision"))
            status = accession.get("status_code")
            if status != "REL":
                raise AntiTagError(
                    f"{ENTRY_ID} is {status}, not REL; a withdrawn or superseded "
                    "entry is not a source"
                )
            if revision != REVISION_PIN:
                raise AntiTagError(
                    f"{ENTRY_ID} is at revision {revision}, pinned at "
                    f"{REVISION_PIN}. The deposited sequence may have changed; "
                    "re-pin deliberately rather than accepting a moved entry"
                )
            declared = (header.get("rcsb_entry_container_identifiers") or {}).get(
                "polymer_entity_ids") or []
            missing = [e for e in BINDER_ENTITIES if e not in declared]
            if missing:
                raise AntiTagError(
                    f"{ENTRY_ID} does not declare entity/entities {missing}; "
                    f"it declares {declared}"
                )

            chains = []
            for entity_id in BINDER_ENTITIES:
                body = _get(f"{ENTITY_URL}/{ENTRY_ID}/{entity_id}")
                poly = body.get("entity_poly") or {}
                sequence = poly.get("pdbx_seq_one_letter_code_can")
                if not sequence:
                    raise AntiTagError(
                        f"{ENTRY_ID} entity {entity_id} carries no canonical "
                        "one-letter sequence"
                    )
                declared_length = poly.get("rcsb_sample_sequence_length")
                if declared_length is not None and declared_length != len(sequence):
                    raise AntiTagError(
                        f"{ENTRY_ID} entity {entity_id} declares "
                        f"{declared_length} residues and carries {len(sequence)}"
                    )
                described = (body.get("rcsb_polymer_entity") or {}).get(
                    "pdbx_description", "")
                chains.append({
                    "entity": entity_id,
                    "description": described,
                    "sequence": sequence,
                    "residues": len(sequence),
                })

            payload = {
                "entry": ENTRY_ID,
                "revision": revision,
                "status": status,
                "tag_system": TAG_SYSTEM,
                "antigen_entity_excluded": ANTIGEN_ENTITY,
                "chains": chains,
                "residues": sum(c["residues"] for c in chains),
            }
            _write_json_atomic(tmp, payload)
            digest = hashlib.sha256(tmp.read_bytes()).hexdigest()
            return {
                "declared_rows": len(BINDER_ENTITIES),
                "observed_rows": len(chains),
                "digest": digest,
                "extra": {
                    "entry": ENTRY_ID,
                    "revision": revision,
                    "tag_system": TAG_SYSTEM,
                    "residues": payload["residues"],
                    "antigen_entity_excluded": ANTIGEN_ENTITY,
                },
            }

        return self.cache.ensure(entry, fetcher)

    def load(self) -> dict:
        """The retrieved binder, or a refusal naming what is absent."""
        entry = next(iter(self.cache_entries()))
        path = self.cache.path(entry)
        if not self.cache.is_valid(entry):
            raise AntiTagError(
                f"no cached anti-tag binder at {path}; run the source fetch"
            )
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)

    def sequence(self) -> str:
        """The binder as deposited: its entities concatenated, nothing edited."""
        payload = self.load()
        return "".join(c["sequence"] for c in payload["chains"])
