# Design decisions

Every entry here is a place where the obvious implementation is wrong. Most were
written after the obvious version had already been built and had already produced
a plausible answer that was not true.

The source carries none of this any more. Each function keeps a one-line docstring
and nothing else, so this file is the only record of *why* the code is shaped the
way it is. Read the entry before changing the line it names.

Line numbers were taken before the strip and will have moved; the surrounding
names have not.

## The two rules underneath most of this

**Missing is a third state.** An absent measurement is never zero, never the
midpoint, never the mean, and never quietly dropped. Almost every status this
codebase defines exists to keep an unmeasured thing from being scored as a
measured one.

**A bound is not a measurement.** Where a number is the best the evidence
supports rather than the thing itself, the code says so and carries the flag
forward, so no reader downstream can mistake one for the other.

A third rule governs the criteria specifically: **every check needs a positive
pin.** A set of purely negative assertions passes completely against a dead code
path. That is not hypothetical here - it is how a stage returned nothing for all
200 targets while every check on it passed.

492 decisions across 48 modules.

## Contents

- **The contract** (9) - `spec.py`, `project.py`
- **Indications** (22) - `indication.py`, `pdac.py`, `breast.py`, `registry.py`
- **Sources and the cache** (140) - `source.py`, `uniprot.py`, `hpa.py`, `gtex.py`, `tcga.py`, `singlecell.py`, `depmap.py`, `genespan.py`, `antibodies.py`, `structures.py`, `domains.py`, `trials.py`, `coverage.py`, `availability.py`
- **Stages** (185) - `stage1.py`, `stage3.py`, `stage4.py`, `routing.py`, `stage5.py`, `stage6.py`, `stage9.py`, `stage10.py`, `stage11.py`, `stage12.py`, `validation.py`, `construct_safety.py`
- **The service** (41) - `pipeline.py`, `server.py`
- **Running it** (62) - `run_all.py`, `bootstrap.py`, `make_artifact.py`, `make_brief.py`, `make_package.py`, `strip_comments.py`
- **The criteria** (79) - `verify_schema.py`, `verify_surface.py`, `verify_ranking.py`, `verify_ranking_final.py`, `verify_pairing.py`, `verify_routing.py`, `verify_binders.py`, `verify_construct.py`, `verify_safety.py`, `verify_developability.py`, `verify_package.py`, `verify_api.py`, `verify_indications.py`


---

# The contract

## `car_pipeline/schemas/spec.py`

**1. spec.py:14** - DatasetStatus has three states rather than the obvious available/unavailable pair: not_configured (no connector exists for this source at all) is kept distinct from unreachable (the connector exists but the data is not readable).

> *Without it:* Collapsing them hides whether the gap is unbuilt or merely unfetched — an operator cannot tell a source that needs a connector written from one that just needs a refetch, and both report identically as missing data.

**2. spec.py:44** - Terminable_risk_ceiling is a second, independent ceiling rather than an adjustment or discount applied to normal_tissue_risk_ceiling for adaptor designs.

> *Without it:* An adaptor design does not make an antigen safer — the adaptor still binds it — it makes the exposure terminable, because activation needs a separately dosed protein. Magnitude and reversibility are different axes, so folding reversibility into the magnitude ceiling scores a design as low-risk when the risk is unchanged and merely stoppable.

**3. spec.py:50** - Terminable_risk_ceiling on DesignConstraints is Optional and defaults to None instead of carrying a platform default like the other constraints.

> *Without it:* Both ceilings are policy inputs, not measurements; a default here would be this code quietly setting clinical policy. Its absence disables the adaptor row instead of admitting adaptor designs against a threshold no clinician chose.

**4. spec.py:70** - Data_availability_score counts only required=True datasets, excluding the optional ones (IEDB immunogenicity, ClinicalTrials.gov prior outcomes) from both numerator and denominator.

> *Without it:* Including optional sources caps an otherwise complete run below 1.0, and it dilutes the denominator so a real blocking gap looks milder than it is — one missing required source among many rows reads as a small shortfall rather than a stop.


## `car_pipeline/schemas/project.py`

**1. project.py:3** - Every model sets ConfigDict(extra="forbid") so unknown fields raise, rather than accepting and ignoring them as pydantic does by default.

> *Without it:* A mistyped field name has to fail loudly: silently accepting `target_antigens` leaves the run in discovery mode while the config appears to carry a supplied target, and nothing downstream can tell the difference between that and a genuine screen.

**2. project.py:117** - A blank or whitespace-only target_antigen is normalised to None by a field validator instead of being kept as the string the form submitted.

> *Without it:* An empty string is truthy as a supplied value, so left as-is it counts as a target: discovery_mode returns VALIDATE and the entire antigen screen is skipped for a project that named no antigen.

**3. project.py:129** - Discovery_mode is a computed property derived from whether target_antigen is None, never a stored field on the input.

> *Without it:* A stored field lets the two disagree — a config can declare mode B while carrying a target, or mode A with none — and there is no way to tell which of the two the run should honour.

**4. project.py:95** - Terminable_risk_ceiling is declared per project and defaults to None; it is never derived from safety_tolerance the way the normal-tissue ceiling is.

> *Without it:* Deriving or defaulting it would have the platform inventing a clinical risk tolerance it has no basis for. Its absence is meaningful — it disables the adaptor row rather than silently authorising an exposure level nobody declared.

**5. project.py:72** - A tissue criticality override requires a rationale of at least 20 characters, which travels into the output header, rather than accepting a bare tier number.

> *Without it:* Overrides move the platform's safety defaults for a tissue. Without the mandatory rationale a relaxed criticality tier is invisible in the output: a reader sees the resulting risk numbers with no record of which safety defaults were moved or why.


---

# Indications

## `car_pipeline/configs/indication.py`

**1. indication.py:24** - Every per-atlas column name is declared in AtlasSchema rather than inferred from the accession.

> *Without it:* None of it is derivable: one GEO submitter chose a name and a later submitter chose a different one. Guessing a column name and getting silence is how an atlas gets read with the wrong annotations.

**2. indication.py:44** - `counts_path` is declared per atlas (reference submission `layers/counts`, CELLxGENE export `raw/X`) instead of assuming one location.

> *Without it:* Reading the wrong one is not an error that surfaces — it silently scores against normalised values as though they were raw counts.

**3. indication.py:49** - `symbol_field` and `ensembl_field` are both declared, because the two submissions have them the opposite way round (reference indexes by symbol with `ensg` beside it; the CELLxGENE export indexes by Ensembl with `feature_name` beside it).

> *Without it:* Assuming either layout reads identifiers out of the wrong var field for the other atlas.

**4. indication.py:41** - Coarse labels not in `compartment_map` become "other" and are reported, rather than dropped.

> *Without it:* Silently dropping unmapped cells shrinks the stromal/immune population a tumour antigen has to beat, without anything saying so.

**5. indication.py:98** - The GEO supplementary bucket is derived from the accession digits instead of carried as a literal.

> *Without it:* The previous literal "GSE202nnn" produced a wrong URL for every accession outside that block — a 404 that looks like an unreachable source rather than a wrong address.

**6. indication.py:76** - A None or empty DepMap lineage makes the dependency component unavailable, not zero.

> *Without it:* Scoring an unmeasured dependency as 0.0 puts a number in the composite that was never measured, and it is indistinguishable from a real measurement of no dependency.

**7. indication.py:81** - `atlas=None` marks the whole screen NOT_USABLE rather than reducing the available weight.

> *Without it:* It is not a missing 0.45 of weight, it is the loss of the only component that rejects stromal and immune genes.

**8. indication.py:55** - `treatment_column`/`untreated_label` are optional; an atlas without a treatment split simply has no untreated subset.

> *Without it:* Fabricating an untreated subset (e.g. treating all cells as untreated) would report a treatment-naive measurement that was never made.


## `car_pipeline/configs/pdac.py`

**1. pdac.py:3** - `target_antigen` is deliberately absent; the null is what selects discovery mode.

> *Without it:* Filling it in turns the screen into a validation of something already assumed — the platform would confirm the chosen antigen instead of searching for one.

**2. pdac.py:42** - Pancreas is overridden to criticality tier 2 for this indication only, with the rationale recorded, while the platform default holds pancreas at tier 1 everywhere else.

> *Without it:* The tumour arises in this organ and the intended population is surgically resected, so expression in normal pancreas is a weaker objection here than the same expression in lung or heart. Applying the relaxation platform-wide would clear pancreas-expressing targets for indications where the organ is still in place.

**3. pdac.py:56** - `terminable_risk_ceiling=0.35` is a policy input fixed before the run and pinned in the spec, distinct from the 0.15 accepted from a T cell that cannot be withdrawn; the pipeline does not derive it.

> *Without it:* This pipeline cannot measure the risk a project should accept from a stoppable exposure. Criterion A9 therefore reports the admitted count across the whole sweep so the choice is visible, and A10 trips if this value ever stops matching the one recorded in the spec — otherwise the ceiling could be tuned to admit more targets with nothing noticing.

**4. pdac.py:86** - The normal-tissue denominator is the bulk GTEx "Pancreas" column, chosen by hand.

> *Without it:* Four GTEx pancreas columns exist and three are cell-sorted fractions; picking one of those as the denominator measures the tumour-versus-normal margin against a sorted cell population rather than the organ. It was a judgement call recorded in stage 3, not an obvious lookup.


## `car_pipeline/configs/breast.py`

**1. breast.py:10** - Breast was taken as the second indication over colorectal and lung because its atlas is a structural drop-in for the existing loader.

> *Without it:* The matrix is CSR in the same layout, the normalisation matches the reference submission (per-cell sum of expm1 ~9,974 against ~10,000), the variable index is 100% Ensembl, and malignancy sits at the coarse annotation level, which is the shape the malignant-cell reader requires. It also brings 113 solid-tissue normal samples against pancreas's 4 — the tumour-versus-normal margin is 28x better supported here than in the indication the platform was designed on.

**2. breast.py:32** - The atlas is the CELLxGENE standardised export, not the GEO supplementary files the reference indication uses.

> *Without it:* The export is a single h5ad with the annotations already in `obs`; the GEO deposit is a matrix plus metadata that would have to be reassembled. 843,892,052 bytes, verified reachable with HDF5 magic at byte 0.

**3. breast.py:38** - The URL is pinned to a specific CELLxGENE asset UUID, and a 404 is treated as a signal rather than a transient failure.

> *Without it:* CELLxGENE re-issues asset UUIDs when a dataset version changes, so this pin is a version pin: a fetch that 404s means the dataset moved, and retrying it as a network blip would silently keep an unpinned, different version of the data in play.

**4. breast.py:46** - Malignancy is read from the coarse annotation `celltype_major`, corroborated by an independent cross-tabulation rather than taken on trust.

> *Without it:* The malignant reader requires malignancy at the coarse level. Cross-tabulating celltype_major against normal_cell_call gives (Cancer Epithelial, cancer) = 24,489 exactly and (Normal Epithelial, normal) = 4,355, which is what confirms the label is the malignant call and not a near-synonym.

**5. breast.py:53** - PVL (perivascular-like) maps to "fibroblast" rather than falling through to "other".

> *Without it:* It is a mural stromal population; sending it to "other" excludes 5,423 stromal cells from the comparison a tumour antigen has to win, inflating malignant_vs_stroma for anything those cells express.

**6. breast.py:69** - Counts_path is "raw/X" and the symbol/Ensembl fields are inverted from the reference submission (symbol in `feature_name`, Ensembl in `_index`).

> *Without it:* Carrying the reference defaults over would read normalised values as counts and index genes by the wrong var field, with no error either time.

**7. breast.py:75** - `untreated_label` is "Naïve" with the diaeresis, verified against the file rather than typed from the paper.

> *Without it:* "Naive" would have matched nothing, and the untreated subset would have come back empty with no error.

**8. breast.py:92** - `tissue_overrides` is deliberately left empty for breast.

> *Without it:* Breast is tier 3 by platform default and nothing about this indication argues for relaxing it. Pancreas needed a tier-2 override only because it defaults to tier 1; inventing an override here to make more targets clear is the move this project exists to refuse.


## `car_pipeline/configs/registry.py`

**1. registry.py:16** - Indications are keyed by the lowercased cancer type as a caller would type it, not by the internal indication key.

> *Without it:* The key a caller sends is free text from a form. Keying on the internal identifier makes every request fail unless the caller already knows the platform's private vocabulary.

**2. registry.py:27** - Short aliases resolve alongside the full names.

> *Without it:* Without them a caller has to reproduce the full oncological name exactly to get a run, and a near-miss reads as an unsupported indication rather than a spelling difference.


---

# Sources and the cache

## `car_pipeline/data/source.py`

**1. source.py:132** - The payload is fsynced, moved into place, and only then is the manifest written — the manifest, not the file, is the commit marker.

> *Without it:* Writing the manifest first (or treating the file's existence as completion) means a run interrupted mid-download reads as complete forever; a truncated file that looks finished is undetectable later, whereas a missing manifest fails in the safe direction.

**2. source.py:63** - _atomic_replace opens the temp file read-write purely to fsync the descriptor before os.replace, rather than renaming straight away.

> *Without it:* Without the flush through to the device, the manifest that blesses the payload can become durable before the payload itself; after a crash the manifest points at bytes that were never written, and every later run trusts it. The handle must be writable for fsync to be permitted, which is why it is reopened O_RDWR.

**3. source.py:55** - Cache validity is keyed on a fingerprint of the whole request — query terms, field lists, release pins — not on the filename or URL.

> *Without it:* Keying on the file alone silently serves data fetched under different terms: change the query or the release pin and the stale payload is returned as if it answered the new request. is_valid also requires a matching manifest_version and an exact byte-size match, so a file edited or truncated after the fact fails instead of passing.

**4. source.py:133** - Declared row counts are compared against observed counts at fetch time, inside commit, and a mismatch deletes the temp file and raises IntegrityError instead of caching it.

> *Without it:* A dropped page in a paginated download is byte-for-byte indistinguishable from a complete one once it is on disk, so no later check can recover it — the source just quietly has fewer rows and every downstream count is plausible and wrong.

**5. source.py:140** - When a source supplies no declared row count, commit prints an explicit 'completeness is unverified' note and records count_verified: false in the manifest, rather than passing over it in silence.

> *Without it:* An entry nobody could check reads exactly like one that passed its check. Silently accepting it means an unverifiable download and a verified one are indistinguishable in the cache afterwards.

**6. source.py:271** - Content-Length and the source-published MD5 are verified inside stream_to_file, before anything is moved into place, and the destination is unlinked on mismatch.

> *Without it:* Verifying after the move (or not at all) commits a truncated or corrupted payload and then blesses it with a manifest; the size check in is_valid would agree with it forever because the manifest records the short length as the truth.

**7. source.py:309** - Capture_headers records the release/version headers the server actually returned on the first page into the manifest, instead of recording only the release the caller pinned.

> *Without it:* Otherwise the release label is an assertion the payload never has to honour — a source that quietly serves a different release than the pin is cached under the pinned name and every span, TPM, or ID derived from it is attributed to the wrong release.

**8. source.py:351** - The pagination 'next' link is extracted with a regex on <target>; rel="next" rather than by splitting the Link header on commas.

> *Without it:* The bracketed target can itself contain commas (field lists, multi-value query params), so comma-splitting chops the URL mid-target: pagination stops early or follows a truncated link, and the run silently ends up with a subset of rows.

**9. source.py:337** - In the paginated writer the header row is written only from page 0 and stripped from every subsequent page, and the row counter counts only post-strip lines.

> *Without it:* Naively concatenating pages interleaves a header line into the middle of the data as if it were a record, and inflates the observed count by one per page, so the observed-versus-declared integrity check passes on a file that is wrong in both content and length.

**10. source.py:180** - Fetchers only write into the temporary path they are handed and return manifest metadata; they are never allowed to move anything into place — commit is the sole writer of the canonical path.

> *Without it:* If a fetcher wrote to the destination directly, an interrupted or failed fetch leaves partial bytes at the canonical path while the previous manifest is still on disk, and the integrity checks that would have rejected the payload never get to run.

**11. source.py:460** - Is_cached requires bool(entries) AND every declared entry valid — the non-empty guard is deliberate, not defensive noise.

> *Without it:* All() over an empty list is True, so a source that declares no cache entries (a stub, or a connector whose entry list failed to build) would report itself fully cached and available. A partially cached source likewise must not read as usable — one missing entry means the stage fails mid-run instead of up front.

**12. source.py:370** - Gzip expansion and zip member extraction run in fixed 1 MiB blocks (CHUNK = 1 << 20) with neither side ever fully resident.

> *Without it:* The obvious read-it-all form dies on the archives this pipeline actually pulls — the coverage report cites an eight gigabyte archive — so the whole-file read is not a slower path, it is a path that does not complete.


## `car_pipeline/data/uniprot.py`

**1. uniprot.py:3** - Requests a seventh field, lipidation (ft_lipid), instead of deciding membrane attachment from transmembrane segments and the localisation text alone.

> *Without it:* GPI-anchored proteins leave no transmembrane segment behind and are not consistently spelled out in the localisation text, so without ft_lipid that entire class of protein is invisible to the attachment gate — silently dropping several of the best known targets.

**2. uniprot.py:9** - Requests an eighth field, chain boundaries (ft_chain), even though the surface filter never consults it and adding it changed nothing about the surface set.

> *Without it:* A mature protein is not always one molecule: where a precursor is cleaved, some resulting chains are released rather than held at the surface, so a binder raised against a released chain meets its antigen in plasma rather than on a cell. The field is not free — the field list sits inside the cache fingerprint, so adding it invalidated the proteome cache and re-fetched all 20,431 entries — but it lets a later stage tell the anchored chain from the shed one instead of treating the precursor as uniformly reachable.

**3. uniprot.py:48** - Pins the UniProt release by validating the release the service reports on the response, rather than by asking for a release in the request.

> *Without it:* The search service always serves its current release and offers no way to request an older one, so a request-side pin is unenforceable. The manifest records the release only as a label, so a service that had moved on would file a different proteome under the pinned name 2026_02 and every count measured against the old release would drift with nothing raising. Bumping RELEASE_PIN also changes the cache fingerprint, so a deliberate bump invalidates and re-fetches instead of silently replacing cache contents.

**4. uniprot.py:376** - Treats a missing X-UniProt-Release header as a hard CacheError instead of proceeding when the service simply did not say which release it served.

> *Without it:* An unlabelled fetch would be cached under the pinned name 2026_02 on no evidence at all, reintroducing exactly the silent substitution the pin exists to catch. check_release is kept separate from fetch so this guard can be exercised without downloading 20,431 entries.

**5. uniprot.py:67** - Matches the plasma-membrane phrases case-insensitively as substrings of the localisation text rather than as exact location values.

> *Without it:* Polarised variants such as "Apical cell membrane" would not match an exact test, silently withholding genuine cell-surface proteins from the admitted set.

**6. uniprot.py:72** - Admission reads only location statements — free-text notes are stripped out by location_statements and decide nothing, in either direction.

> *Without it:* Across all fourteen entries where a plasma-membrane phrase appears ONLY in a note, the note means four different things a substring test reads identically: that the protein is there ("Located on cell surface microvilli."), that it is NOT there ("Integral membrane protein not detected at the cell membrane."), that it merely passes through ("Cycles via the cell surface and endosomes upon lumenal pH disruption."), or something about lipids ("Preferentially binds to cardiolipin relative to other common cell membrane lipids."). Matching notes admits the negation; dropping notes discards the assertion — neither direction is safe.

**7. uniprot.py:177** - An entry whose only plasma-membrane evidence sits in a note is recorded in a third state (outward_note_only) — not admitted, but enumerated by name in the output rather than dropped quietly.

> *Without it:* Resolving these either way is wrong: admitting them risks an unreachable target, which is the dangerous direction, while dropping them quietly hides genuinely ambiguous annotation. The set is small enough to audit — 14 entries in the pinned release — so it is reported instead of resolved.

**8. uniprot.py:105** - The outward-facing topological note "Extracellular" is matched by exact equality on the note value, not by a substring test.

> *Without it:* "Lumenal" and "Perinuclear space" sit on the same side of the bilayer but inside an organelle; a looser test admits them as outward-facing, putting proteins that are unreachable from outside the cell into the surface set.

**9. uniprot.py:342** - Gate 2 requires positive evidence of outward facing rather than excluding a list of intracellular compartments; the compartment list is used only for reporting, never for admission.

> *Without it:* An exclusion-list gate admits the multi-pass proteins of internal compartments, which are topologically outward-facing yet unreachable from outside the cell.

**10. uniprot.py:110** - Withheld anchored proteins are split into "placed in a named compartment" (internal_anchored) and "nothing on record says where" (compartment_unresolved) and reported apart, instead of pooled as one rejected group.

> *Without it:* Pooling counts the second group as judged-and-rejected when it is only held out for want of evidence — 534 unresolved entries in the pinned release would be silently scored against, as opposed to the 1,362 that were actually measured and placed internally.

**11. uniprot.py:213** - _count_extracellular_residues returns None when no segment length could be read, instead of returning the running total of zero.

> *Without it:* A protein nobody annotated would be scored as measured and found tiny — exactly the imputation the accessibility component is required to avoid. Zero and unmeasured must stay distinguishable.

**12. uniprot.py:205** - A TOPO_DOM segment with uncertain bounds ("<", ">", "?") contributes nothing to the residue total and does not set the measured flag, rather than being counted after stripping the markers.

> *Without it:* The segment exists but its length does not; counting it invents a length, and letting it set measured on its own would turn an entry with no readable segment into a reported measurement.

**13. uniprot.py:225** - _bound returns the position together with an uncertainty flag instead of stripping the hedge markers and returning the bare number.

> *Without it:* "<37" means "somewhere at or before 37", not 37; returning 37 turns a hedge into a measurement that no caller can tell apart from an exact one. The number is kept as the best estimate and the flag beside it lets Chain.exact refuse a rule that needs a boundary that was never exact.

**14. uniprot.py:264** - Chain.contains returns False when either bound is missing, and Chain.length returns None, rather than substituting a default for the missing bound.

> *Without it:* Defaulting a missing start or end silently answers a containment question the annotation never supported, placing a residue inside or outside a chain on a guessed boundary.

**15. uniprot.py:277** - Parse_chains keeps a chain whose bounds could not be read, with empty bounds, instead of skipping the malformed row.

> *Without it:* Dropping it makes a cleaved precursor look like an uncleaved one — the difference between a protein held at the surface and one released into plasma — and understates how many pieces a precursor is cut into.

**16. uniprot.py:165** - An empty chains list is defined to mean the entry carries no chain annotation at all, and is deliberately not equated with the protein being a single chain.

> *Without it:* Reading "no annotation" as "one chain" asserts an uncleaved precursor on absent evidence, so a protein that is actually cut and partly shed would be treated as uniformly surface-anchored.

**17. uniprot.py:285** - Chain retains the annotation's chain identifier alongside the coordinates instead of storing start and end only.

> *Without it:* The stage that picks between chains has to be able to say which one it picked; a start and end alone name a range rather than a chain, so the choice becomes untraceable.

**18. uniprot.py:314** - Parse_row pads and truncates each row to len(FIELDS) rather than to a fixed hard-coded width.

> *Without it:* A row normalised to the old width drops the newly added column silently, and the drop surfaces only as an empty chain list — indistinguishable from a protein that genuinely has no chain annotation, and only after the full re-fetch of 20,431 entries had already been paid for.

**19. uniprot.py:494** - The self-check expectations are this code's own measured output against the pinned release, not figures carried forward from a prior run.

> *Without it:* The previous expected figures were an estimate that had never been this code's output; they sat 0.46% above the surface count and were read as a real discrepancy. The current values (entries 20,431; surface 3,466; single_pass 1,446; multi_pass 1,884; gpi_anchored 136; internal_anchored 1,362; compartment_unresolved 534; outward_note_only 14) are reproducible now that the release is enforced on the response, so a difference means something actually changed.


## `car_pipeline/data/hpa.py`

**1. hpa.py:3** - The atlas is pinned to release 23 rather than tracking the latest release.

> *Without it:* Per-tissue staining and pathology tables were withdrawn from the per-dataset downloads after v23; later releases ship only the consolidated gene table with no per-tissue calls at all, so a "newer" pin silently loses every per-tissue staining level. All four files come from the same release so the gene set and tissue vocabulary stay consistent with each other.

**2. hpa.py:12** - Placeholder level rows ("N/A" and empty) are dropped at parse instead of being counted as staining data, but the gene's symbol, accession and localisation call are still kept.

> *Without it:* Counting placeholder rows as coverage inflates the surface-protein staining figure by roughly a fifth (expected 1,945 of 3,402 matched surface proteins with staining). Dropping the whole gene instead would be the opposite error — absence of staining is not absence of the gene.

**3. hpa.py:46** - Gradient calls ("Ascending", "Descending") are given rank 1, the same as "Low", rather than being sorted off the end of the ordinal scale.

> *Without it:* Ranked below "Not detected" (0), a gut-enterocyte gradient call would read as the protein being absent when the call actually reports it is present — understating expression, the dangerous direction for a dataset whose whole job is flagging on-target off-tumour risk.

**4. hpa.py:60** - "Not representative" rows are dropped as withdrawn rather than mapped onto the scale.

> *Without it:* It is a withdrawn call: unscored, and not a detection. Scoring it as 0 would assert the protein was looked for and not found in that tissue; scoring it as positive would assert a detection the atlas has retracted.

**5. hpa.py:215** - A level string not in LEVEL_RANK is dropped and counted as withdrawn, instead of defaulting to rank 0 via `LEVEL_RANK.get(text, 0)`.

> *Without it:* A default of 0 turns any unrecognised vocabulary — a new level term in a future file, a stray whitespace variant — into a positive assertion of "Not detected" in that tissue, which is the finding that clears a target rather than flags it.

**6. hpa.py:89** - `peak_level` returns None when a gene has no staining rows, and takes the max over tissues rather than an average.

> *Without it:* Returning 0 would make a never-stained gene indistinguishable from one stained and found negative everywhere; averaging over tissues would dilute a single high-staining tissue below the level that matters, since one strongly staining normal tissue is the whole risk signal.

**7. hpa.py:100** - Required column names are resolved up front and a missing one raises KeyError, instead of tolerating an absent index.

> *Without it:* A missing column left as an absent index produces an empty table rather than an error, and an empty table is indistinguishable from a source that simply has nothing to say — a silently staining-free atlas passes every downstream safety check.

**8. hpa.py:178** - The Uniprot cell is split on commas and only the first accession is kept.

> *Without it:* The column can carry several accessions; used whole, the cell ("P12345, Q9Y6…") matches no surface accession at all and the gene drops out of the accession index silently. The first entry is the one the atlas treats as canonical (expected 19,300 of 20,162 genes carrying an accession).


## `car_pipeline/data/gtex.py`

**1. gtex.py:3** - The join route (direct symbol vs. Ensembl bridge) is recorded on each TissueProfile at join time rather than being inferred afterwards.

> *Without it:* After the fact a renamed gene and a directly matched one look identical, so the bridge count (expected 29 of 3,432 matched) could never be recovered — and that route is what a later rejection criterion counts.

**2. gtex.py:138** - Profiles are resolved by looping over the surface records in the protein's own order, not by driving the loop from the file's rows.

> *Without it:* Driving it from the file would make the recorded join route depend on which row happened to be reached first, so the same protein could be marked symbol-matched or bridge-matched depending on file ordering, and the bridge count a rejection criterion depends on becomes non-deterministic.

**3. gtex.py:93** - A protein with no row in the baseline is omitted from the returned profiles entirely rather than given a row of zeros.

> *Without it:* Not measured and measured at zero are different findings, and only one of them is reassuring: a zero row would present an unmeasured protein as silent in all 68 normal tissues, which is exactly the profile that clears the safety filter.

**4. gtex.py:46** - `silent_everywhere` requires every tissue to be below the threshold (default 1.0 TPM) rather than testing a summary statistic across tissues.

> *Without it:* A mean or median across 68 tissues is dragged below 1.0 by the many silent tissues even when one tissue expresses strongly, and that one expressing normal tissue is the entire safety concern.

**5. gtex.py:134** - The .gct header's own declared row and column counts are checked against what was actually read, and a mismatch raises IntegrityError.

> *Without it:* Checking costs nothing and is the only thing standing between a truncated download and a gene total that looks entirely reasonable — expected 59,033 genes x 68 tissues; a file cut short still parses and still yields per-tissue medians.


## `car_pipeline/data/tcga.py`

**1. tcga.py:3** - The cohort is downloaded in batches of 30 files rather than as one whole-cohort request.

> *Without it:* The whole cohort as a single request is roughly three quarters of a gigabyte that either arrives complete or fails complete, with no way to resume; one network blip discards the entire transfer instead of one recoverable batch.

**2. tcga.py:8** - The gene axis of every downloaded file is compared against the first file's axis before any values are stacked, and a mismatch raises IntegrityError instead of stacking.

> *Without it:* After stacking, columns are positional: two files disagreeing on gene order would misalign the 60,660-column matrix silently, and nothing downstream could detect it — every tumour/normal median after that point is read off the wrong gene.

**3. tcga.py:13** - The single metastatic sample is kept as its own category instead of being folded into the tumour or normal group.

> *Without it:* Folding it in would misstate that group's median, and the cohort is small enough that one sample moves it — 183 samples split 178 primary tumour / 4 solid normal / 1 metastatic, so the normal group is four samples wide.

**4. tcga.py:102** - The project id is embedded in the cache entry KEY as well as the filename (`cohort__TCGA-PAAD`), not just in the payload filename.

> *Without it:* Manifest paths are derived from the key alone, so namespacing only the payload leaves two cohorts sharing one manifest — the same overwrite bug with a longer fuse: TCGA-BRCA data validated against a TCGA-PAAD manifest that still looks valid.

**5. tcga.py:89** - Callers ask for the unqualified kind ("cohort") and `_entry` performs the project join in one place, raising KeyError when nothing matches.

> *Without it:* With the naming scheme spread across call sites, a call site that missed the join would silently match no entry and return an empty/absent cache slot rather than failing; the error surfaces at the lookup instead.

**6. tcga.py:58** - Rows whose gene id starts with `N_` are skipped at parse rather than parsed as genes.

> *Without it:* The leading rows of each STAR-Counts file are alignment tallies (N_unmapped, N_multimapping, …), not genes; parsed as genes they enter the matrix as enormous plausible-looking TPM values on a shared column position.

**7. tcga.py:163** - The file index checks the source's own declared pagination total against the number of hits received and raises rather than trusting the returned list.

> *Without it:* A truncated or paginated-short index yields a cohort built from fewer files than the project actually has; the run completes and every median is computed over a silently smaller sample set (expected 183 samples).

**8. tcga.py:209** - A batch of exactly one id is treated as the raw file itself; only multi-id batches are opened as a tar.gz archive.

> *Without it:* The GDC data endpoint returns the bare file for a single-id request and an archive otherwise, so unconditionally untarring fails on the final batch whenever the cohort size is 1 mod 30.

**9. tcga.py:271** - Each downloaded batch asserts that the archive held exactly as many files as were requested before any of them are parsed.

> *Without it:* An archive short a member would otherwise drop those samples from the matrix while the sample list and the counts printed downstream still look plausible; the mismatch is caught at the batch rather than discovered as an unexplained sample count.

**10. tcga.py:286** - A sample type outside the known three keeps its own raw label instead of being folded into a catch-all bucket.

> *Without it:* Folded, it would vanish from every per-category count while still inflating the sample total, and the two would stop adding up silently; the __main__ check for categories outside {primary_tumour, solid_normal, metastatic} exists to make that visible.

**11. tcga.py:350** - Surface proteins are matched by symbol first and only then through the Ensembl identifier bridge, with the route recorded per accession; a protein with no column is simply absent from the result.

> *Without it:* Giving an unmatched protein a column of zeros would make "not measured in this cohort" read as "measured at zero in the tumour" — an unexpressed-looking, safe-looking target. Expected 3,435 of the surface set match; the rest must stay absent.


## `car_pipeline/data/singlecell.py`

**1. singlecell.py:16** - A gene reading zero across every cell type is treated as evidence about the assay, never as evidence against the protein — this source separates compartments and is never allowed to refute a target.

> *Without it:* This assay drops transcripts that bulk measurement finds abundantly present: CEACAM5's peak across all groups here is essentially zero while bulk puts it near 300 transcripts and 409x normal. A zero-rejects rule would delete exactly the strong targets.

**2. singlecell.py:3** - Uses a tumour GEO series rather than the obvious atlas for this indication, and records the substitution explicitly in the stage 1 dataset list instead of leaving it implicit.

> *Without it:* Every pancreatic collection in the obvious atlas is normal, islet or diabetes. Reading it produces compartment means with no malignant compartment at all, so malignant-vs-stromal purity is computed against normal tissue and looks like a real answer.

**3. singlecell.py:46** - The dropout floor is DROPOUT_EPSILON = 0.001, not an exact-zero test; below it a group mean is called a capture failure rather than absence.

> *Without it:* Genes that bulk measurement puts in the hundreds of transcripts read here at a ten-thousandth, not at zero, so an exact-zero test misses the case entirely. Measured across the surface set: 267 genes read silent everywhere at 0.0, 357 at 0.001, 534 at 0.01.

**4. singlecell.py:38** - No dataset-specific constants live in the module: series, archive name, column names, category values and the compartment map are all declared on an AtlasSchema in the indication config and reached through self.atlas.

> *Without it:* A module global naming one submission is exactly what made this loader single-indication — a second atlas then needs an edit to library code rather than a declaration, and every column reference silently keeps pointing at the reference submission's names.

**5. singlecell.py:59** - Author branches not named in the compartment map are reported as OTHER rather than being forced into one of the five compartments.

> *Without it:* Forcing an unrecognised branch into the nearest compartment contaminates that compartment's mean with cells of a different kind, and the resulting purity ratio still looks like a clean number. OTHER holds 28,671 cells in the reference run — enough to move any compartment it was folded into.

**6. singlecell.py:66** - The stromal/immune comparator is the peak over fibroblast, immune and endothelial, not their average — one compartment expressing the gene is enough to disqualify it.

> *Without it:* Averaging dilutes a gene that is high in exactly one non-tumour compartment: HLA-DRA and CD74 come out at 0.07x malignant-over-peak against genuine targets at 193x (CEACAM6), 43x (CLDN18), 35x (MUC1), 20x (MSLN). Averaged, an immune-restricted gene reads as tumour-enriched.

**7. singlecell.py:92** - Peak_group uses np.nanmax and returns NaN only when the whole column is NaN, so groups holding no cells are skipped rather than participating in the maximum.

> *Without it:* A plain max propagates the NaN that empty groups carry, so every gene's peak comes back NaN and the entire surface set reads as unmeasured — or, if empty groups were zero-filled instead, an uncaptured cell type is scored as a measured zero.

**8. singlecell.py:137** - The group-means cache fingerprint carries derived_version (now 2) and the epsilon alongside the series and file, so a change in the derivation invalidates the artifact.

> *Without it:* Version 1 stored the identifier column as raw category codes rather than decoding them; without a version bump that cached artifact is reused and no lookup can ever match it — every protein joins as absent with no error anywhere.

**9. singlecell.py:166** - CART_NO_MATRIX_FETCH lets a deployment refuse to materialise the matrix at all, raising CacheError with instructions instead of fetching.

> *Without it:* On a container that ships only the derived summaries, fetching downloads 2.6 GB and expands it to 8.3 GB onto an in-memory filesystem; the instance is killed mid-job and reports nothing at all.

**10. singlecell.py:175** - The offline guard is evaluated before the not-gzipped passthrough branch, not after it.

> *Without it:* With the guard below that branch, a deployment that forbids fetching the matrix would happily download an uncompressed one instead — 844 MB onto an in-memory filesystem, the exact outcome the variable exists to stop.

**11. singlecell.py:192** - An archive whose name does not end in .gz is returned as the matrix directly rather than being run through the gzip expander.

> *Without it:* A CELLxGENE export is already an h5ad; feeding it to gzip.open fails on the first read instead of having nothing to do, so a perfectly usable atlas is unloadable.

**12. singlecell.py:246** - Obs/var columns are decoded through their category table (both the modern codes+categories group and the older shared __categories table) rather than read raw.

> *Without it:* Reading such a column raw yields the integer codes rendered as text, which look like perfectly good values and match nothing. The failure is entirely silent: every lookup misses, and a missed lookup is indistinguishable from a value the source never carried.

**13. singlecell.py:261** - Decoded categorical columns are built as a text-dtype array, not an object array, and negative codes become the empty string.

> *Without it:* An object array cannot be reloaded from the cache without allow_pickle, i.e. arbitrary deserialisation — which is not a thing a data cache should ever need; load() opens the npz with allow_pickle=False and would simply fail on it.

**14. singlecell.py:267** - Gene symbol and Ensembl id are each resolved through a declared field name that may be the index or a named column, instead of assuming the index holds symbols.

> *Without it:* The two submissions are mirror images — the reference indexes by gene symbol with Ensembl in a column, the CELLxGENE export indexes by Ensembl with the symbol in feature_name. Assuming one layout gives the other atlas Ensembl ids as 'genes', and every symbol join misses.

**15. singlecell.py:279** - The raw-counts matrix is looked up along the declared counts_path and a missing path raises by name, listing the keys actually present — there is no fallback to X.

> *Without it:* Falling back to X is the failure that does not announce itself: X is log1p(CP10K), so normalised values would be consumed as though they were counts and every downstream count would be wrong but plausible.

**16. singlecell.py:702** - Detection is defined on the raw integer counts layer, so counts is what gets read, rather than thresholding the normalised X.

> *Without it:* X is per-cell normalised log1p(CP10K), so a fixed threshold on it is not a fixed count threshold — raw depth across malignant cells runs from 96 to 9,642, meaning the same X value means a different number of transcripts in a shallow cell than a deep one.

**17. singlecell.py:318** - After reading the identifier column, the loader checks that at least half the values start with ENSG and raises with a sample of what it actually read.

> *Without it:* A bridge built from values of the wrong kind (codes, symbols, blanks) matches nothing and reports every protein as absent. The check makes it loud at load rather than invisible at join; the error quotes the count, e.g. 'N of 22,164 recognisable'.

**18. singlecell.py:343** - Before mapping, every level-3 cell type is checked to sit inside exactly one level-1 branch, and any straddling type aborts the build by name.

> *Without it:* If a cell type spanned two branches, every compartment figure below it would be a blend of two things — a malignant mean quietly averaged with stromal cells, producing a purity ratio that is wrong and unremarkable-looking.

**19. singlecell.py:374** - An atlas with no treatment column, or one whose split lacks the declared untreated label, is read anyway with an empty untreated subset — the untreated rows are present and empty rather than refused or filled.

> *Without it:* Absent is a third state, not a zero. Raising rejects a usable atlas; substituting zeros tells a reader the untreated cells were measured and silent, when in fact none were captured.

**20. singlecell.py:406** - Group means are means of np.expm1(X), i.e. the log1p is reversed per value before summing, not means of the stored values.

> *Without it:* X is log1p(CP10K); averaging the stored values gives a mean of logs, which is not the log of the mean and compresses exactly the high-expressing cells that make a target. The compartment marker checks (KRT19 malignant 14.1 vs fibroblast 0.02, COL1A1 18.3, PTPRC 16.7, VWF 16.2) only land on those figures in reversed space.

**21. singlecell.py:411** - Per-cell reversed totals are accumulated blockwise and stored alongside the means, purely as a scale check on the reversal.

> *Without it:* The transform is reversible only if the totals land where the stated normalisation says they should — median expected ~9,600 for CP10K. Without the check, a wrongly-scaled or already-linear X still produces means that look like plausible expression values.

**22. singlecell.py:449** - Normalise() writes NaN for groups with zero cells rather than dividing by a substituted denominator.

> *Without it:* A group with no cells has no mean. Zero-filling reports it as measured and silent, which is a different and much more reassuring statement than no cells of that type having been captured — and it drags any peak/average taken over compartments toward a fabricated zero.

**23. singlecell.py:549** - Per-cell malignant counts are built as their own cache entry keyed by a gene-set digest, not derived from the group means.

> *Without it:* The group-means artifact consumes the cell axis inside its accumulation loop and stores only 78 x 22,164 group means, so no conjunction over cells (co-expression, per-patient prevalence) can be recovered from it; and rebuilding it would invalidate the ranking stage for no reason.

**24. singlecell.py:569** - Wanted gene columns are found with a full-width lookup table over all var names, not by searching within each row's index array.

> *Without it:* The column indices in this file are stored in descending order within a row, so anything relying on the ascending order this sparse format usually carries (searchsorted-style membership) matches nothing and reports every gene as absent, silently.

**25. singlecell.py:590** - The donor column name is declared on the schema and both its absence from the schema and its absence from obs raise a named CacheError listing available columns.

> *Without it:* Reading a hardcoded name gave a bare TypeError three frames later when the column was absent, which named nothing — the operator sees a type error, not 'this atlas has no patient column, so per-patient prevalence cannot be measured'.

**26. singlecell.py:605** - The malignant builder guards the treatment column exactly the way the group-means reader does, tolerating a plainly-stored column or a missing untreated label.

> *Without it:* Unguarded, an atlas whose treatment column is stored plainly or whose declared label is absent builds its group means fine and then dies here — the same late, mislocated failure shape as the donor-column bug this replaced.

**27. singlecell.py:652** - Before writing into the uint16 count matrix, values are checked against np.iinfo(np.uint16).max and an overflow raises.

> *Without it:* Without the check a count above 65,535 wraps silently to a small number, so the most highly expressed gene in a cell reads as near-absent and every downstream detection rate for it is wrong but plausible.

**28. singlecell.py:662** - Savez_compressed is given an open file handle, not the temp path.

> *Without it:* Passing the path makes numpy append its own .npz suffix, so the cache commits a file that is not the one it just wrote — the entry validates against an empty or stale temp file.

**29. singlecell.py:730** - MalignantCells.subset slices the columns already loaded rather than triggering a fresh derivation for the smaller gene set.

> *Without it:* A smaller gene set is a slice of this one, not a different derivation; re-deriving streams 8.3 GB again to answer a question the loaded matrix already contains.

**30. singlecell.py:738** - Subset records requested genes that have no column in `missing` instead of rejecting them, mirroring what the builder does.

> *Without it:* Some requested genes have no row in the matrix at all, so a subset that refused them would behave differently from the derivation it is a subset of — the same gene list would succeed at build time and fail on a narrowing.

**31. singlecell.py:759** - The atlas join is symbol-first with an Ensembl-identifier bridge for renamed genes, and the route taken (JOIN_SYMBOL vs JOIN_ENSEMBL_BRIDGE) is returned with the column.

> *Without it:* A symbol-only join drops every gene the atlas records under an older or newer name and reports it as having no row; this join reaches the two heaviest components in the score. An unrecorded route also cannot be audited — the run prints the bridged count for exactly that reason (3,311 surface proteins matched, 185 with no row).

**32. singlecell.py:773** - Ensembl keys are stripped of their version suffix (split on '.') before being indexed, and the first occurrence wins.

> *Without it:* Versioned ids (ENSG00000105388.12) never compare equal to the unversioned ids the bridge is queried with, so every bridged lookup misses and the renamed genes the bridge exists to rescue stay reported as absent.

**33. singlecell.py:857** - Genes are ranked by the largest value any cell-type group reaches (nanmax over group means), not by a population mean across all cells.

> *Without it:* A mean over every cell is dominated by the nuclear transcripts this assay always captures, which say nothing about any cell type in particular; the sanity check is that the acinar enzymes CTRB1, CPA1 and PRSS1 come out at the top of 22,164 genes, which a population mean does not deliver.

**34. singlecell.py:8** - Nothing in this module ever holds the whole matrix: the archive is streamed on download, streamed again on expansion in 4 MiB blocks while hashing and fsyncing, and read back in row blocks of 8,192 for both the group-means accumulation and the malignant-cell pass.

> *Without it:* The obvious implementation, reading the h5ad in whole and grouping by cell type, is not a slower path but one that does not complete: 8.3 GB expanded against the memory available. Anyone reinstating a whole-matrix read, or raising the row block to go faster, gets an OOM kill mid-aggregation rather than an error naming the cause.

**35. singlecell.py:12** - The authors' annotations are used as given, with one editorial act: the two author immune branches are merged into one compartment, and every cell-type to branch to compartment mapping is printed line by line so that act is visible in the run output rather than buried in a dict.

> *Without it:* The compartment axis is otherwise the authors' annotation plus one silent re-annotation nobody can see. A reader comparing the run against the source publication finds one fewer immune branch with nothing explaining where it went, and a future collapse of two branches that do NOT describe the same compartment lands the same way. COMPARTMENT_ORDER holds six labels against a longer list of author branches, so the merge is not inferable from the constant alone.

**36. singlecell.py:116** - The registry key for this source stays the class name; the accession is per-run detail reported through series_name, not baked into the dataset identity. Cache entries are namespaced by atlas tag separately.

> *Without it:* Putting the accession in the dataset identity makes the same source a different dataset for every indication, so the Stage 1 row for the single-cell atlas changes name whenever the atlas does and can no longer be tracked as one required source across indications.

## `car_pipeline/data/depmap.py`

**1. depmap.py:3** - Files are taken from the archived Figshare release rather than the live DepMap portal.

> *Without it:* The portal serves newer releases from behind an interactive verification step that an unattended run cannot pass, so a portal fetch hangs or returns a challenge page in an automated run; the release (24Q4) has to be pinned anyway.

**2. depmap.py:7** - Download targets are resolved by file NAME from the article listing on every fetch; the numeric file identifier is never stored in the repository.

> *Without it:* Numeric identifiers change whenever a file is re-uploaded, so a hardcoded id either 404s or, worse, resolves to a different file — breaking or corrupting the fetch silently. Names survive a re-upload.

**3. depmap.py:178** - An empty cell in CRISPRGeneEffect.csv becomes NaN, never 0.

> *Without it:* Zero is a real and meaningful score here — it says the knockout did nothing — so writing it in place of a missing screen converts "nobody looked" into "looked and found no dependency", a claim the data does not support. With the overall median effect at -0.034, injected zeros also sit right on top of the real distribution and are invisible in aggregate.

**4. depmap.py:60** - Dependency accounting is NaN-aware throughout: `median_effect` uses nanmedian and returns NaN for an all-unscreened gene, and `dependent_lines` is reported against `screened_lines`, not against the lineage total.

> *Without it:* Of 68 models in the Pancreas lineage, 21 were never screened and only 47 have data; dividing a dependency count by 68 understates every gene by a third (KRAS is dependent in 43/47, not 43/68), and a plain median over a column containing NaN returns NaN for genes that do have data.

**5. depmap.py:40** - Only the derived per-lineage matrix is namespaced by lineage; the raw DepMap CSVs stay shared, unnamespaced cache entries.

> *Without it:* The raw files are pan-cancer and identical for every indication, so namespacing them re-downloads the same large CSVs per lineage; conversely, leaving the derived matrix unnamespaced lets a Pancreas-derived matrix be served to a Lung run as a valid cache hit — a complete, plausible dependency table for the wrong lineage.

**6. depmap.py:168** - Gene column labels are truncated at the bracketed identifier (`c.split(" (")[0]`) rather than used as written.

> *Without it:* The header reads "KRAS (3845)", so an untrimmed label matches no gene symbol at all and the symbol-only join returns zero of the expected 3,169 surface matches — an empty intersection that reads as "no surface protein is a dependency".

**7. depmap.py:222** - The join to this source is symbol-only; the Ensembl identifier bridge used for TCGA and GTEx is deliberately not applied.

> *Without it:* This source carries no stable gene identifier, so a bridge lookup has nothing to key on and would either match nothing or match through a stale symbol mapping; the honest consequence is that 327 surface proteins have no row here and are absent rather than bridged (3,169 matched of the surface set).

**8. depmap.py:93** - The published md5 from the article listing is passed to the downloader and verified, rather than accepting whatever bytes arrive.

> *Without it:* A truncated CRISPRGeneEffect.csv still parses as valid CSV and still yields a lineage matrix — just with fewer models or fewer gene columns — so without the checksum the only symptom is a screened-line count that is quietly below the expected 47.


## `car_pipeline/data/genespan.py`

**1. genespan.py:3** - Gene span is measured and reported alongside co-expression figures rather than corrected out of the per-cell detection rate.

> *Without it:* Detection rate tracks gene length harder than it tracks expression: rank correlation with genomic span is +0.68 against +0.20 for bulk tumour expression, and the span effect holds inside every quartile of expression (the atlas was quantified against a pre-mRNA reference, so intronic reads count and intronic content scales with span). Without the span reported beside it, a co-expression overlap that is high only because both genes are long is read as biologically high.

**2. genespan.py:15** - GENCODE span is not declared as a Stage 1 required dataset and gets no connector registration — it annotates, nothing gates on it.

> *Without it:* Registering it like the other sources makes the pipeline block on, and score availability against, a dataset no stage consumes: an annotation-only source drags the availability score down and can stop a run that has every dataset it actually needs. A Stage 1 row is owed only if a later stage ever gates on a span-derived quantity.

**3. genespan.py:61** - When a symbol appears more than once in the annotation, the longest annotated span wins (max), not the first, last, or averaged entry.

> *Without it:* A symbol on more than one contig is a real thing in this annotation; taking whichever row was seen first (or averaging them) understates the span and makes a long gene look short, so the length confound it exists to expose goes unflagged. The maximum is the conservative reading — it is the figure the capture artefact would scale with.


## `car_pipeline/data/antibodies.py`

**1. antibodies.py:5** - Chain roles (heavy, light, antigen) are read from the curated structure summary rather than inferred from each entity's free-text description.

> *Without it:* Without the curated table the heavy/light/antigen assignment has to be guessed from entity description text, which the specification named as the weak link of the whole approach — a wrong guess silently mislabels which chain is the binder and which is the target.

**2. antibodies.py:8** - The summary is fetched from the site's REST route, not written off as unreachable because an earlier probe of the front end came back empty.

> *Without it:* The front end that serves the CSV is a JavaScript application; an earlier probe read that shell as the source being unreachable and would have dropped a live source. The application has a REST interface and this URL is one of its routes.

**3. antibodies.py:12** - Named therapeutics are admitted as usable binders on their variable-region sequences alone, without requiring a solved structure.

> *Without it:* Restricting binders to deposited structures discards every clinically staged antibody whose sequence is published but whose structure is not solved — a construct is built from a sequence, so requiring a structure throws away usable binders for no gain.

**4. antibodies.py:17** - The compound Target field is split on `;` (bispecific antigens) and `/` (synonyms) and matched by exact token, never by substring.

> *Without it:* Substring matching puts MUC16 and MUC18 into the MUC1 bucket and matches `CLDN1` inside `CLDN18`, and this table carries CLDN1 and CLDN18 as separate targets. Measured field values include `CEACAM5/CD66e`, `MUC1/PEM/EMA`, `CLDN18;CD3E` and `IAP/CD47;CLDN18`.

**5. antibodies.py:76** - Entry identifiers are normalised to the four-character short form before joining, stripping the `pdb_0000` prefix, and the extended id is retained on each instance for traceback.

> *Without it:* The curated summary keys on the extended form (`pdb_00004f3f`) while the structure search returns the short form (`4F3F`); joining one directly against the other matches nothing, and matches nothing silently — every target gets an empty candidate list that reads exactly like "no binder exists". It did exactly that for all 200 targets until a positive known answer was added to catch it.


## `car_pipeline/data/structures.py`

**1. structures.py:3** - Structures are retrieved by exact UniProt accession cross-reference, never by full-text search on the gene symbol.

> *Without it:* Measured while writing the specification: a full-text query for one pool target returned 369 entries whose top hits were a bacterial RNA chaperone, a sulfur transferase and a photosystem supercomplex — an HTTP 200 carrying valid JSON and a plausible integer that was entirely spurious. Name matching is how this project has produced its worst answers.

**2. structures.py:11** - An empty result list is returned only for HTTP 204 with a zero-length body; the obvious repair of catching the JSON parse error and returning "no structures" is rejected and every other shape raises.

> *Without it:* A no-hit query answers 204 with an empty body, not 200 with a count of zero, so parsing the body raises. Swallowing that exception makes a broken query indistinguishable from an honest absence, and downstream records the empty list as "the literature holds nothing for this target" — a claim about biology manufactured out of a transport failure.

**3. structures.py:119** - Both directions of the zero contract are checked: a 204 that carries any bytes raises, and any non-204 status with an empty body raises.

> *Without it:* Only 204 may be empty. Accepting an empty body on a 200, or a 204 with content, lets a malformed or partial response pass as a validated zero for that accession.

**4. structures.py:94** - Transport failures (dropped connection, timeout, URLError) are retried a fixed 3 times with a 2.0-second escalating backoff and then raised — never converted into an empty list.

> *Without it:* Returning an empty list on a dropped connection is recorded downstream as "the literature holds nothing" for that target, turning a network blip into a biological conclusion that looks identical to a real negative.

**5. structures.py:133** - When the reported `total_count` exceeds the rows actually returned, the call raises instead of returning the first page.

> *Without it:* A page cap of 500 rows that quietly truncates undercounts both the entries and the candidate binders drawn from them, and the result looks like a protein with fewer structures rather than like a truncated read.

**6. structures.py:156** - Each entry summary records an `is_model` flag set when any experimental method contains THEORETICAL.

> *Without it:* Without the flag a computed model is reported as retrieved experimental evidence, presenting a prediction beside a measurement — which §1 forbids.


## `car_pipeline/data/domains.py`

**1. domains.py:3** - Every residue range is read from the fetched entry's own feature annotation; a part whose feature is missing raises instead of falling back to a remembered coordinate.

> *Without it:* A hard-coded coordinate that is close enough to look right produces a construct that assembles, translates and reports cleanly while being the wrong protein segment — the failure mode this project keeps finding. No coordinate in the module is transcribed from memory.

**2. domains.py:10** - Provenance is split into `proteome` (accession + range + release pin, re-derivable) and `synthetic` (designed literal, named), and synthetic parts are never given a fabricated accession.

> *Without it:* Blurring the two classes lets a designed linker or skip peptide be presented as database-derived evidence, so a reader cannot tell which parts are re-derivable from a source and which were written down by a designer.

**3. domains.py:16** - Parts are fetched one accession at a time rather than by widening the existing proteome query.

> *Without it:* Widening the proteome query invalidates that cache and forces a re-fetch of roughly twenty thousand entries to obtain eight sequences.

**4. domains.py:46** - The hinge is not read from a "hinge" feature (none exists); its length is a stated design choice of 45 residues and its position is derived as the segment immediately preceding the annotated transmembrane start, with a raise if it does not fit.

> *Without it:* There is no hinge annotation to look up, so the obvious move is a remembered literal range — which drifts from the actual entry. Deriving the end from the annotated transmembrane start keeps the stalk anchored to the entry, and the explicit check stops a 45-residue stalk from silently running off the N-terminus when it cannot fit before the transmembrane segment.

**5. domains.py:168** - `_slice` converts the annotation's one-based inclusive coordinates with `sequence[start - 1:end]` rather than slicing with the stated numbers directly.

> *Without it:* Feature ranges are stated one-based inclusive; a direct Python slice drops the first residue of every part and still yields a plausible-looking sequence of nearly the right length, so the error never announces itself.

**6. domains.py:218** - Caspase-9 is taken from `card["end"] + 1` to the end of the sequence — the CARD is removed, and its boundary is read from the entry rather than assumed.

> *Without it:* Including the CARD gives a constitutively active caspase-9 instead of a dimeriser-inducible safety switch, so the switch would fire without the drug. An assumed CARD boundary cuts the domain in the wrong place while still producing a sequence.

**7. domains.py:222** - FKBP12 is sliced from its annotated Chain feature rather than taken as the whole entry sequence.

> *Without it:* The mature protein starts after the initiator methionine, so the full entry sequence carries a leading residue that is not in the part. This was the one part bypassing the feature lookup that every other part uses.

**8. domains.py:243** - The adaptor receptor's anti-tag binder declares a size (240 aa) with an empty sequence, and consumers branch on `supplied` / `declared_residues` instead of on sequence length.

> *Without it:* The binder targets a tag, not a human protein, so it is not retrievable from the proteome, and no anti-tag antibody exists in the cached structural set either — 0 of 21,914 SAbDab entries name fluorescein as antigen. Writing a plausible-looking scFv would place an unverified binder inside a construct that reads as designed, the single thing this stage exists to prevent. The declared 240 aa is the median of the 30 real scFvs Stage 5 retrieved (705-765 bp, median 720 bp), which is what the payload budget actually needs.

**9. domains.py:256** - Stage 6 propagates `binder_supplied=False` for any construct using the adaptor binder and withholds the amino-acid sequence rather than emitting a fabricated one.

> *Without it:* Without the flag the construct emits as a complete designed sequence, and a downstream reader cannot tell the one unverified binder from the retrieved ones.

**10. domains.py:86** - Part.described defines what it means for a part to say where it came from per provenance class rather than uniformly: synthetic needs only a name, structure only an accession, proteome an accession plus both bounds. It is the property criterion K4 stands on.

> *Without it:* The obvious uniform rule, requiring an accession with a residue range for everything, rejects every structure-derived part, because a part retrieved from PDB entities carries an accession and no residue range at all -- the anti-tag binder is exactly that case -- so K4 would trip on a correctly sourced construct. Loosening it the other way to accession-only lets a proteome part with no range pass K4. The criterion's own spec wording is narrower than the implementation, so the code cannot be reconstructed from the criterion either.

## `car_pipeline/data/trials.py`

**1. trials.py:3** - The registry count is reported as the number of studies that *mention* the antigen, never as the number of studies testing a binder against it.

> *Without it:* The query is free text over the registry, so reporting the count as trials-against-the-antigen is a claim the registry cannot support. This project has already been bitten once by a free-text search returning a plausible integer that was entirely spurious.

**2. trials.py:9** - The terminated/withdrawn/suspended tally is surfaced as a read-before-dosing signal, not reported as evidence of a safety problem.

> *Without it:* Stopped trials have many causes other than harm; presenting the tally as proof of a safety problem overstates it, while dropping it hides that something already happened to people.

**3. trials.py:16** - The stopped, phase and CAR tallies are computed only over the studies on the single returned page, and every row carries `returned` and a `truncated` flag saying so rather than presenting the tallies as registry-wide.

> *Without it:* `total` is registry-wide but the tallies are not; treating them as complete understates terminations for any high-volume antigen. Paginating an antigen with 26,605 studies just to count its terminations is not proportionate to what the number is used for, so the tallies are declared a floor instead.

**4. trials.py:22** - The query term is the bare gene symbol only, and a zero is documented to mean "no study mentions this symbol", never "this antigen is untried".

> *Without it:* Synonyms are not searched and the undercount is large rather than marginal: measured live, `CLDN18` returns 3 studies while `Claudin 18.2` returns 156. Reading a low count as clinical novelty would invert the actual picture for that antigen.


## `car_pipeline/data/coverage.py`

**1. coverage.py:123** - The evidence class records the best tier of measurement available and never what the measurement said: staining that came back clean is PROTEIN_CONFIRMED, in the same class as staining that came back strong.

> *Without it:* Downgrading negative results makes 'nobody looked' and 'somebody looked and saw nothing' indistinguishable, and only one of those is reassuring — a protein with clean staining would be dumped into DATA_INSUFFICIENT alongside proteins with no protein-level data at all.

**2. coverage.py:158** - Nothing is dropped: every protein in the surface set leaves with a class, unmeasured ones landing in DATA_INSUFFICIENT and being broken down rather than filtered out.

> *Without it:* An unmeasured protein is not a refuted one; discarding them removes the only signal that the group exists and breaks the sum-of-classes == surface-set check (expected 1945 PROTEIN_CONFIRMED + 1491 RNA_SUPPORTED + 60 DATA_INSUFFICIENT, with the unresolved 60 splitting 14 polymorphic immune loci / 13 retroviral envelopes / 33 other).

**3. coverage.py:11** - Join routes are recorded into join_paths at the moment the join happens (accession hit vs symbol fallback vs ensembl bridge), never reconstructed afterwards by comparing symbols.

> *Without it:* Reconstruction gets the transcript baseline wrong: its own key is an identifier rather than a symbol, so a row reached by its primary key is misread as one reached by a fallback. The 'bridged' count — rows reachable only after the symbol failed — comes out inflated, and that route is what a later rejection criterion counts.

**4. coverage.py:58** - Polymorphic immune loci are matched with a fully anchored pattern over the real naming scheme (HLA-.+, KIR[23]D\w*, TR[ABGD], TR[ABGD][CVDJ]\d*) instead of a leading-prefix test.

> *Without it:* A bare prefix test on 'TRA' also captures triadin, the TRAF-interacting proteins and the arginine transporter, and 'KIR' captures the kirre-like family — none of which are polymorphic immune loci. The failure is nearly invisible: the count can still come out right while the test is wrong, because a miscategorised gene only shows up if it happens to land in this class.

**5. coverage.py:36** - 'Silent in every normal tissue' is a per-tissue test — every tissue below SILENT_TPM = 1.0 — not an average/median across tissues and not an exact-zero test.

> *Without it:* An aggregate over the tissue profile lets a gene with one strongly expressing normal tissue pass as silent, which is exactly the on-target/off-tumour case the flag exists to catch; the threshold sits at 1 TPM rather than zero. The report expects 472 measured-and-silent proteins, so an aggregate or exact-zero rule changes that count without looking broken.

**6. coverage.py:16** - The coverage report is built deliberately without the two heaviest sources. It takes the surface set, the tissue atlas indexes, the normal baseline and the tumour join, and reads neither heavy archive.

> *Without it:* Wiring the heavy sources in for completeness makes every coverage run expand an eight gigabyte archive to produce counts it cannot change: the evidence classes and join paths come from the four light sources alone. Unrecorded, the omission reads as an oversight and the next person adds them back, paying the decompression on every run for identical numbers.

## `car_pipeline/data/availability.py`

**1. availability.py:1** - Dataset status is resolved against the local cache manifests, never by probing the network.

> *Without it:* Reachability answers the wrong question: a server that responds but whose data was never fetched cannot feed a stage, and a source already fetched works whether or not anything is reachable right now. A network probe reports both backwards, so a run starts that cannot complete (and stalls on an offline machine that has everything it needs).

**2. availability.py:8** - Only manifests are read here; the deep integrity check — row counts, gene axes, checksums — is deliberately deferred to load time.

> *Without it:* Running the deep check at status time attributes any failure to the whole source rather than to a specific file, so a single bad member reports the entire dataset as unavailable with nothing pointing at which file to refetch.

**3. availability.py:25** - The two binder sources declared by Stage 1 are deliberately left out of the CONNECTORS map so they resolve to NOT_CONFIGURED, rather than being registered to a connector that does not exist.

> *Without it:* Adding a name before its connector exists reports the source as merely uncached — a fetchable thing someone just has not fetched — when in fact no code exists to fetch it. 'No connector' and 'connector whose data will not read' are different findings and must not collapse into one.

**4. availability.py:40** - Connector keys are cross-checked against KNOWN_DATASET_NAMES at import and any orphan raises RuntimeError, instead of the two name lists being trusted to agree.

> *Without it:* A connector registered under a name no stage emits leaves that dataset reported as having no connector at all while the connector sits unused. Both halves look plausible on their own — the dataset row reads not_configured, the connector class exists — so nothing surfaces the typo without the explicit check.

**5. availability.py:60** - AVAILABLE requires every entry the source declares to be present with a matching manifest; a partial cache resolves to UNREACHABLE.

> *Without it:* Treating any present file as availability lets a source missing one of several entries pass the pre-run gate, and the stage then dies partway through with the run already committed.

**6. availability.py:88** - The merged binder dataset row was split into two rows even though it lowered the availability score from 0.857 to 0.750 (6 of 8 blocking datasets available).

> *Without it:* Nothing became less available at the split — the merged row was never connected either — so counting two unconnected sources as one understated the gap and made 0.857 read as better coverage than exists. The numerator moves only when a connector is actually built.
**7. availability.py:82** - The release pins live beside the connectors that hold them, and the per-indication ones are marked as such.

> *Without it:* Six of the nine pins are module constants shared by every indication and three come from the indication config, so a flat list would read as though the cohort and the atlas were platform-wide facts. A package identifies the exact data a candidate was screened against, and half of that data is chosen per indication.


---

# Stages


## `car_pipeline/stages/stage1.py`

**1. stage1.py:46** - The single-cell atlas dataset row is named generically ("Single-cell tumour atlas") instead of carrying the reference accession/series it was originally declared with; which series was actually read is reported per run rather than asserted in the registry.

> *Without it:* Embedding the accession ties the dataset registry to one submission: a second indication that reads a different series no longer matches the registered name and is reported as having no connector at all, while its connector sits there unused.

**2. stage1.py:57** - The binder sources are declared as two separate rows — PDB for structures and SAbDab for antibody chain/numbering annotation — rather than one merged binder-data row, even though both were declared together while both were assumed unreachable.

> *Without it:* A single row cannot carry two statuses. The structures and the antibody annotation over them are separate services that connect, cache and fail independently, so merged, one of them being available would report the other as available too.

**3. stage1.py:34** - The dataset list is mode-dependent, not static: the screening sources (UniProt membrane topology, GDC TCGA bulk expression, DepMap dependency) are emitted only in DISCOVER mode, while the shared sources (HPA, single-cell atlas, GTEx) are emitted in both because they describe the antigen's behaviour rather than nominate it.

> *Without it:* A supplied target means there is nothing to screen. Emitting the screening rows anyway makes them blocking (required=True) entries in a validation run, dragging data_availability_score down for data the run never needed and making a complete run look gap-ridden.

**4. stage1.py:71** - KNOWN_DATASET_NAMES is exported as a frozenset built from all three dataset tables so the status resolver can cross-check its own registry against the names this stage can actually emit, instead of each side independently keeping its own name strings.

> *Without it:* Without the cross-check a name that drifts on one side fails silently: the resolver has no entry matching the emitted name and quietly reports a dataset that is configured and readable as having no connector.

**5. stage1.py:86** - The project id appends a random 8-hex uuid suffix to the cancer-type slug and second-resolution timestamp rather than relying on slug plus timestamp alone.

> *Without it:* The timestamp alone collides on back-to-back builds in the same process — two distinct projects share one id, and whichever writes second overwrites the first's outputs.

**6. stage1.py:23** - The construct size budget subtracts BACKBONE_OVERHEAD_KB = 1.2 from the vector payload limit rather than using the declared payload limit (default 4.7 kb) directly as max_construct_kb.

> *Without it:* The regulatory and structural elements the construct carries before any binder or signalling domain is added occupy 1.2 kb. Using the raw limit hands later stages a plausible-looking 4.7 kb budget instead of the real 3.5 kb, so designs pass the size check and do not fit the vector.

**7. stage1.py:117** - Allowed_car_formats is built by filtering CARFormat.AUTO out of the enum rather than passing the whole enum through.

> *Without it:* AUTO is a sentinel meaning a later stage picks the format, so it cannot also be one of the things that stage picks between — leaving it in lets the selector "choose" AUTO and return a request-to-choose as its answer.

**8. stage1.py:128** - Build_spec deep-copies the ProjectInput onto the spec instead of storing the caller's object by reference.

> *Without it:* A later stage writes a discovered target back onto the spec's inputs; with a shared reference that write mutates the shared indication config itself, so every project built after it in the same process silently starts in validation mode with the previous run's antigen.

**9. stage1.py:155** - _resolve gates on importlib.util.find_spec — whether the resolver module exists — and then imports it unguarded, rather than wrapping the import in try/except ImportError.

> *Without it:* Catching import failures turns a broken dependency anywhere in the connector chain into a clean report that every dataset is not_configured: a wrong answer that looks like a valid one. The find_spec form only short-circuits when the resolver genuinely has not been built yet, where every dataset defaulting to not_configured is correct.


## `car_pipeline/stages/stage3.py`

**1. stage3.py:6** - The three numbers this stage emits per protein — attractiveness, normal-tissue risk, evidence confidence — are never combined into a single score; risk is a gate and confidence is a description of how much measurement stands behind the other two.

> *Without it:* Folding risk into a composite lets a high tumour-side score buy off an unmanageable normal-tissue reading, and folding in confidence lets a well-measured mediocre target and a barely-measured strong one land on the same number. Attractiveness is tumour-side only by construction.

**2. stage3.py:68** - The margin denominator floor (c3_baseline_floor_tpm = 0.1 TPM) is listed inside SATURATION with the other free parameters rather than kept as a fixed constant, so the twelfth rejection criterion perturbs it along with the saturation points.

> *Without it:* Treating it as a constant exempts it from the sensitivity sweep even though, as the comment states, "it moves the top of the ranking as hard as any saturation point does" — the sweep would then certify a ranking whose most influential parameter was never varied.

**3. stage3.py:510** - The tumour/normal fold change floors its denominator at the detection limit (0.1 TPM) instead of treating a zero normal reading as an unbounded ratio.

> *Without it:* Dividing by zero or near-zero "awards a perfect margin to proteins absent from the tumour as well, where there is no margin because there is nothing there" — a protein at 0 TPM in tumour and 0 in normal would top the c3 component (saturation c3_fold = 64). Flooring keeps a genuine not-found-in-normal reading favourable while leaving the numerator to decide whether anything is actually present.

**4. stage3.py:878** - Source disagreement is measured as a departure from the systematic median log-ratio offset between the two normal denominators (tolerance FOLD_DISAGREEMENT = 2.0x), not as a departure from parity — the offset is computed first over all proteins carrying both folds, then each protein is tested against that.

> *Without it:* The population baseline and the cohort adjacent-normal describe different tissue and disagree by construction and by a fairly consistent amount, so "a parity test flags everything and discriminates nothing" — every protein comes back flagged and the flag carries no information.

**5. stage3.py:82** - Field_elevation (baseline fold / cohort fold, i.e. how much more antigen sits in tumour-adjacent pancreas than in healthy population pancreas) is reported on every row and never enters any score.

> *Without it:* The cohort's adjacent normal comes from a resected pancreas that is typically inflamed and often already carries precursor lesions, so it expresses tumour antigens before any malignancy is involved; scoring the ratio would confuse where a margin comes from with how large it is. It is retained as an annotation because what stays in the patient after resection is field tissue, so an antigen elevated here is one the therapy will meet in tissue that is not being removed.

**6. stage3.py:121** - Organs the platform added beyond the reference criticality table (vascular, eye, mucosa, connective) are tracked in PLATFORM_ADDED_ORGANS and marked with [+] in the header rather than silently merged into ORGAN_TIERS.

> *Without it:* Without the mark a reader has to trust the whole tier table equally and cannot tell which criticality assignments are inherited from the reference and which were supplied here — including the tier-1 calls on vascular and eye that can single-handedly fail a target through the worst-organ maximum.

**7. stage3.py:126** - Cultured fibroblasts, EBV-transformed lymphocytes and 'N/A' are excluded from risk entirely rather than mapped to some organ.

> *Without it:* They are not normal tissue. Mapping them to a nearest organ (connective, marrow_and_blood) lets a cell-line artefact set the criticality-weighted worst-organ risk for a target and fail it on tissue that does not exist in a patient.

**8. stage3.py:130** - Both tissue→organ tables (BASELINE_ORGANS, ATLAS_ORGANS) are exhaustive exact-match dictionaries; no substring or keyword test is used anywhere in the mapping, and an unmatched label falls through to the RiskModel.fall_through audit set instead of being guessed.

> *Without it:* Substring matching is what "put the adrenal at kidney criticality and folded the kidney into the brain" — 'Adrenal_Gland' matching on renal and 'Kidney_Cortex' matching on cortex. Both mistakes move an organ between tier 1 (w=1.0) and tier 2 (w=0.6) and silently rewrite the risk gate for every protein.

**9. stage3.py:276** - Staining levels are placed on the transcript axis by measurement — every organ carrying both a staining call and a transcript value contributes one observation, each level becomes the median of its own population, and that median is then scored through the same continuous function the transcript side uses — rather than asserting an evenly spaced ordinal (0/0.33/0.67/1.0).

> *Without it:* The two sources describe the same organs in different units, so "until they sit on one axis the worst-organ maximum is comparing incomparable numbers" — the max in worst_organ would pick whichever source happened to have the larger arbitrary units, not the riskier organ.

**10. stage3.py:293** - CalibrationCurve.score returns exactly 0.0 for level 0 (Not detected) instead of the calibrated median TPM that level's paired population actually sits at.

> *Without it:* The calibrated value for level 0 is a central estimate of a wide distribution, so using it as a risk term means "an observation of absence raises measured risk — which it cannot." Scoring zero rather than dropping the organ still leaves the organ measured, so the per-organ maximum lets a positive transcript reading stand rather than being cancelled.

**11. stage3.py:340** - Calibration pairs are formed exactly the way risk is computed — level = maximum across the organ's cell types, transcript = maximum across the tissues mapping to that organ — not from raw per-tissue pairs.

> *Without it:* Calibrating the same aggregated quantity that gets scored is what lets the curve absorb whatever inflation the double maximum introduces; calibrating on raw pairs and then scoring on maxima leaves a systematic upward bias in every risk score with nothing to correct it.

**12. stage3.py:1082** - The curve publishes its own quality — per-level n, Q1/median/Q3, a monotonicity check, and a rank-based separation statistic where 0.50 means the two levels say nothing about each other — plus the explicit header line that "the scale is real but weak; a one-level difference is not decisive."

> *Without it:* A calibrated median printed alone reads as a hard scale. The level populations overlap heavily, so a reader who treats a Low-vs-Medium difference as decisive is reading noise; the separation numbers are the only thing in the output that says how much of the ordinal is real.

**13. stage3.py:408** - The c3 margin denominator is the single bulk Pancreas column, not a median across the four pancreas entries the baseline names.

> *Without it:* Three of the four entries are cell-sorted fractions (Pancreas_Acini, Pancreas_Islets, Pancreas_Mixed_Cell) and one is bulk; the tumour side of the ratio is bulk, so a median across all four "mixes measurement types on one side of a ratio" and produces a fold change whose magnitude is partly an artefact of sorting.

**14. stage3.py:412** - The risk gate deliberately keeps reading all four pancreas entries and takes the worst of them, even though the margin component uses only the bulk one.

> *Without it:* For safety the question is "whether any pancreatic compartment carries the antigen, not what the organ averages" — averaging the sorted fractions into the bulk value hides an antigen concentrated in islets or acini behind an organ-level mean and clears a target that would be hit in that compartment.

**15. stage3.py:753** - A margin_label naming a column absent from the baseline raises KeyError rather than leaving the denominator unset.

> *Without it:* "Silently having no denominator would make the margin component unmeasured for every protein at once, and the evidence floor would then quietly drop a third of the universe" — c3 carries weight 0.25, so losing it pushes many proteins under MINIMUM_MEASURED_WEIGHT = 0.40 and they vanish from the ranking as below-floor rather than as a configuration error.

**16. stage3.py:477** - A malignant single-cell value at or below DROPOUT_EPSILON makes c1 unmeasured (value None, note 'below capture threshold'), never a score of zero.

> *Without it:* At or below the silence threshold "this assay is reporting its own capture failure, not the protein" — scoring it zero penalises a protein for a dropout and, worse, counts its 0.25 weight as measured, so the composite is renormalised as if the evidence were there.

**17. stage3.py:491** - C2 returns unmeasured when either the malignant value or the stromal peak sits at or below DROPOUT_EPSILON, instead of computing the ratio.

> *Without it:* "A ratio against a denominator this assay failed to capture is not a specificity measurement, however large it comes out" — a stromal dropout divides into an enormous apparent malignant-vs-stroma ratio that saturates the component (c2_ratio = 50) on nothing but a capture failure.

**18. stage3.py:554** - For the lipid-anchored class with no residue annotation, the residue term is excluded and a class default is used (1.0 with plasma-membrane confirmation, 0.6 without) rather than reading the zero annotated outward residues as zero accessibility.

> *Without it:* GPI-anchored proteins "report zero annotated outward residues by construction — there is no transmembrane segment to annotate around" — reading that as a measurement scores an entire, entirely accessible protein class at zero on the accessibility component.

**19. stage3.py:560** - The anchor-class default applies only when residues is None; a protein carrying both a lipid anchor and a transmembrane segment keeps its measured ectodomain size.

> *Without it:* The by-construction-zero argument "only holds while there is no annotation." A handful of proteins carry both, and for those the ectodomain really was measured — substituting the class default "would throw away the better evidence of the two" and hand a measured-small ectodomain the full 1.0 anchor score.

**20. stage3.py:557** - A non-anchored protein with no ectodomain annotation is unmeasured for c5 — never imputed with a default or a mean.

> *Without it:* "A protein nobody annotated must not outrank one measured and found small." Imputation converts an absence of annotation into a positive accessibility score and also counts c5's weight toward the 0.40 evidence floor.

**21. stage3.py:578** - C6 is unmeasured when the gene was never screened (screened == 0) or the effect is NaN, rather than scoring the missing dependency as zero.

> *Without it:* Zero is a real c6 value meaning 'no essentiality', so an unscreened gene would be reported as measured-and-non-essential; its 0.05 weight would count toward the evidence floor and toward confidence on a screen that never ran.

**22. stage3.py:629** - Per-organ expression scores are factored out of compute_risk into per_organ_scores so the pairing stage can take a conjunction across two proteins before the maximum, instead of reimplementing the tissue mapping.

> *Without it:* "A second implementation of this would be a second place for the tissue-mapping bugs to live" — the two copies drift and the pairing stage's risk stops meaning the same thing as the single-target risk it is compared against.

**23. stage3.py:656** - Where both sources cover an organ, the per-organ score is the maximum of atlas staining and baseline transcript, never the minimum or a mean.

> *Without it:* "A 'not detected' staining call must not cancel a positive transcript measurement for the same organ. That would understate risk" — the minimum lets one negative antibody call clear a target that transcript data says is expressed in a tier-1 organ.

**24. stage3.py:669** - Normal-tissue risk is the criticality-weighted maximum over organs (tier weights 1.0 / 0.6 / 0.3), not the mean.

> *Without it:* "One unmanageable organ disqualifies a target however clean the rest are" — averaging over ~21 organs dilutes a full-strength brain or heart signal to near nothing and clears a target that would be lethal in one tissue.

**25. stage3.py:851** - Worst_organ returns None when no source measured the protein in any organ, and `cleared` requires risk is not None — so unmeasured risk fails the gate.

> *Without it:* "Undefined risk is not low risk. It fails." Defaulting a missing risk to 0.0 clears every protein with no normal-tissue coverage at all, which are precisely the ones nothing is known about.

**26. stage3.py:703** - _confidence sums the weights this run actually scored with (the wts argument), not the module-level WEIGHTS defaults, and adds an evidence-class bonus of 0.3 / 0.15 / 0.0 on top of measured*0.7.

> *Without it:* On a perturbed-weight run the module defaults would produce "a confidence describing a different weight set from its own composites" — the reported confidence and the composite next to it would then be computed under two different models, with no sign of it in the output.

**27. stage3.py:738** - A criticality override naming an organ absent from ORGAN_TIERS raises KeyError listing the known organs, instead of being ignored as a no-op.

> *Without it:* "An override naming an organ that does not exist would do nothing at all, while the output header went on reporting it as an applied relaxation complete with its rationale. A safety default that only appears to have been changed is worse than one that was never touched."

**28. stage3.py:789** - Cohort medians are guarded: prevalence and tumour median are only taken when the primary-tumour group is non-empty, and every median passes through _finite so NaN/inf become absent rather than values.

> *Without it:* "An empty group has no median, and a not-a-number median is not a measurement of zero. Either would otherwise flow into the margin as a real reading and be scored as a target with no enrichment" — a NaN or 0 tumour median silently becomes a measured c3/c4 rather than a gap.

**29. stage3.py:870** - The bridged audit flag is `row.bridged or cell_bridged` — the single-cell atlas's Ensembl-bridge join counts as a bridge too, not just the coverage row's.

> *Without it:* The cell atlas join "feeds the two heaviest components" (c1 at 0.25 and c2 at 0.20, 0.45 of total weight), so omitting it "would exempt the most consequential bridge from the criterion that polices them" — a symbol-failed, accession-bridged join could carry nearly half a protein's score without ever being audited.

**30. stage3.py:910** - Tier_rank is assigned within each evidence class (PROTEIN_CONFIRMED, RNA_SUPPORTED, DATA_INSUFFICIENT) with unscored rows sorted last, rather than one global ranking on composite.

> *Without it:* Composites are renormalised over measured weight, so a protein scored on the bare 0.40 minimum is on the same 0–1 scale as one scored on all six components; a single global sort puts a data-poor protein at rank 1 above fully measured, protein-confirmed targets.

**31. stage3.py:956** - The configuration hash covers the tissue mapping tables themselves — BASELINE_ORGANS, ATLAS_ORGANS, EXCLUDED_LABELS — not only the tier assignments and overrides.

> *Without it:* "Adding or moving a single label shifts thousands of risk scores. A hash that did not move with it would let two materially different experiments compare as the same one, which is worse than publishing no hash at all."

**32. stage3.py:970** - The margin denominator column name is included in the configuration hash payload.

> *Without it:* "Swapping it moves every fold change and therefore the whole ranking, so a hash that ignored it would let two different experiments compare as one" — two runs against different normal denominators would be indistinguishable in the reproducibility header.

**33. stage3.py:977** - The measured atlas-level calibration curve is hashed into the configuration, rather than being treated as a fixed property of the code.

> *Without it:* It is "the measured curve, not an assumed scale" — derived from the pinned data releases, so it belongs to the experiment. "Two runs that calibrated differently must not compare as the same one," and without it a data-release change that moved every staining-derived risk score would leave the hash unchanged.

**34. stage3.py:1003** - The header prints the saturation set, weights, margin denominator and calibration actually passed into this run, falling back to module defaults only when nothing was supplied.

> *Without it:* Printing the module constants "describes the module, not the run that produced the output" — a perturbed sensitivity run would emit a header claiming the default parameters beside numbers produced by different ones, and the printed configuration hash would be the only surviving hint.


**35. stage3.py:597** - Attribution reconstructs the three reductions from the same inputs rather than being recorded as `compute_risk` runs.

> *Without it:* Instrumenting the scoring path makes the explanation a by-product of the thing it explains, so the two cannot disagree and the reconstruction proves nothing. Recomputing independently is what gives T1 its force: the attribution and the risk are derived twice from the same measurements, and the criterion asserts they meet to within 1e-12.

**36. stage3.py:604** - An organ whose two arms score equally is reported `TIED`, and a target reaching its maximum on several organs lists all of them.

> *Without it:* Both reductions are a maximum, and a maximum silently picks a winner among equals. `worst_organ` returns whichever organ dict iteration reached first, which is an artefact of insertion order and not a fact about the target. Reporting one organ where several tie would attribute the verdict to evidence that only shares the credit; 314 targets reach their maximum on more than one organ.

**37. stage3.py:246** - The attribution payload carries its numbers unrounded, unlike every other payload in this file.

> *Without it:* The record has to reconstruct the risk from itself - that is the whole claim of §3. Rounding a TPM to four places changes the score the baseline curve returns from it, so a reader recomputing from the served record would get a number that disagrees with the one served beside it, in the fifth decimal, for no reason they could see.

**38. stage3.py:660** - `RiskInputs` bundles the five objects attribution needs so one target can be explained on demand.

> *Without it:* The alternative is attributing all 3,466 targets during the screen and carrying the result on `Ranked`, which puts 65,077 organ rows in memory to serve the one gene a reader asked about. The facility is a read of state that already exists; making it a stored output would make it something the pipeline maintains.

## `car_pipeline/stages/stage4.py`

**1. stage4.py:13** - Per-organ risk scores are read from the ranking stage rather than recomputed here, so pairing a target with itself reproduces its single-antigen risk exactly.

> *Without it:* A second implementation of the tissue mapping is a second place for the mapping bugs to live, and the two copies drift silently; the self-pairing identity is what would have caught the drift.

**2. stage4.py:30** - The pool is the top 200 of the tumour-side ranking with risk ignored entirely — risk is not used to pre-filter pool membership.

> *Without it:* Risk is the thing pairing exists to fix, so filtering on it first leaves the stage unable to reach the risky targets a dual-antigen design is for; the stage would only ever be offered targets that already cleared alone.

**3. stage4.py:37** - A cell counts as carrying an antigen at >=1 captured molecule, but the whole analysis is also reported at thresholds 2 and 3 rather than committed to one number.

> *Without it:* The ordering of pairs is not stable across detection thresholds, so a single hard-coded threshold reports a ranking that is an artefact of the threshold and looks like a result.

**4. stage4.py:42** - COVERAGE_FLOOR is 0.02, set against the measured range rather than chosen as a round number.

> *Without it:* The best pair of known targets in this atlas reaches only 0.047 of malignant cells, so a floor at 0.05 or 0.10 eliminates every pair the stage exists to evaluate — the stage returns nothing and looks like it found nothing. Paired with criterion P16, which trips if even 0.02 admits nothing.

**5. stage4.py:50** - Patients with fewer than 100 malignant cells are excluded from the per-patient coverage denominator, and the exclusions are reported rather than dropped quietly.

> *Without it:* A proportion cannot be estimated from three cells; including such patients produces per-patient fractions of 0 or 1 that swamp the patient-fraction floor. Only 29 of 43 patients clear this, so silently dropping 14 would misstate what the number is over.

**6. stage4.py:54** - When fewer than 10 cells are positive for either member, the pair is marked unmeasured rather than scored zero coverage.

> *Without it:* Below this a positive fraction is not a measurement. Scoring it zero turns an assay capture failure into evidence against the target — a single-cell zero must never reject a target.

**7. stage4.py:97** - Combined pair risk is the minimum of the two antigens per organ, then the maximum across organs — not the maximum per organ.

> *Without it:* An AND gate only fires where both antigens are present, so the organ's risk is bounded by whichever antigen is scarcer there. Taking the maximum reduces the pair to its more dangerous member and ignores the architecture entirely, so no pair ever looks safer than its worst antigen and dual design appears worthless.

**8. stage4.py:104** - The per-organ minimum is used knowing it assumes perfect within-organ overlap (maximal co-expression), because that is pessimistic for safety and optimistic for coverage at the same time.

> *Without it:* Neither source carries a joint distribution over cells, so no measurement of within-organ co-expression exists to be had; any other choice would be inventing one. The two error directions coinciding is what makes this a legitimate bound to gate on rather than a guess.

**9. stage4.py:111** - Organs measured for neither member are left out of every variant rather than filled with a default.

> *Without it:* Filling them breaks the identity that pairing a target with itself reproduces its own single-antigen risk, which is the stage's only self-check against the tissue mapping; the ranking stage already tolerates organs nobody measured.

**10. stage4.py:129** - For an organ where only one member was measured, the missing member is assumed present (conservative keeps the known score; pessimistic charges 1.0) rather than the organ being credited as an absence or dropped.

> *Without it:* The pair's whole safety claim is an absence, so an unobserved absence cannot be credited to it. Crediting it lets a pair clear the ceiling on organs nobody looked at for one of its antigens.

**11. stage4.py:80** - The pessimistic variant (unresolved organs at full criticality) is computed and reported but is never the gate.

> *Without it:* It answers 'would this pair still clear if the unmeasured antigen turned out to saturate the organ', which is what 'clearance depends on an unresolved organ' means as a number. Gating on it would reject pairs for missing data; hiding it would let a reader miss that the clearance is contingent.

**12. stage4.py:142** - Pair risk is rounded to 4 decimals, the same precision the ranking stage stores single-antigen risk at.

> *Without it:* Comparing a full-precision pair risk against a rounded single risk reports a disagreement at the fifth decimal that is arithmetic rather than substance — a self-paired target would appear not to reproduce its own risk.

**13. stage4.py:203** - `escape` (1 - f_ab) is documented as a floor on the escape population, not an estimate.

> *Without it:* This assay drops transcripts, so the measured intersection understates the true one and this number overstates the true escape population; read as an estimate it makes every AND gate look worse than it is.

**14. stage4.py:231** - All pairwise double-positive counts come from one boolean matrix product in float32, cast to int64 after.

> *Without it:* Float32 counts integers exactly only far above this atlas's cell count, which is what makes the fast path safe; a larger atlas would silently return approximate intersection counts that still look like integers.

**15. stage4.py:265** - `p_b_given_a` is named for the antigen whose cells are the denominator: P(B|A) describes A's coverage, not B's.

> *Without it:* The obvious reading attributes the conditional to the wrong antigen, so a reader reports the retained fraction of the wrong member — a plausible-looking number about the other gene.

**16. stage4.py:439** - Per-patient double-positivity is accumulated as a count of patients at or above the coverage floor, rather than pooling all cells across patients.

> *Without it:* A pair that is double-positive in half the patients and absent in the rest pools identically to one that is uniform across every patient, and those are different products — the first treats half the cohort and reports the same number as the second.

**17. stage4.py:417** - `evaluate` raises if any pool member has no normal-tissue risk at all, even though the pool deliberately ignores the value of risk.

> *Without it:* The pool ignores the value of risk, not whether it exists. A member with no risk makes every pair containing it unresolvable, and silently so — 199 other pairs come back unexplained.

**18. stage4.py:305** - Pair confidence is the minimum of the two members' confidences (then discounted), never an average, and is never combined with risk or coverage into one score.

> *Without it:* A pair cannot be better evidenced than the least evidenced antigen in it; averaging lets a well-characterised antigen carry an unmeasured partner to a confident-looking pair. The two discounts — fraction of the organ union resolved for both members, and 0.75 when co-expression was not measurable — are the things pairing itself adds to the question.

**19. stage4.py:322** - `cleared` tests the conservative combined value, not the optimistic one.

> *Without it:* Using the optimistic value lets an unresolved organ clear a pair — the pair passes the ceiling precisely because nobody measured the organ that would have failed it.

**20. stage4.py:339** - `rescued` is a condition, not a magnitude: a member is rescued only when its own risk is above the ceiling and the pair's is not.

> *Without it:* A large risk movement that does not cross the ceiling buys nothing clinically, but reported as a statistic it reads as a success; the delta still reports how far it got, separately.

**21. stage4.py:365** - Coverage does not gate admissibility — only risk clearance and measurability do; `f_ab` never selects.

> *Without it:* Over the pool, f_ab's rank correlation with genomic span is +0.68 against +0.20 with bulk expression, and the confound reaches the joint quantity itself (+0.63 for f_ab against +0.08 for expression) and survives stratification by expression. A threshold on f_ab therefore admits and rejects partners substantially on how long their genes are. What is given up is stated: nothing now stops a pair with negligible overlap being recommended, so f_ab and its span percentile are reported per pair.

**22. stage4.py:375** - `coverage.measured` is still required for admissibility even though the coverage floor does not gate — measurability and magnitude are kept as separate questions.

> *Without it:* They are different questions: whether co-expression was observable at all versus whether it cleared a number. Dropping the measurability check alongside the floor would let a pair nobody could measure be recommended.

**23. stage4.py:388** - `build_pool` keeps exactly one entry per gene symbol, deduplicating accessions.

> *Without it:* Several symbols carry more than one accession; two pool members naming the same gene would pair with each other and report a perfect intersection that means nothing — the top-ranked 'pair' would be a gene with itself.

**24. stage4.py:511** - Every decision carries its own accession and its partner's, rather than leaving a consumer to look them up from the gene symbol.

> *Without it:* A symbol is not an identity in this proteome — several symbols carry more than one accession and build_pool keeps exactly one. A consumer re-deriving the accession from the symbol would sometimes pick the other one and would never be told.

**25. stage4.py:518** - Each decision records its position in the pool as Stage 4 ordered it.

> *Without it:* Without it the ordering has to be trusted across the persistence round trip instead of checked; a reordered artifact reads as correct.

**26. stage4.py:522** - The routing architecture and reason are carried on every decision, including the ones that route nowhere.

> *Without it:* A target that would have gone to an unbuilt architecture is a different finding from one that goes nowhere at all, and the reason string is the only thing that distinguishes them; collapsing them understates what the missing architectures are worth.

**27. stage4.py:532** - Span percentiles use 5 buckets (quintiles).

> *Without it:* Fewer fails to separate a 7 kb gene from a 1.1 Mb one, so the percentile still carries the span confound it exists to remove; more slices the pool so finely that a bucket holds too few pairs to rank within.

**28. stage4.py:537** - The minimum tumour-expression threshold is applied to the partner only, not to the target — the asymmetry is deliberate.

> *Without it:* A target earns its place through the tumour-side composite, which already scores expression and prevalence, while a partner is chosen purely on risk and would otherwise be rewarded for being absent. An AND gate fires only where both antigens are present, so a partner absent from the tumour contributes nothing to killing it while contributing everything to the pair looking safe.

**29. stage4.py:545** - Partner eligibility is measured on bulk tumour transcript level, not on the per-cell measure.

> *Without it:* The per-cell measure is confounded with genomic span (§6.5b); using it here would reintroduce the exact artefact the coverage-does-not-gate decision exists to work around, so partners would be admitted on gene length again.

**30. stage4.py:549** - PARTNER_MIN_TUMOUR_TPM is 5.0 rather than 3.0, and the reason is concentration of partner choice, not safety margin.

> *Without it:* One protein sat far below every other candidate at 0.0277 against 0.2272 for the next lowest. At a 3 TPM threshold the lowest eligible partner still leads the next by 0.0489, so a single gene wins as partner for nearly every target; at 5 TPM the leaders cluster within 0.0036 of each other and no single gene can win for every target. The threshold sits at roughly the pool's 8th percentile and retains 182 of 200, so it is not fitted to exclude two named genes.

**31. stage4.py:628** - A gene with no tumour-expression measurement is ineligible as a partner; missing is not treated as a pass.

> *Without it:* Absence of evidence that the partner is on the tumour is exactly the case the threshold exists for; treating missing as a pass puts the missing-is-a-third-state rule the wrong way round and admits the least characterised partners preferentially. Passing None disables the filter entirely, which is for measuring what the filter does, not for running without it.

**32. stage4.py:644** - When no tolerances are supplied, routing is disabled and every decision reports NOT_CONFIGURED — no default ceiling is fabricated.

> *Without it:* `Ranked` carries no ceiling, so inventing one here would have set it to 0.0 and quietly routed nothing to CONVENTIONAL. A caller that has not been updated keeps its old answers instead of getting a full sweep of false negatives.

**33. stage4.py:664** - Admissible pairs are ordered by how far under the ceiling the combined risk sits (then partner name for determinism), not by co-expression.

> *Without it:* The previous ordering was by f_ab — the better question — but f_ab is confounded with genomic span, so the 'best' partner was substantially the one with the longest gene. This is a weaker criterion honestly measured rather than a stronger one measured on an artefact.

**34. stage4.py:571** - Pairs whose members have no span on record get span_geomean_kb and span_percentile left as None rather than assigned a middle value.

> *Without it:* A default mid-percentile reads as 'typical overlap for genes this long', which is a measurement the pair does not have; the percentile is reporting only and must not manufacture reassurance.

**35. stage4.py:604** - Ranks within a span bucket average across ties instead of using positional order.

> *Without it:* A block of identical f_AB values would otherwise be given an ordering the data does not support, and the resulting percentiles would separate pairs that are indistinguishable.

**36. stage4.py:732** - The ADAPTOR outcome is reached only after the target has cleared neither alone nor paired, and the risk number is carried through unchanged — only the ceiling it is compared against differs.

> *Without it:* This is the row that recovers targets the old blind gate discarded. Rewriting the risk instead of the ceiling would make the recovery look like a safety improvement when the antigen is no safer, only the exposure stoppable.

**37. stage4.py:749** - The `coverage_below_floor` counter counts pairs that cleared risk and were measured but fell under the reported floor, and is kept separate from `unmeasured`.

> *Without it:* The floor no longer excludes anything, so a non-zero count here is information about the pair, not a reason it was rejected — folded into `unmeasured` it reads as a rejection reason that no longer exists. A pair nobody could measure and a pair measured and found thin are different facts.

**38. stage4.py:764** - A dedicated `partner_ineligible` counter records pairs rejected only because the partner is not expressed enough in the tumour.

> *Without it:* Without this row a target excluded entirely on partner eligibility is persisted with every counter at zero — a rejection with no stated reason, indistinguishable from a bug.

**39. stage4.py:771** - The UNRESOLVED-versus-NO_DESIGN test uses `coverage.measured`, not `coverage_ok`.

> *Without it:* The coverage floor no longer selects anywhere else in the stage (§6.5b); leaving it here decides NO_DESIGN against UNRESOLVED on a threshold the stage has stopped applying, so the outcome disagrees with the admissibility rule that produced it.

**40. stage4.py:775** - Partner eligibility is applied to the salvageable set as well, not only to the admissible set.

> *Without it:* A target whose only salvageable pairs run through partners the tumour-expression gate rejects is not salvageable: resolving the missing organ would not make that partner usable, so reporting UNRESOLVED promises a design that cannot follow.

**41. stage4.py:835** - The configuration hash covers the SELECTION_RULE string and the partner TPM threshold, not just the numeric floors.

> *Without it:* Without it a run from before coverage was removed from selection hashes identically to one after, and `read_decisions(expect_stage4_hash=...)` accepts the old artifact as current — the one thing carrying the hash is meant to stop.

**42. stage4.py:833** - Retired parameters (coverage floor, patient fraction floor) stay in the configuration hash even though neither selects any more.

> *Without it:* They still bound the reported coverage numbers, so a run with different floors produces different reported numbers and must not hash as the same experiment.

**43. stage4.py:850** - Decisions are persisted payload-first, manifest-second, both atomic, with the manifest acting as the commit marker.

> *Without it:* Nothing downstream could otherwise read the decisions without re-running the stage, and re-running re-derives the very numbers under question. A payload with no manifest beside it is a run that died mid-write; without the commit-marker discipline it reads as a complete result.

**44. stage4.py:860** - The manifest carries a version that is bumped when the payload shape changes, and the reader refuses a mismatch.

> *Without it:* An artifact written under an older layout otherwise reads as current and fails somewhere further away, where the cause is no longer visible.

**45. stage4.py:872** - `failed_on` is written as None when empty rather than as an empty mapping.

> *Without it:* An absent mapping and a mapping of zeros mean different things — no terminal outcome versus a terminal outcome with no counted cause — and collapsing them hides the second case.

**46. stage4.py:896** - The criteria outcomes are written into the manifest and `write_decisions` refuses to write without them.

> *Without it:* These decisions are currently produced by a run that stops on five tripped criteria; an artifact that did not say so would be read as a result. The consumer is expected to refuse the payload when anything is tripped, which is why the outcome is stored beside the data and not in a log.

**47. stage4.py:908** - The criteria argument is validated before anything on disk is touched, ahead of the manifest unlink.

> *Without it:* Raising after the unlink would destroy a previously valid artifact to punish a bad call, turning a caller's mistake into data loss.

**48. stage4.py:918** - The old manifest is deleted before the new payload is written, not after.

> *Without it:* A crash between the two writes then leaves an unblessed payload, which the reader refuses. In the other order it leaves a stale manifest blessing a new payload — a blessed mismatch the reader cannot detect.

**49. stage4.py:959** - `read_decisions` returns the manifest alongside the rows, and refuses unusable artifacts unless `allow_unusable=True` is passed explicitly.

> *Without it:* The rows alone do not say whether they may be read as a result and a caller handed only the rows cannot find out; today they may not, because the writing run stops on five tripped criteria. Wanting them anyway — to inspect why the run stopped — is legitimate, so it has to be said out loud.

**50. stage4.py:965** - The payload digest is re-derived from the rows on read rather than trusted from the manifest.

> *Without it:* A truncated or hand-edited payload that still parses as JSON is exactly the failure this guards, and it is completely silent without the check.

**51. stage4.py:999** - The recorded stage3/stage4 hashes are actually compared against the caller's expectations, not merely stored.

> *Without it:* Recording them and never checking them lets an artifact from a different ranking be read as current — precisely what carrying the hashes was supposed to prevent.

**52. stage4.py:125** - Combined risk per organ is min(a, b), the scarcer of the two members, and only then the maximum across organs.

> *Without it:* Taking the maximum per organ reduces the pair to its more dangerous member and ignores the architecture entirely - which is the whole point of pairing. A gate that needs both antigens present is bounded in each organ by whichever is scarcer there. The minimum assumes perfect overlap within the organ: pessimistic for safety and optimistic for coverage at the same time, which is what makes it the right thing to gate on. It is a bound, not a measurement - neither source carries a joint distribution over cells, so within-organ co-expression cannot be measured at all. The max across organs comes after, because the patient has every organ.

**53. stage4.py:128** - An organ where one member is measured and the other is not is recorded as unresolved, and the unmeasured member is assumed present rather than absent.

> *Without it:* The pair's entire claim is an absence, so crediting it with an absence nobody observed lets an unmeasured organ manufacture safety. Organs measured for neither member stay outside both mappings and contribute to nothing, which preserves the identity that pairing a target with itself reproduces its own risk.

**54. stage4.py:141** - Pair risk is rounded to the same precision the ranking stage stores single-target risk at.

> *Without it:* Comparing a full-precision pair risk against a rounded single risk reports a disagreement at the fifth decimal that is arithmetic rather than substance, and a reviewer cannot tell that from a real one.

**55. stage4.py:778** - `write_decisions` takes the tolerances the run used and threads them through the configuration hash, and records them in the manifest as a structure rather than leaving them to be inferred from a reason string.

> *Without it:* Decision 13 against `routing.py` puts the declared ceilings in the hash so that two runs under different policy cannot compare as one. `write_decisions` computed the hash without them, which was harmless only while the artifact it wrote was never routed. The moment one is, a routed and an unrouted set hash identically and each reads as the other's cache. The manifest entry is what a criterion can test structurally; matching on `route_reason` text would make the check a string comparison against a message.

**56. stage4.py:798** - The manifest tallies every outcome the stage can emit and every outcome actually present, not a fixed list of four.

> *Without it:* The list predated the adaptor outcome and did not include it, so a routed set of 200 rows reported as 195. A tally whose categories are enumerated separately from the thing being tallied drops whatever is added later, silently and while still summing to a plausible number.


## `car_pipeline/stages/routing.py`

**1. routing.py:4** - The risk profile selects the architecture first and the ceiling follows from it, instead of one risk ceiling applied at Stage 3 before any architecture is known.

> *Without it:* The old order is backwards: the architectures that exist to make a risky target tolerable were only ever offered to targets that had already cleared without them. 199 of 200 pool members died at that gate.

**2. routing.py:10** - Two ceilings, persistent and terminable, are kept strictly separate and never blended into one relaxed number.

> *Without it:* An adaptor does not make the antigen safer — the adaptor still binds it — it makes the exposure terminable, because activation needs a separately dosed protein with a finite half-life. Magnitude and reversibility are different axes; a blended ceiling reports a target as safer when only its stoppability changed.

**3. routing.py:17** - Routing never rewrites a risk number — the target's Stage 3 risk is carried through unchanged and only the ceiling it is compared against changes.

> *Without it:* Substituting the receptor's risk for the target's makes the gate vacuous, because every adaptor receptor binds the same tag and so looks equally harmless; every target would clear.

**4. routing.py:27** - Architectures are tried in a fixed order of increasing product complexity — one receptor and one product, then two receptors, then two manufactured products — and the first admitting architecture wins.

> *Without it:* This is a stated preference, not a tuned score. Scoring architectures instead would let a two-product design win over a conventional one that was already admissible, adding manufacturing complexity for no risk benefit.

**5. routing.py:36** - Architectures the spec names but this stage does not build (AND_NOT, SWITCHABLE) are reported by name with a reason instead of being dropped.

> *Without it:* A target that would have routed to an unbuilt row is a different finding from one that routes nowhere; collapsing the two understates exactly what the missing architectures are worth.

**6. routing.py:49** - AND_NOT is not built on the pairing output: the pair source cannot supply an exclusion antigen.

> *Without it:* Pairing selects for tumour co-expression and an inhibitory CAR needs the opposite relation — reusing the pair partners as exclusion antigens would pick precisely the antigens most present on the tumour, inhibiting the design where it should fire.

**7. routing.py:51** - SWITCHABLE is not built because the FKBP12 in this build is wild-type.

> *Without it:* A rapamycin ON-switch and the mandatory rimiducid suicide switch would then answer to the same drug, so turning the product on and killing it are the same action. The fix is a point mutation, which has no provenance class in this pipeline.

**8. routing.py:79** - `Tolerances.terminable` is optional with no default; its absence disables the adaptor row rather than being filled in.

> *Without it:* Both ceilings are policy inputs, not measurements — the platform cannot derive a clinical risk tolerance from expression data. A default would be this code quietly setting clinical policy, and every adaptor routing would inherit a number nobody declared.

**9. routing.py:114** - A target with no measured risk routes to NO_ARCHITECTURE with a stated reason rather than being compared against any ceiling.

> *Without it:* Unmeasured is not safe. A None risk compared against a ceiling either errors or passes, and passing routes an uncharacterised antigen to a conventional persistent design.

**10. routing.py:125** - The AND gate is admitted against the same persistent ceiling as a conventional design, not a relaxed one.

> *Without it:* An AND gate makes activation conditional; it does not make the exposure terminable — the T cell still self-amplifies and cannot be withdrawn. Moving the ceiling for it would grant reversibility credit the architecture does not provide.

**11. routing.py:136** - When no terminable ceiling is declared, the target is reported as NOT_CONFIGURED rather than NO_ARCHITECTURE — the deferral is evaluated after the unbuilt-row checks.

> *Without it:* 'No ceiling was declared for this project' and 'the risk exceeds every declared ceiling' are different findings; merging them reports a policy gap as a scientific rejection and the missing declaration is never noticed.

**12. routing.py:161** - `sweep` reports how many targets the adaptor row admits across a whole range of terminable ceilings instead of defending one value.

> *Without it:* The terminable ceiling is a number this pipeline cannot measure (criterion A9), so a single chosen value asks the reader to trust it; the sweep shows what any choice would have bought.

**13. routing.py:177** - The declared ceilings are folded into the Stage 4 configuration hash via `configuration_payload`.

> *Without it:* Two runs under different tolerances route targets to different architectures and are different experiments; without the ceilings in the hash they compare as one and a cached artifact from other policy reads as current.


## `car_pipeline/stages/stage5.py`

**1. stage5.py:8** - The structure route and the sequence route are reported apart and never summed or collapsed into one binder count, and a target with a sequence-route binder but no deposited structure is not NO_BINDER.

> *Without it:* Collapsing the two discards the variable-region sequences, which are the more useful half — the thing a construct is actually built from — and would mark a target that has a named therapeutic with heavy and light sequences as having no binder at all.

**2. stage5.py:11** - A deposited complex is confirmed as antibody-containing by curated chain annotation, not by reading the entity description text.

> *Without it:* Deciding from free-text entity descriptions makes the antibody call on prose rather than on the annotated heavy/light chain assignment, so entries are admitted or rejected on wording; the curated annotation is also what supplies the heavy/light/antigen chain identifiers the candidate is built from.

**3. stage5.py:20** - Affinity is the literal string NOT_CONNECTED for every candidate and is never a number, and nothing is ranked on it.

> *Without it:* The reason is measured, not assumed: the curated structure summary carried affinity in a previous release and does not in this one, and the bioactivity database returns zero records for both validation molecules and both validation targets. 'We did not rank on affinity' and 'we could not' are different statements and only the second is true; printing a number, or ranking as though one existed, asserts the first.

**4. stage5.py:27** - The stage does not re-rank targets and carries Stage 4's emission order unchanged; binder availability never reorders anything.

> *Without it:* Binder availability tracks how much attention a protein has had, so ordering by it lets the literature choose the targets — exactly what the discovery stages exist to avoid. The output would look like a ranked shortlist while encoding popularity.

**5. stage5.py:72** - Every candidate's isoform is reported as ISOFORM_UNRESOLVED rather than defaulting to the canonical isoform.

> *Without it:* Neither route can say which isoform a binder engages: the therapeutic table records only the gene, and the structure route would need the deposited antigen sequence aligned against each isoform. Filling in the canonical isoform prints a specific, checkable-looking claim that was never determined.

**6. stage5.py:97** - Entries that exist for the accession but carry no antibody instance are counted into entries_without_antibody instead of being silently skipped.

> *Without it:* Without the count, 'the protein has structures' and 'the protein has a binder' become indistinguishable downstream — the gap between the two claims is the point, and a silent skip erases it.

**7. stage5.py:175** - Entries whose coordinates are computed rather than measured are excluded from the structure route, and counted, even when they do contain an antibody.

> *Without it:* A computed model is not retrieved evidence; admitting one reports a binder that was predicted, not observed. Counting the exclusion is what keeps the structures-versus-binders gap reported downstream exact.

**8. stage5.py:101** - Model exclusions are kept in a separate counter from no-antibody entries rather than merged into one 'entries dropped' number.

> *Without it:* 'No antibody in it' and 'not an experiment' are different reasons to discard an entry; one combined count cannot be read back as either, and the reader cannot tell whether a target lacks binders or merely lacks experimental depositions.

**9. stage5.py:214** - The retrieval is cached to disk as payload-then-manifest, with any existing manifest deleted before the payload is rewritten, rather than re-running retrieval for each downstream stage.

> *Without it:* Retrieval is one network call per pool member and about five minutes; re-running it per stage is slow, loses everything to a single dropped connection, and hits the external service again for nothing. Writing the payload while the old manifest is still in place would leave a stale manifest blessing a half-written payload if the run dies mid-write.

**10. stage5.py:283** - Cache reuse requires the Stage 4 configuration hash to match as well as the gene set; stage4_hash is mandatory and omitting it means 'retrieve', not 'trust the disk'.

> *Without it:* Changing a criticality override or the risk ceiling changes the Stage 4 configuration hash while leaving the top-200 pool identical, so a gene-set check alone hands back binders retrieved under the previous configuration and records the wrong provenance in every artifact that chains from this one — a silent error with no visible symptom.

**11. stage5.py:294** - Only the 'no manifest' CacheError falls back to retrieval; any other cache integrity failure is re-raised rather than caught and retried.

> *Without it:* A pool or configuration change is ordinary and should just retrieve, but a truncated or digest-mismatched payload is an operator's problem: swallowing it turns it into an unexplained five-minute pause and leaves the bad artifact on disk for the next stage to crash on.

**12. stage5.py:310** - A payload that passes its digest but cannot be rebuilt into records falls back to retrieval with a printed message, not silently.

> *Without it:* That failure means the manifest version was not bumped when the record shape changed. Retrying quietly hides a layout-versioning bug that will keep invalidating every cache write, and the operator sees only unexplained slowness.

**13. stage5.py:348** - Read_binders checks the manifest layout version in addition to the SHA-256 digest, and refuses a payload no manifest blesses.

> *Without it:* A payload written under an older record layout is internally consistent and passes the digest check perfectly, so digest-only validation reconstructs records with missing or misinterpreted fields; a payload with no manifest beside it came from a run that never finished.


## `car_pipeline/stages/stage6.py`

**1. stage6.py:32** - The DNA is reverse-translated with one fixed codon per amino acid and is explicitly a map, not an ordering sequence — it is not codon-optimised.

> *Without it:* Codon optimisation makes the round trip non-deterministic and destroys the exact domain boundaries the segment map depends on; worse, an optimised sequence reads as something you could order and synthesise, which this is not. The table's human-frequent choices are incidental — reproducibility is the point.

**2. stage6.py:55** - Residues the codon table cannot encode are returned by assemblable() and hard-fail the assembly, rather than being mapped to a fallback codon.

> *Without it:* One retrieved therapeutic carries lowercase residues, and the obvious `CODON.get(residue, "NNN")` turned those into an ambiguous codon that translated back as a mismatch — a construct whose DNA did not encode its own protein, produced with no error raised.

**3. stage6.py:195** - Best_binder picks the shortest usable sequence-route binder (and filters out any whose residues cannot be encoded), not the first or the most clinically advanced.

> *Without it:* The budget verdict is only meaningful as a bound: the smallest binder that fits is the most favourable reading, so BUDGET_EXCEEDED means no available binder fits. Picking an arbitrary or most-advanced binder reports an over-budget construct when a shorter binder would have fit, and admitting an unencodable candidate would fail assembly later.

**4. stage6.py:44** - Constructs are built only for outcomes Stage 4 recommended (SINGLE, DUAL, ADAPTOR); a target with a usable binder but a rejected outcome gets NO_CONSTRUCT with the outcome named in the reason.

> *Without it:* Assembling for a rejected target presents a finished design for something upstream says is not designable, and a reader cannot be expected to carry that caveat alongside a printed sequence.

**5. stage6.py:224** - ADAPTOR is included in the buildable outcomes despite having no supplied binder sequence, and the adaptor branch runs before the missing-binder guard.

> *Without it:* It is a real routed architecture that assembles; excluding it, or letting it fall through to the 'no binder' guard, drops a valid design because its receptor binds the tag rather than the antigen — the antigen specificity lives in a separately manufactured adaptor that is not in this vector, which is the whole reason it fits the budget.

**6. stage6.py:165** - When a part declares a length but no residues, the assembler returns an empty protein and empty DNA while still emitting the full segment map, instead of padding the missing residues with filler.

> *Without it:* The layout and the length are real but the residues are not supplied and must not be invented — a padded protein would be read downstream as a designed sequence. For the adaptor receptor this is the only such part: no anti-tag antibody exists in the cached structural set.

**7. stage6.py:143** - Segment coordinates still advance by a sequence-less part's declared size, and the construct's length comes from declared_bp (segment spans plus the stop codon) when there is no DNA to measure.

> *Without it:* If the unsupplied part contributed zero, every downstream segment boundary would shift and the construct would measure short against the 3500 bp budget — a plausible-looking length that under-counts the receptor and could turn a BUDGET_EXCEEDED into a BUILDABLE.

**8. stage6.py:288** - A dual design splits the signal across two receptors — activation (CD3zeta) on one, costimulation (4-1BB) on the other — rather than putting both binders on one fully signalling receptor.

> *Without it:* With both scFvs on one complete receptor, either antigen alone gives a full signal and the design is no longer an AND gate; the split means neither antigen alone produces a complete signal.

**9. stage6.py:321** - Partner_binder_name is filled only when the outcome is DUAL, not whenever row['partner'] is present.

> *Without it:* Stage 4 fills `partner` on the SINGLE branch too, so reading it unconditionally labels a single-arm construct with a second binder name that its amino acid sequence does not contain — a wrong annotation that looks entirely plausible.

**10. stage6.py:256** - The missing-binder verdict says 'no binder with a usable variable-region sequence' rather than 'no binder'.

> *Without it:* Assembly needs a sequence so it reads the sequence route only; a target with a structure-route binder and no sequence lands here, and calling that 'no binder' directly contradicts Stage 5, which reported one.

**11. stage6.py:336** - This stage's configuration hash covers Stage 5's hash, not Stage 4's, and the 3500 bp budget is carried from Stage 1 rather than recomputed.

> *Without it:* Passing the upstream-of-upstream hash leaves a Stage 5 change invisible to this stage's identity, which is the whole reason the hash is carried — a stale construct set would pass as current. Re-deriving the budget locally lets it drift from Stage 1's without anything reporting a conflict.


## `car_pipeline/stages/construct_safety.py`

**1. construct_safety.py:1** - Every finding carries a basis, CODON_INVARIANT or MAP_SPECIFIC, and the two are reported side by side rather than merged.

> *Without it:* Stage 6's DNA is a reverse translation under one fixed codon per residue, so every leucine is CTG and every serine AGC. A nucleotide finding on that map is a property of the encoding and not of anything anyone would order: repeats are inflated because identical peptides produce identical nucleotides, and whether GT falls at a position is decided by the codon table. Reporting the four flat would state properties of an arbitrary encoding as properties of a therapeutic, and every number would look like real sequence analysis.

**2. construct_safety.py:66** - The repeated-part detector reads the domain map, not the sequence.

> *Without it:* The dominant recombination hazard in these constructs is architectural: a dual design carries the CD8A leader, hinge and transmembrane segment twice by construction. Read from the segments it is exact and survives any encoding; read from nucleotides it would be one more map-specific string match, and a codon-optimised sequence would hide it by diversifying the copies.

**3. construct_safety.py:243** - The arm gates nothing and mutates nothing.

> *Without it:* Turning a count of splice motifs into a block sets a tolerance, and this platform has no outcome data to set one from. Rewriting a codon to remove a site would make Stage 6 emit a sequence that is no longer the map its own specification says it is.

**4. construct_safety.py:20** - Every threshold is a named constant fixed in the specification before the run, and stated as conventional rather than measured.

> *Without it:* A threshold chosen after seeing the counts is the defect this repository has recorded twice, and there is nothing here to calibrate against. The criteria test the detector against planted known answers in both directions; none of them tests a threshold.

## `car_pipeline/stages/stage9.py`

**1. stage9.py:5** - Off-tumour risk is carried verbatim from Stage 3 and never recomputed here, even though the gate is the place that consumes it.

> *Without it:* A second implementation of the tissue mapping would be "a second place for its bugs to live" — the gate could block or admit on a risk number that disagrees with the one Stage 3 published for the same target.

**2. stage9.py:8** - The passing verdict is named PASSES_STATED_CHECKS rather than SAFE.

> *Without it:* Passing means only that three specific questions failed to show a problem; a verdict called SAFE reads as a safety claim the stage never made, and everything it does not check (epitope immunogenicity is not even connected) silently becomes part of the claim.

**3. stage9.py:28** - ORIGIN_STEMS is ordered longest stem first and matched in that order, not in dictionary/alphabetical order.

> *Without it:* `-xizu-` ends in `-zu` and contains `-xi-`; a shorter-first scan calls a chimeric/humanised binder plain "chimeric" or "humanised". Every stem in the table is also a suffix of another one (`o`, `u` match almost anything), so order is the whole rule.

**4. stage9.py:46** - An unrecognised or absent stem returns ORIGIN_UNKNOWN — origin is never guessed from sequence and never defaulted to human.

> *Without it:* The stem is "a CONVENTION, not a measurement": a molecule re-engineered after naming keeps its original stem, and a binder with no INN has no stem at all. Defaulting the no-stem case to "human" turns an unanswered question into a clean immunogenicity result; the flag text itself has to say the finding "is a convention and not a measurement".

**5. stage9.py:123** - Origin is read from every named binder and the distinct set is kept, instead of taking one binder (e.g. the alphabetically first) as the target's binder.

> *Without it:* Measured: picking one "would have read CLDN18 as ORIGIN_UNKNOWN from Ciletatug while Zolbetuximab, the spec's own example, sits in the same set and is chimeric" — the immunogenicity flag vanishes on the spec's own worked case.

**6. stage9.py:74** - SafetyRecord keeps a `binder_origins` list of every distinct origin, not a single origin field.

> *Without it:* "A target with both a human and a murine binder is not described by either alone" — collapsing to one value reports a mixed binder set as uniformly safe or uniformly foreign, and the reader cannot see which binders drove the verdict.

**7. stage9.py:131** - The single reported `binder_origin` is the foreign one when any binder is foreign, rather than the first or most common origin.

> *Without it:* The conservative reading is used because the design may use the foreign binder; taking the first origin lets a target with one human and one murine binder report "human" and drop the anti-CAR immunogenicity flag entirely.

**8. stage9.py:77** - `epitope_immunogenicity` is hard-set to NOT_CONNECTED instead of being computed or defaulted.

> *Without it:* "Not a lookup" — answering it needs a k-mer scan of the variable region against the bulk epitope table, which this stage does not do. Any default value (empty, none-found, low) reads as a negative result from a check that was never run.

**9. stage9.py:83** - Stopped-trial tallies carry a `trials_truncated` flag and the reason text appends "(a floor: tallied over one page)" rather than reporting the count as a total.

> *Without it:* The tallies can cover fewer studies than the registry holds, so "2 trials terminated" is a lower bound presented as a complete count; a reader comparing targets would treat an under-counted target as the cleaner one.

**10. stage9.py:101** - A binder is anything usable from either route — named sequence candidates plus structure candidates with an identifier — not sequences only.

> *Without it:* Stage 5 defines a structure-route binder as a binder; counting only sequences here contradicts it, so "a target with a solved complex and no named therapeutic would read as ungateable" (NO_GATE) instead of being gated on the binder it actually has.

**11. stage9.py:136** - A missing Stage 3 risk is BLOCKED, not skipped, defaulted to zero, or passed through.

> *Without it:* Stage 3 §"Undefined risk is not low risk" is explicit; treating a missing measurement as a pass "would invert it" — a target nobody measured would emerge with the same verdict as one measured clean.

**12. stage9.py:147** - The risk comparison uses the ceiling the target was actually routed against (`route_ceiling`, e.g. the terminable one), not the persistent ceiling passed into the gate.

> *Without it:* Gating every target on the persistent ceiling "would re-apply the blind gate that routing exists to replace, and would block an adaptor design for the very risk its architecture was chosen to carry" — Stage 4a's architecture choice would be silently overruled at the gate.

**13. stage9.py:164** - Admission under a terminable ceiling appends an explicit reason (which also makes the verdict FLAGGED) instead of passing silently.

> *Without it:* "Recorded, not silent": clearing a stoppable-exposure ceiling is a different claim from clearing the persistent ceiling of 0.15, and without the note a target admitted at the looser tolerance is indistinguishable in the output from one that cleared the strict one.

**14. stage9.py:174** - A target with no usable binder gets NO_GATE with "no binder, so there is nothing to gate", not PASSES.

> *Without it:* With no binder there is nothing to check for origin or immunogenicity, so the checks trivially raise no reason; falling through to PASSES_STATED_CHECKS reports checks that were never actually run as having been passed.

**15. stage9.py:204** - The configuration hash folds in the ceiling, the origin stem table, and the NOT_CONNECTED epitope marker, not just the upstream hash and gene list.

> *Without it:* "The ceiling decides every BLOCKED verdict, so a run at a different tolerance must not hash the same as this one" — otherwise two runs with opposite block/pass verdicts are indistinguishable by hash and a cached result is reused across tolerances.


## `car_pipeline/stages/stage10.py`

**1. stage10.py:6** - Developability scores the Stage 5 binder variable regions, not the Stage 6 constructs that are the actual product.

> *Without it:* "Stage 6 produces zero buildable constructs and scoring those would score nothing" — the obvious target of the analysis yields an empty report every run.

**2. stage10.py:9** - Liabilities are reported as separate counts and values and are never summed into a single developability score.

> *Without it:* "One number would let a strong liability be averaged away by four weak absences" — a binder with an odd cysteine and three aggregation-prone regions scores the same as a mildly imperfect one, and nothing here predicts manufacturing failure anyway.

**3. stage10.py:23** - One fixed pKa table (C 8.5, D 3.9, E 4.1, H 6.0, K 10.5, R 12.5, Y 10.1, N-term 9.7, C-term 2.3) is written down in the module rather than taken from a library default.

> *Without it:* "The exact values differ between published sets", so an implicit table makes the isoelectric point non-reproducible: the same sequence yields a different pI across environments, and the pI-near-formulation-pH flag flips with it.

**4. stage10.py:61** - PI is found by 100-step bisection bounded to pH 1..14 and always returns a value, rather than being undefined for sequences with no ionisable side chain.

> *Without it:* "A sequence with no ionisable side chain still has termini, so the point always exists inside those bounds" — returning None/0 for such a sequence would drop it from the pI-window flag instead of scoring it.

**5. stage10.py:77** - The N-glycosylation sequon scan excludes N-P-S/T — proline at the X position is not a sequon.

> *Without it:* "The proline exclusion is the part of this rule most easily dropped, and dropping it silently inflates every count" — the sequon flag fires on motifs that are not glycosylated, and the inflated count looks like an ordinary result.

**6. stage10.py:100** - Aggregation-prone hits are merged: a new start is only recorded once it clears the previous start by the full 7-residue window, instead of one hit per qualifying window.

> *Without it:* A single hydrophobic patch satisfies the mean-hydropathy>=1.0 test at many consecutive offsets, so counting every window reports one region as several and the "N region(s)" figure scales with patch length rather than with the number of patches.

**7. stage10.py:143** - Flags are (kind, detail) pairs with the kind as a stable key, not a single prose sentence.

> *Without it:* "The kind is a stable key so the report can group without parsing prose back out of a sentence" — grouping by sentence text breaks the moment a threshold or number inside the message changes.

**8. stage10.py:169** - Cysteines are reported as parity, "odd" or "even" — never as an unpaired count, and never "unpaired: 0".

> *Without it:* "Pairing cannot be read from sequence, and an even count is only the absence of a guarantee, not evidence of pairing" — reporting zero unpaired cysteines states a disulfide assignment the sequence does not support.

**9. stage10.py:179** - When no binder carries a sequence, `assess` returns the NOTHING_TO_SCORE status, not just an empty row list.

> *Without it:* "A table with no rows reads as 'nothing had liabilities', which is the opposite of 'nothing was examined'" — the empty table is a clean bill of health for binders that were never scored.


## `car_pipeline/stages/stage11.py`

**1. stage11.py:6** - Ranking is a Pareto front over four objectives; there is no weighted sum and no invented weights.

> *Without it:* "The objectives are not commensurable": a weighted total "would let a design with an unmanageable safety margin be rescued by a strong tumour score" — precisely the trade Stage 3 refuses by keeping its three numbers apart.

**2. stage11.py:28** - Every pool member is attributed to the first gate it fails, in pipeline order, rather than counted against each gate it fails.

> *Without it:* "The counts partition the pool rather than overlapping" — independent counting makes the attrition figures sum to more than the pool size, and a target blocked on tissue risk is also counted as having no binder and no construct.

**3. stage11.py:38** - ADAPTOR is included in RECOMMENDED alongside SINGLE and DUAL.

> *Without it:* Measured: leaving it out "attributed eight routed designs to 'no design recommended' — a gate they never failed". ADAPTOR is a real recommendation made against the terminable ceiling, so its omission moved eight designs into an attrition bucket for a failure that did not happen.

**4. stage11.py:23** - Survivors that all lack a binder sequence get their own status, RANKED_AWAITING_BINDER, instead of RANKED or NO_DESIGN_REACHES_THE_END.

> *Without it:* "Nothing here can be ordered, and a plain RANKED would promise otherwise"; NO_DESIGN_REACHES_THE_END would equally misreport, because designs did get through. Both obvious statuses give a caller a false answer to a different question.

**5. stage11.py:55** - `binder_supplied` is a separate field rather than folded into `survived`.

> *Without it:* "'Nothing got through' and 'designs got through with one part outstanding' are different answers and a caller has to be able to tell them apart" — folding the missing sequence into survival deletes designs that passed every gate from the survivor count.

**6. stage11.py:61** - Cleanliness is the negated liability flag count, so higher is better like every other objective.

> *Without it:* Pareto dominance tests `x >= y` on all four objectives; feeding raw flag counts in makes the most liability-laden binder dominate the cleanest one, and the front is populated by exactly the wrong designs while every number still looks plausible.

**7. stage11.py:125** - ADAPTOR outcomes are exempted from the "no binder retrieved" gate rather than tested for a binder like every other architecture.

> *Without it:* "An adaptor receptor binds a tag, not the antigen, so it needs no target binder. Failing it here would attribute it to a gate its architecture does not have to pass" — valid adaptor designs disappear into an attrition bucket.

**8. stage11.py:131** - "Assembled" is tested as having segments (a layout and a length), not as carrying residues.

> *Without it:* "An adaptor's binder sequence is declared unsupplied rather than absent, and calling that 'no construct assembled' would report a design that exists as one that does not" — the design is counted as an assembly failure instead of a survivor awaiting a binder.

**9. stage11.py:164** - The configuration hash includes the RECOMMENDED tuple, not just the upstream hash, genes and gate names.

> *Without it:* "Which outcomes count as a recommendation decides the whole attribution, so a run that admits ADAPTOR must not hash as one that does not" — otherwise the run that misattributed eight designs is hash-identical to the corrected one.


## `car_pipeline/stages/stage12.py`

**1. stage12.py:1** - The package assembles and never recomputes. Every number in it is carried from the stage that measured it.

> *Without it:* An assembler that re-derives anything becomes a second opinion on a number the platform already has, and the two can disagree without either being wrong. Carrying means the only defect the stage can have is a lossy copy, which is what criteria Q4 and Q5 exist to catch: they re-assert the round trip and the risk reconstruction on the packaged copy rather than trusting the source.

**2. stage12.py:51** - Every gap the package declares carries a probe the verifier executes, and a gap whose claim is a judgement carries none and says so.

> *Without it:* The gaps section is the one place the package asserts something of its own, and an assertion nobody recomputes is exactly how the construct narrative came to say three single-antigen targets were one. A probe makes the claim falsifiable: build Stage 7 and forget the table, and Q6 trips. Marking the nine judgements as unprobed keeps the split between verified and asserted visible instead of implying the whole table was checked.

**3. stage12.py:236** - A gap that is a property of this run rather than of the code is recomputed from the run, not declared in the table.

> *Without it:* Stage 10 covering none of the shipping designs is true today because all five survivors are adaptors. Writing that into a static table would make it a fact about one run that nobody recomputes, and it would survive the day it stopped being true. Recomputed, it disappears by itself the moment one shipping design carries a sequence-route binder.

**4. stage12.py:257** - No section, key or placeholder is emitted for Stage 7 or Stage 8; their absence lives only in the gaps section.

> *Without it:* A `structural_report: null` reads as computed-and-found-nothing, which is the strongest available claim about a stage that does not exist. Absence stated once with its reason is a different sentence from a null field, and Q9 enforces the difference.

**5. stage12.py:394** - Where Stage 5 retrieved nothing but the construct carries a binding domain, the binder section names the route that supplied it.

> *Without it:* Every surviving design in the worked indication is an adaptor, so its binder section reads `NO_BINDER` beside a construct that plainly has a binding domain. The two are consistent - no antigen-specific binder was retrieved, and the receptor binds a tag instead - but a reader seeing only the verdict would conclude the construct has no binder, which is false.

**6. stage12.py:449** - Where Stage 10 scored nothing for a candidate, the developability section says why rather than showing an empty table.

> *Without it:* Stage 10 reads Stage 5 sequence-route binders. An adaptor's binding domain comes from a deposited structure, so no shipping design is scored, and an empty table reads as a clean sheet. It is the opposite: no developability figure in this platform describes the binder these constructs carry.

---

# The service

## `car_pipeline/api/pipeline.py`

**1. pipeline.py:33** - An indication with no single-cell atlas is refused as NOT_USABLE rather than scored on the components that remain and served with the mean renormalised over them.

> *Without it:* Dropping the atlas costs C1 and C2 and still scores 3,399 of 3,466 targets, filling the top of the pool with immunoglobulin, TCR and MHC-II genes, because malignant_vs_stroma (C2) is the ONLY component that rejects stromal and immune expression. Losing 0.45 of weight is survivable arithmetic; losing the only discriminator against stroma is not, and renormalising hides it by rescaling what remains as the whole score. A number that looks like an answer is worse than a refusal, because only one of the two gets checked.

**2. pipeline.py:5** - A run that produces no buildable design returns what each stage measured and is recorded as complete; nothing maps an empty result onto an error.

> *Without it:* Turning emptiness into a failure discards the finding — the emptiness IS the result, and something specific and measured stopped each design.

**3. pipeline.py:71** - A missing tumour cohort or normal-tissue denominator raises immediately instead of letting the loaders fall back to their defaults.

> *Without it:* None means "not connected", but every loader reads it as "use the default", which is the reference indication's pancreas cohort — so the indication would be screened against another one's data and reported under its own name.

**4. pipeline.py:107** - A DepMap lineage that fails to load degrades the run and names the source in `unavailable`, rather than propagating the exception.

> *Without it:* A lineage with no screened cell lines reaches np.vstack and raises "need at least one array to concatenate"; six of the 36 lineages in the cached model table qualify. Escape resistance is 0.05 of the weight, so ending the run costs far more than the missing component does.

**5. pipeline.py:131** - The atlas-less refusal returns before the pairing stage, not at the end of run() where it was originally written.

> *Without it:* Deciding it 55 lines after pairing had already dereferenced the absent atlas made the refusal unreachable: the real observed behaviour was an AttributeError, not the refusal the code appeared to implement.

**6. pipeline.py:180** - Malignant-cell loading in the pairing stage constructs SingleCellSource(indication.atlas) explicitly rather than a bare SingleCellSource().

> *Without it:* The bare constructor silently paired one indication's targets against another indication's single-cell data; the only reason it was ever visible is that the artifact now carries its accession in the filename.

**7. pipeline.py:199** - Both routing tolerances come from the project spec; a project with no declared terminable ceiling gets no adaptor row at all rather than a default one.

> *Without it:* Inventing a terminable ceiling here would route designs to an adaptor architecture against a tolerance no clinical policy ever set, while the run reported it as the project's own.

**8. pipeline.py:190** - Gene-span annotation is wrapped so a failure is swallowed; it annotates and never fails the run.

> *Without it:* Letting the span source raise ends a completed screen over a decoration none of the scores depend on.

**9. pipeline.py:211** - Stage5 binder retrieval is keyed on the Stage 4 hash, not the Stage 3 hash.

> *Without it:* The Stage 4 hash is the slot this artifact records and what every downstream hash chains from; passing the Stage 3 hash corrupts the provenance of a cache three other drivers read.

**10. pipeline.py:282** - Mode A (validate one supplied target) runs the entire Mode B screen and reads the target out of it, instead of taking a cheap direct path to that one protein.

> *Without it:* A separate code path lets the two modes drift, and the first symptom is the platform reporting one risk in discovery and a different one in validation for the same protein on the same evidence. The cost of agreement-by-construction is ranking the whole proteome, about twenty seconds.

**11. pipeline.py:315** - A symbol absent from the reviewed proteome returns UNKNOWN_TARGET, while a symbol present but not surface-accessible returns UNSUITABLE; the two are never collapsed.

> *Without it:* They are different findings: the first is a question about the input, the second is a real rejection on topology (a CAR binder has nothing to engage). Merging them reports a typo as a scientific verdict, or a genuine topology rejection as a bad input.

**12. pipeline.py:357** - A target that was scored and ranked but never carried into the pairing pool gets its own verdict NOT_ASSESSED, not UNSUITABLE.

> *Without it:* Collapsing it into UNSUITABLE told a caller reading only the status that the platform had assessed and rejected a target it never assessed at all. "No architecture was routed" means "it was not among the candidates considered", which is a different statement from "no architecture fits".



**13. pipeline.py:264** - The packages are built after the run dict exists and are folded into it, rather than being assembled from arguments.

> *Without it:* The package's input is the whole run, so passing its parts as arguments means listing eleven of them and keeping that list correct. Building the dict first and adding the packages to it makes the stage's dependency exactly what it is: everything the run produced.
## `car_pipeline/api/server.py`

**0. server.py:210** - The pair view's docstring names the set it actually returns, which is every measured pair, not the admissible ones.

> *Without it:* It said "the admissible pairs". The handler filters on `coverage.measured` - 19,110 pairs - while `admissible` is a different property that holds for 290. The claim was invisible at the default `limit=50`, because the rows are sorted by combined risk ascending and the first fifty happen to be admissible; raise the limit and the response fills with pairs the docstring says are not there. The response envelope was always honest, reporting `evaluated`, `measured` and `returned` by name, which is what makes the docstring the only wrong part. An adversarial panel refuted this flag on the strength of the default response, which is the observed-output trap in miniature.


**1. server.py:3** - A screen is submitted as a job and polled, rather than answered in the request.

> *Without it:* A run reads a 9 GB single-cell matrix, evaluates 19,900 pairs and makes a network call per pool member; a synchronous endpoint would time out in every client.

**2. server.py:8** - A run that assembles nothing answers HTTP 200 with a `status` and a `reasons` list — never 404, never 500, never a bare empty list.

> *Without it:* An empty list reads as "we looked and there was nothing to say"; the truth is that something specific and measured stopped each design, and the caller needs that rather than the absence of it.

**3. server.py:40** - Completed results are capped at MAX_RESULTS = 8 and evicted oldest-first, never the run in flight.

> *Without it:* Each result holds the ranked surface proteome and 19,900 evaluated pairs. The endpoint is open and the single instance has nowhere to shed load, so without a cap a caller creating projects in a loop takes the process out of memory — at which point the job table dies with it and every outstanding poll answers 404.

**4. server.py:145** - A finished result is popped and re-inserted into RESULTS rather than assigned in place.

> *Without it:* Assigning in place keeps the dict's original insertion position, so the oldest-first eviction would drop the most recently completed run instead of the least recently completed one.

**5. server.py:99** - Exactly one run at a time across the whole process, not one per project.

> *Without it:* Per-project was not enough and could not be fixed by deployment: a run is a detached thread and the request returns 202 in milliseconds, so a request-concurrency limit never sees two runs overlap. Two different projects each start a thread, both write the shared binder cache, and its writer replaces the payload before the manifest — a reader landing between those two steps gets an artifact nothing blessed.

**6. server.py:70** - An empty-string target_antigen is rejected with 400 instead of being treated as absent.

> *Without it:* Absent selects discovery mode; empty is a caller who meant to supply something. Coercing empty to absent silently runs a discovery screen for someone who asked a validation question.

**7. server.py:75** - Create_project resolves the indication at creation time so an unconfigured cancer type is refused immediately.

> *Without it:* Deferring the check to run time means the caller learns only after a multi-minute run — or receives another indication's results under their own name.

**8. server.py:188** - `usability` and `unavailable` are attached to every collection view rather than living behind their own endpoint.

> *Without it:* They were computed and returned by the pipeline and served by nothing, so an HTTP caller received a fully-formed ranking with no indication that a component was unmeasured — the "served with a caveat" outcome the pipeline refuses, except the caveat was not served either. A caller reading /targets should not have to know to ask elsewhere whether the ranking beneath it is supported.

**9. server.py:223** - Targets_view sorts before slicing, using the same key the pool uses.

> *Without it:* Stage3.rank() emits in universe order, not ranked order, so an endpoint called /targets would return unranked rows under a ranked name — and the first N would not be the top N. Using the pool's key keeps the two orders in agreement.

**10. server.py:265** - Pairs_view sorts on `(combined is None, combined or 0.0)` rather than the obvious `combined or 9` fallback for missing risk.

> *Without it:* `or 9` treats a combined risk of exactly 0.0 — the safest pair reachable, and a real measurement, since a per-organ score of 0.0 is measured — as the worst, sorting it off the page.

**11. server.py:281** - The span-confounding caveat is emitted as a field on every pair row, not as a footnote or a note in the docs.

> *Without it:* Per-cell detection tracks genomic span (rho +0.68) more strongly than expression (+0.20), so a reader who sees coverage_f_ab without the percentile beside it reads a span artefact as a coverage measurement.

**12. server.py:307** - Constructs are reported in three states — COMPLETE, AWAITING_BINDER, BUDGET_EXCEEDED — instead of counting only designs that fit and carry a sequence.

> *Without it:* These are genuinely different answers: a design that fits with its sequence is finished, one that fits but whose binder sequence was never supplied is waiting on one part, one that does not fit is stopped by the budget. Counting only the first, as this view used to, reported two designs where there were ten.

**13. server.py:322** - Over-budget designs are returned in the rows rather than filtered out.

> *Without it:* The overage is the measurement; dropping the row leaves a caller unable to see how far past the budget the design ran.

**14. server.py:331** - `binder_supplied` is False where a part declares a size but no sequence; the residues are not invented to fill the field.

> *Without it:* The length and the domain map are real measurements; fabricating residues to complete the record produces a construct that looks orderable and is not.

**15. server.py:352** - Every reason string is computed from this run's own numbers, with no hardcoded fallback figure in the prose.

> *Without it:* Prose carrying a remembered number presents it as a measurement of the run being served. The zero-design wording is likewise only emitted when there are no survivors, so it is never printed over a run that produced designs.

**16. server.py:394** - When every fitting design lacks a binder sequence the status is BUILDABLE_AWAITING_BINDER, not BUILDABLE.

> *Without it:* Reporting plain BUILDABLE promises something that cannot be ordered yet: the adaptor receptor binds a tag, no anti-tag binder exists in the connected sources, so its size is declared and its sequence is not invented.

**17. server.py:423** - Validation_view on a discovery run returns status NOT_APPLICABLE with an explanation, rather than an empty verdict.

> *Without it:* An empty verdict reads as "unsuitable" — the platform would appear to have rejected a target the project never supplied.

**18. server.py:464** - A run with no survivors says that nothing arrived to be ranked, rather than returning an empty ranking.

> *Without it:* An empty ranking reads as "nothing ranked highly"; the true statement is that nothing reached the ranking stage at all, and every drop upstream is a measurement rather than a failure of the stage that made it.

**19. server.py:170** - RunNotComplete (409) and NotFound (404) are separate exceptions with separate messages.

> *Without it:* "No completed run yet" is a statement about the job and is fixed by POSTing a run; "the run completed and does not contain this identifier" is not. Telling the second caller to re-run sends them round a loop that cannot help.

**20. server.py:601** - Both handlers end in a bare `except Exception` that still writes an HTTP response.

> *Without it:* Letting an exception escape the handler closes the connection, and a closed connection is the one failure a client cannot interpret.

**21. server.py:642** - The job dict is copied to a snapshot inside _LOCK before being serialised.

> *Without it:* The worker adds a "trace" key on failure, so serialising the live dict can raise "dictionary changed size during iteration" — a poll that fails precisely when the job it describes has failed.

**22. server.py:610** - `?limit=` is parsed from the query string and clamped to max(1, min(500, n)), with a malformed value falling back to the default.

> *Without it:* The parameter was accepted in the view signature and never passed, so it silently did nothing; unbounded it lets one request serialise the whole ranked proteome.

**23. server.py:700** - The default bind address is 127.0.0.1; 0.0.0.0 must be requested explicitly.

> *Without it:* A 0.0.0.0 default means a laptop demo is also serving the internet.

**24. server.py:705** - PORT is read from the environment after argparse runs, with --port defaulting to None, rather than being converted while building the parser.

> *Without it:* Converting the environment value in the parser default means an empty or malformed PORT raises an uncaught ValueError before argparse ever runs, so an explicit --port could not override it — the flag that exists to fix the problem becomes unreachable.

**25. server.py:241** - No count of criteria appears in a reason the service cannot compute.

> *Without it:* The pairing reason read "Stage 4 closed at 10 of 14 criteria with four documented limitations." The service does not run the suite and cannot know either number, so both were transcribed once and left there. A criterion was later added: the stage now closes at 10 of 15 with five limitations, and the sentence had been wrong on both figures ever since, on a live endpoint. Naming the limitations instead of counting them says more and cannot go stale.

**26. server.py:517** - Risk attribution hangs off the existing evidence trail instead of becoming its own endpoint.

> *Without it:* Two endpoints would report the same risk from different code, and the number they disagreed about would be the one under dispute. The evidence trail already answers per gene and already carries `risk` and `risk_organ`; the attribution belongs beside them, where a reader sees the verdict and its grounds in one response.
**27. server.py:496** - The package is served whole at one endpoint and per candidate at another, and both carry the gaps section.

> *Without it:* A per-candidate endpoint that omitted the gaps would serve eight sections that read as twelve, which is the failure the whole stage exists to avoid.


---

# Running it


## `run_all.py`

**0. run_all.py:29** - The HTTP surface runs once per indication, and a stage carries its own arguments.

> *Without it:* Every no-survivor path in the service is dead code under the pancreatic indication, which always has surviving designs. The refusal reasons, the attrition explanation and the conservative-backup counts are reached only by an indication that returns nothing, and one of them held a NameError that answered HTTP 500 for months. Testing the surface against the indication that never enters those states is the gap that let it live. Two stages then share one script, so the log filename takes the stage number to break the tie rather than the second silently overwriting the first.

**1. run_all.py:4** - Each stage's verifier runs as its own subprocess rather than being imported and called in-process.

> *Without it:* A stage that dies takes only its own interpreter with it and the run carries on, so a late failure still reports every earlier stage instead of losing the lot.

**2. run_all.py:8** - `--fresh` deletes only the derived artifacts (data/stage4, data/stage5); the raw source caches under data/ are never touched.

> *Without it:* The 11 GB of matrices and atlases are the input to the run, not part of what it recomputes; deleting them turns a re-run into a multi-hour re-download.

**3. run_all.py:33** - Stage order is dependency order, and Stage 5 must run before 6 and 9.

> *Without it:* Stage 5 is the only verifier that exercises the binder retrieval route live; 6 and 9 read the cache it blessed, which is the path the API takes. Run standalone against a warm cache, those two cannot see a dead retrieval route — which is why Stage 5 retrieves for real even though it costs five minutes.

**4. run_all.py:49** - The multi-indication verifier runs last.

> *Without it:* It runs BOTH indications and then asserts that neither one's artifacts moved, so it must see the caches in the state every earlier stage left them; run earlier, it would be asserting against caches later stages then rewrite.

**5. run_all.py:58** - Known-tripping criteria are listed in ACCEPTED with the decision that accepted each, instead of the runner failing on any trip.

> *Without it:* Five criteria trip on every run, so a runner that failed on any trip would be red forever and a genuine regression would land in a report nobody reads as changed. A criterion tripping that is not on the list is a regression.

**6. run_all.py:274** - The accepted list is checked in both directions: an entry that STOPS tripping is reported too, including identifiers the run no longer emits at all, while a stage that reported nothing is excluded from the check.

> *Without it:* A stale entry silently grants an exemption nothing needs, which is how an accepted-limitations list turns into a place to hide things. R13 is scheduled to disappear — it was withdrawn and replaced — so checking only identifiers the run still emits would let its entry sit there forever, pre-approving whatever later reuses the name. A silent stage is excluded because it cannot distinguish "gone" from "never ran".

**7. run_all.py:69** - Stage 3's R13 was withdrawn rather than fixed by reweighting.

> *Without it:* The two populations differ by construction (breadth 51 vs 7) and the max-over-sources gate rewards being unmeasured, so no combination rule reaches the 5x limit. It was replaced by R13-prime, which clears.

**8. run_all.py:74** - Coverage (f_AB) is reported beside a span-matched percentile and removed from partner selection rather than used to gate pairs.

> *Without it:* F_AB tracks genomic span (+0.68) more than expression (+0.20), so selecting partners on it selects long genes, and the resulting coverage number looks like a measurement of co-expression.

**9. run_all.py:76** - A cleared pair containing HLA-A is recorded rather than filtered out.

> *Without it:* The pool is not curated by hand; removing the embarrassing hit by name would hide that the scoring admits it, and the same hole would stay open for every unnamed equivalent.

**10. run_all.py:77** - An unmeasured antigen is kept as a third state rather than imputed.

> *Without it:* 48.6% of cleared pairs stop clearing if an unmeasured antigen saturates its organ. That sensitivity is the honest cost of not imputing; imputing a value would replace it with a clearance rate that reads as measured.

**11. run_all.py:83** - The terminable ceiling stays at the spec's 0.35 even though a pre-registered pin (A6) expected MSLN to route to an adaptor.

> *Without it:* MSLN matches that row's condition in words — serious normal-tissue expression — but its measured risk 0.6366 is nearly twice the declared ceiling of 0.35. Admitting it needs a ceiling near 0.65, which also admits about 120 others: a clinical policy decision, not a code change. A9 reports the whole sweep so the trade stays visible.

**12. run_all.py:94** - The criterion regex accepts both indentation spacings used across the verifiers.

> *Without it:* Matching only one spacing silently reports zero criteria for half the stages — stages that ran and passed would appear to have verified nothing.

**13. run_all.py:99** - A separate `_CHECK` pattern parses the spec verifier's named checks alongside the numbered-criteria pattern.

> *Without it:* Without it the spec verifier's 31 checks showed a count in the summary and nothing underneath it, which is the one shape a verification report must not have.

**14. run_all.py:103** - The surface verifier's two summary patterns are compiled with re.M.

> *Without it:* Without it the `^` anchors only match at the start of the whole captured output, so the surface stage parsed as zero criteria while exiting 0 — a stage that passed and reported nothing.

**15. run_all.py:128** - Counts are read out of each verifier's own report rather than recomputed from the parsed criteria.

> *Without it:* The verifier is the authority on its own criteria; a second count living in the runner is free to disagree with it, and then neither number can be trusted.

**16. run_all.py:122** - `Stage.ok` requires total > 0, and a stage that reported no criteria at all is called out separately from a failing one.

> *Without it:* An empty verification is not a passing one: exit 0 with nothing parsed means the verifier never checked anything, and without this it would be counted as clear.

**17. run_all.py:178** - A stage timeout is recorded as exit 124 with the partial output kept, rather than being allowed to hang the run.

> *Without it:* A stalled socket in one verifier would otherwise block every remaining stage; the timeout is generous against the slowest stage (binder retrieval, ~320 s) but finite at 1800 s, so one stalled stage fails alone.

**18. run_all.py:187** - Each stage's transcript is written to disk as soon as that stage finishes, not buffered until the end of the run.

> *Without it:* Buffering every transcript until all stages finish means an interrupt during stage 5's five minutes discards the nine that already succeeded.

**19. run_all.py:193** - The API section of the report quotes the API verifier's own closing block instead of recomputing it.

> *Without it:* A recomputed summary is free to disagree with what the API actually returned, and the report is supposed to be evidence of that output.

**20. run_all.py:301** - Preflight() checks the launching interpreter for the declared dependencies and stops before any stage runs.

> *Without it:* Every stage runs under sys.executable, so launched with the wrong Python all ten fail identically on the first import and the report says "no criteria parsed" ten times instead of naming the one cause.

**21. run_all.py:308** - The dependency list is parsed from requirements.txt, skipping option lines, rather than written out in run_all.py.

> *Without it:* A list written here is free to drift from what the project actually declares, so the preflight would pass while a real dependency was missing; an option line (`-r`, `--index-url`) is not a module and must not be probed as one.

**22. run_all.py:326** - Preflight emits two different remedies depending on whether the current interpreter already IS the venv interpreter.

> *Without it:* Telling someone already running the venv interpreter to run the venv interpreter is not advice; the real fix in that case is pip install -r requirements.txt.

**23. run_all.py:402** - A structured JSON run record is written beside the Markdown report for downstream generators.

> *Without it:* A regex pinned to render()'s f-strings would break on a formatting tweak and, worse, could half-match and render empty sections — a published page that looks complete and is missing results.

**24. run_all.py:426** - "Crashed" is defined as a non-zero exit with criteria parsed and nothing tripped, kept distinct from a non-zero exit caused by a tripped criterion.

> *Without it:* A verifier exits non-zero BECAUSE a criterion tripped, so a non-zero exit alongside tripped criteria is the expected shape, not a fault. Conflating them produces a message that is simply untrue on every normal run.


## `bootstrap.py`

**1. bootstrap.py:14** - The packaged release carries the derived `group_means.npz` instead of asking every clone to rebuild it.

> *Without it:* Deriving it streams the whole 8.3 GB single-cell matrix and took 2 h 19 min on the machine that first built this cache, against about 15 minutes for every download combined; the archive skips that by carrying the 5.7 MB result.

**2. bootstrap.py:20** - Bootstrap is standard-library only except inside `--from-sources`, which imports the pipeline lazily.

> *Without it:* A fresh clone has no dependencies installed, and bootstrap is exactly the tool you reach for then. Importing the pipeline at module level makes it unrunnable in the situation it exists for.

**3. bootstrap.py:50** - HEAVY is the single authority on archive exclusions; `_members()` iterates it rather than repeating the patterns.

> *Without it:* Two copies of the exclusion list drift, and then what is printed as skipped is not what is actually skipped.

**4. bootstrap.py:56** - BUILD_ONLY marks the single-cell archive and matrix as payloads a deployed cache is expected NOT to have, and their absence is reported as intended rather than as a gap.

> *Without it:* A served run never opens either — measured, not assumed. Reporting them as missing marks every correctly deployed cache incomplete, and the 8.3 GB matrix and 2.6 GB archive would be re-downloaded for nothing.

**5. bootstrap.py:61** - The trials cache is listed as "deferred" rather than missing, because it is fingerprinted by the antigen list it covers.

> *Without it:* It has no meaning until a pool exists and is built during the first screen; listing it as merely missing would report a complete cache as incomplete forever.

**6. bootstrap.py:85** - Tumour-side caches are enumerated per registered indication instead of as flat entries named "tcga" and "singlecell".

> *Without it:* The old flat list treated each as a single thing, so a clone provisioned the reference indication and reported itself complete — "10/10 sources usable" — while a second registered indication had nothing.

**7. bootstrap.py:100** - One `_sha256` implementation is shared by the packager and the verifier.

> *Without it:* Two copies are free to diverge, and a divergence there produces a checksum mismatch on a perfectly good archive with no way to tell which side is wrong.

**8. bootstrap.py:113** - `_state` opens each manifest and checks the payload file it names, instead of globbing for manifests.

> *Without it:* A manifest only marks that a fetch finished. Globbing for manifests reports a cache as complete while the payload beside it has been deleted — the failure this tool is most likely to be asked about.

**9. bootstrap.py:147** - `_indications()` swallows any import failure and returns an empty list.

> *Without it:* Bootstrap runs before dependencies are installed on a fresh clone, so it must degrade to the shared-source report rather than dying on an import — otherwise the first command a new machine runs is the one that cannot run.

**10. bootstrap.py:277** - An unreadable or empty checksum sidecar is treated as a mismatch (exit 2), not as a crash.

> *Without it:* A truncated checksum file is exactly the state this machinery exists to survive, so it must not be the thing that raises.

**11. bootstrap.py:294** - A checksum mismatch refuses to unpack rather than unpacking and warning.

> *Without it:* A truncated transfer produces a cache that reads as present and answers with the wrong data, which is worse than having no cache at all.

**12. bootstrap.py:307** - `tar.extractall(..., filter="data")` is used with no compatibility fallback for older interpreters.

> *Without it:* The filter is what refuses absolute paths and parent traversal in member names. An unpacker that silently drops that guard on an older Python is worse than one that refuses to run there; this project runs on 3.13.

**13. bootstrap.py:344** - `from_release` passes require_checksum=True, so a downloaded archive with no sidecar is refused.

> *Without it:* Unpacking 298 MB straight off the network unverified is the exact failure mode the checksum exists for; the local `--from-archive` path may proceed unverified, a download may not.

**14. bootstrap.py:320** - The release download goes through the `gh` CLI and reuses the clone credential rather than inventing a token.

> *Without it:* The repository is private and the asset inherits that, so an unauthenticated fetch returns 404 — the right answer but an unhelpful one, indistinguishable from a missing release.

**15. bootstrap.py:386** - The rebuild calls `fetch()` on each source, and names TCGA's cohort build and the single-cell group means explicitly as separate steps.

> *Without it:* `fetch()` is the method that populates a cache and several sources have no `load()` at all (GTEx exposes only match_surface, which needs the surface proteome). For the two derived artifacts the download alone leaves the expensive part undone, so a cache that looks fetched is not usable.

**16. bootstrap.py:392** - Rebuild steps are ordered cheapest-first.

> *Without it:* A broken network then fails in seconds rather than after the 2.6 GB download.

**17. bootstrap.py:403** - Per-indication rebuild steps are generated in a loop over the registry, each lambda binding its indication with a default argument.

> *Without it:* This used to be three unparameterised calls, so a clone provisioned the reference indication only and the multi-indication stage failed on a fresh machine while every other stage passed. Without the `i=ind` default binding, every queued lambda would close over the last indication instead of its own.

**18. bootstrap.py:428** - An interrupt returns 130 and says nothing is corrupted, because each artifact is committed with its manifest only once complete.

> *Without it:* Without the commit-on-completion rule a partial file is mistaken for a finished one and a re-run resumes on top of it; with it, a re-run resumes rather than restarts the 2 h 19 min step.

**19. bootstrap.py:434** - A failed rebuild step prints the full traceback, not just the exception message.

> *Without it:* A failure two hours into the last step with one line of context means re-running two hours to find out where it was.


## `make_artifact.py`

**1. make_artifact.py:22** - Stages 7 and 8 are rendered as a named gap rather than omitted, because they are absent from the pipeline, not from the page.

> *Without it:* A reader counting 1 to 11 and finding nine entries cannot tell a stage that was never built from two dropped out of the report. Showing the gap is the information.

**2. make_artifact.py:110** - A tripped criterion with a recorded decision behind it renders differently from one without.

> *Without it:* Rendering them identically is exactly the silencing this page claims not to do: an accepted limitation and an unexplained regression would look the same to whoever reads the page instead of the run.

**3. make_artifact.py:147** - The number that got through is read from the attrition chain's own last remaining value, never restated as a constant.

> *Without it:* This number has already changed once. A page asserting the old one while the run reports the new one is precisely the failure the page exists to surface.

**4. make_artifact.py:162** - The ceiling sweep is lifted out of its criterion line and given its own section.

> *Without it:* It is the one number this pipeline cannot measure. What each setting buys belongs in front of whoever is entitled to choose it, not buried inside a criterion line that reads as passed.

**5. make_artifact.py:201** - The summary line is written from the accepted list rather than asserting that every trip is a recorded limitation.

> *Without it:* One of them is not. The accepted list exists to make that sentence true; hardcoding the claim makes the page lie exactly when a genuine regression appears.

**6. make_artifact.py:499** - The report path is an argument, not a constant, because the runner accepts --report.

> *Without it:* A fixed path renders a previous run's page while printing a confident success line about the current one.


## `car_pipeline/stages/validation.py`

**1. validation.py:219** - The two design classes are defined from the architecture table, not from what survived: conservative is the conventional single-antigen receptor carrying a binder with clinical precedent, advanced is an architecture the spec lists as non-conventional, meaning a gated dual or an adaptor.

> *Without it:* Both labels are the kind that get filled in by whatever is available. Defining conservative as "the safest thing in the pool" labels a dual conservative on a run where no single survives, and defining advanced as "the rest" labels a plain receptor advanced. The definitions have to be fixed against the architecture table so the honest answer stays reachable, which is the answer this pool actually gives: no conservative design exists here, reported rather than filled.

**2. validation.py:237** - Every quantity in the no-conservative-backup explanation is counted from the constructs, including the ones that read as fixed background.

> *Without it:* The sentence carried two hardcoded claims inside prose that otherwise looks computed. It said "the only single-antigen target, NPSR1" where three were recommended - MSLNL, NPSR1 and ZPLD1 - and it asserted "the two dual designs that assemble are both over the payload budget" where zero duals assemble and no construct is over budget at all. Both were true when written and neither was ever recomputed. This is the standing correction from the stale-pin finding, in output prose rather than in a criterion: state the property, derive the subjects from the run. It is worse here than in a criterion, because a criterion that goes stale eventually trips, and prose never does - it is the platform's own front-page explanation of its result, and it was confidently wrong to every reader.

## `make_brief.py`

**1. make_brief.py:26** - The targets the brief works through are named here, in the report generator, and nowhere in the pipeline.

> *Without it:* The attribution facility has to stay target-agnostic - that is the platform's first rule and T8 enforces it on the code path. But a brief with no worked example explains nothing, and the example has to be some particular target. Putting the selection in the generator keeps both true: the pipeline knows no gene, and the report names two.

## `make_package.py`

**1. make_package.py:265** - One document per candidate, and the gaps section is repeated in full in each rather than referenced.

> *Without it:* A reader receives one candidate's package, not the set. A pointer to a shared appendix is a pointer they may not follow, and the whole purpose of the section is that what the package cannot tell them travels with what it can.

## `strip_comments.py`

**1. strip_comments.py:41** - The strip works from a token stream, not a regular expression over lines.

> *Without it:* A pattern matching the comment character corrupts every string that contains one, and this repository has several: URL fragments, format strings, and the attribute-doc markers. Only the tokeniser knows which occurrence starts a comment and which sits inside a literal.

**2. strip_comments.py:52** - An f-string's middle token is re-escaped before being written back, doubling every brace.

> *Without it:* The tokeniser hands back the *unescaped* text, so a source doubled brace arrives single. Written out verbatim it stops being a literal brace and becomes a format field - which silently rewrote a regex fragment matching a twelve-character identifier into one matching a repetition count. Nothing about the file would have looked wrong.

**3. strip_comments.py:58** - The column for that token is advanced by the re-escaped text rather than by the reported end column.

> *Without it:* The tokeniser ends the token at the first brace of an escape pair, so the reported column understates what was consumed and every following token on the line looks displaced by one - inserting a space inside the literal.

**4. strip_comments.py:64** - A line-break token's own text is dropped, because the row-gap arithmetic already emits one.

> *Without it:* Appending both double-spaces the entire file. The tree comparison cannot catch it: blank lines never reach the parser, so the guard passes while every file silently doubles in length.

**5. strip_comments.py:84** - Docstrings are shortened through the parse tree, not by matching triple quotes.

> *Without it:* A docstring can contain triple quotes, and a triple-quoted string that is not the first statement is not a docstring. Only the tree tells the two apart.

**6. strip_comments.py:120** - Two blank lines survive before a top-level definition and one anywhere else, and none directly under a colon or inside an open bracket.

> *Without it:* Deleting a comment paragraph otherwise leaves a hole the exact size of the comment - a blank line as the first line of a block, or a gap in the middle of a literal, in places no one would write one.

**7. strip_comments.py:137** - No file is written whose parsed tree changes, with docstrings normalised out of both sides first.

> *Without it:* This is the only failure the tool can have, and it is silent: an edit that removes a statement rather than a comment produces a file that still parses and still imports. The guard has already refused a write twice - once for the brace escaping above, and once more before that.

**8. strip_comments.py:65** - A docstring is reduced to its opening paragraph *joined* onto one line, not to its first physical line.

> *Without it:* Where an author wrapped a summary sentence across two lines, taking the first physical line deleted the tail and left a fragment ending on a dangling article. Two shipped that way and were found only by reading them: "so the" and "the smallest that fits is the". The guard at decision 7 cannot see this class at all - it blanks every docstring before comparing trees, which is exactly what makes the comparison safe for code and blind to prose. Joining loses nothing: paragraphs after the first are still dropped, which is the stated purpose, but no sentence is ever cut mid-clause.

**9. strip_comments.py:186** - Every summary the tool joins is listed by file and line at the end of the run.

> *Without it:* The tree guard is structurally incapable of reporting docstring changes, so a rewrite that shortens prose leaves no trace anywhere. Printing what was joined is the only signal that the tool touched a summary, and it is what makes decision 8's behaviour reviewable rather than merely better.

**10. strip_comments.py:16** - The skipped trees are matched as resolved paths under the project root, not by directory name anywhere in the path.

> *Without it:* The skip list names `data` for the cache at `data/`. Matching the name against every path component also matched the source package `car_pipeline/data/`, so thirteen modules and three hundred comments were never visited, and the check that confirmed the sweep reused the same predicate and confirmed its own blind spot. The two entries look interchangeable and are not: `reports` is safe under either rule and `data` is not, so the next person shortening this back to a name test sees no failure until a package happens to share a label with a cache. The general form of the error is written up in `specs/verification-sharing-assumptions.md`, instance 7.


---

# The criteria

## `verify_schema.py`

**1. verify_schema.py:28** - Unknown fields are rejected, not absorbed.

> *Without it:* Absorbing them means a misspelled field silently does nothing and the caller is told the request succeeded.

**2. verify_schema.py:38** - A blank antigen from a form is treated as an absent one, and a supplied antigen flips the mode.

> *Without it:* An empty string is what a browser sends for an untouched field. Treated as a value it selects validation mode against a target nobody named.

**3. verify_schema.py:120** - The returned spec must not alias the shared indication config.

> *Without it:* An aliased spec lets one run mutate the configuration every later run reads.

**4. verify_schema.py:124** - Back-to-back builds must not collide.

> *Without it:* Two specs built in the same tick sharing an identifier would silently overwrite each other's results.


## `verify_surface.py`

**1. verify_surface.py:12** - The headline figures are measured against the pinned proteome release rather than reconstructed, with a tolerance kept only because the pin is enforced on the response rather than the request.

> *Without it:* A release bump is caught by the fetch, but nothing here should fail on the last digit of a count. The exact validation sets below remain the real test.

**2. verify_surface.py:32** - The rejects are asserted by name, not by count.

> *Without it:* These entries were admitted for a while because a plasma-membrane phrase appeared in a free-text note rather than a location statement. A regression here is silent: the count moves by fourteen and nothing else looks wrong.

**3. verify_surface.py:42** - Two named entries stand for the note-only class: one describes transient trafficking, the other lipid binding.

> *Without it:* Neither is a statement about where the protein rests, and a substring test over notes reads both as evidence that it is on the surface.

**4. verify_surface.py:71** - The subdivision of the withheld set is reported, not gated.

> *Without it:* The boundary rests on a compartment vocabulary that had to be inferred rather than read. No protein's fate depends on which side it lands, so gating on it would turn an inference into a rejection.

**5. verify_surface.py:140** - The topology criterion explicitly excludes the note-only state.

> *Without it:* A bare not-outward test is now satisfied by the note-only state as well, so without the exclusion the criterion stops isolating the gate it names.

**6. verify_surface.py:165** - Every entry in the third state is listed by name.

> *Without it:* Fourteen is small enough to read, and a silent change in the set is exactly what a count alone hides.


## `verify_ranking.py`

**1. verify_ranking.py:104** - The two axes are put on one scale by measurement before anything is scored, and the curve is reported and hashed as part of the run.

> *Without it:* Asserting the correspondence between a staining axis and a transcript axis bakes in a conversion nobody measured, and hides it from every reader of the result.

**2. verify_ranking.py:123** - Every accession carrying a symbol is tested, not just one.

> *Without it:* Several symbols carry more than one accession. Keeping only the last lets a second entry under the same name go untested, which is exactly how a ubiquitous protein clears the ceiling unnoticed.

**3. verify_ranking.py:173** - Hash stability across processes is checked here rather than assumed.

> *Without it:* A hash seeded by anything process-local defeats itself silently: every run looks reproducible within itself and nothing matches between runs.

**4. verify_ranking.py:303** - A run that scored nothing at all must trip the criteria, not raise.

> *Without it:* The degenerate run is the state these criteria exist to catch. Raising turns the detection into a crash and loses every other criterion's result with it.

**5. verify_ranking.py:338** - The offset criterion is tested against the measured systematic offset, not against parity.

> *Without it:* A parity test flagged all 25 and told a reader nothing.

**6. verify_ranking.py:350** - A fold built on a floored denominator reads differently from one that was measured.

> *Without it:* It is a lower bound, not a measurement, and presenting the two identically invites a conclusion the evidence does not support.

**7. verify_ranking.py:409** - Component drift is reported on every run, not only when it trips.

> *Without it:* A composite drifting toward one component stops being a multi-criteria score long before it reaches the rejection threshold. Only the trend shows that.

**8. verify_ranking.py:355** - T5 counts the winning organs each arm decided, separating those carrying a non-zero score from the zero-scored ties, and requires both arms to be present among the former.

> *Without it:* T1 to T4 are identities: an attribution correct about every target it describes satisfies them, including one that describes none and one that only ever sees a single arm. T5 is the criterion that would catch a facility exercised on half its inputs. The separation matters because a protein absent from every organ scores zero everywhere, so every one of its organs ties for the maximum; counted raw that is 4,893 baseline winners, of which 2,769 carry a score above zero. Reporting only the raw figure would claim more evidence than the run contains.

**9. verify_ranking.py:414** - T8 reads the attribution functions' string constants out of the parse tree and intersects them with the run's own gene symbols.

> *Without it:* "No gene is seeded" is the platform's first rule, and an attribution facility written around one target would break it quietly - the code would look general and behave specially. Checking literals structurally rather than by grep means the criterion cannot be satisfied by spelling a symbol differently, and taking the symbol set from the run rather than a list means it covers every gene the screen actually saw.

**10. verify_ranking.py:339** - T4 recomputes each organ's staining presence from the atlas entry rather than reading it off the attribution it is checking.

> *Without it:* As first written T4 asked whether a staining score appeared where the protein arm was unmeasured, and whether a winning arm reported itself absent. Both are impossible by construction: `attribute_risk` only records a staining reading when it found one, and only names an arm that supplied a value. The criterion could not fail. Deriving the stained organs independently from `atlas_gene.staining` gives it something to disagree with, and it does: blinding the attribution to the protein arm trips it on 34,623 organ rows.


## `verify_ranking_final.py`

**1. verify_ranking_final.py:85** - The front check uses a synthetic case where one candidate is dominated on every objective.

> *Without it:* Real data may contain no strictly dominated pair, so the criterion would pass without ever exercising the dominance test.

**2. verify_ranking_final.py:27** - Binder records are read through the same hash-gated path its four sibling verifiers use, rather than straight out of the cache slot.

> *Without it:* `read_binders` validates the payload against its own digest and nothing else - not the Stage 4 hash, not the gene set. The slot is shared, and its last writer is the deliberately degraded run inside the multi-indication verifier, whose pool differs by two genes. Read against that state this verifier reported 6 of 6 clear with nothing able to notice, and only the order of the suite kept it from mattering on a run that reuses derived artifacts.

**3. verify_ranking_final.py:59** - The gate is given the constructs, for the reason recorded against `verify_safety.py`.

> *Without it:* The same defect and the same silence: an adaptor receptor's binder is invisible to the origin check unless the constructs carrying it are passed in.


## `verify_pairing.py`

**0. verify_pairing.py:163** - The declared tolerances are passed to `decide`, so the artifact this verifier persists is the decision set the platform ships.

> *Without it:* Passing none is a positive instruction to disable routing, so no adaptor row can exist in the artifact. Five verifiers read that file as the Stage 4 result, and the adaptor is the architecture every surviving design in the worked indication uses. The construct stage assembled nothing, the safety gate never reached its terminable branch, and the final ranking reported that no design reaches the end while the service returned five. The file is named for the stage's result and has to hold it.

**1. verify_pairing.py:36** - A known-answer check sits on the per-cell derivation, keyed to a gene whose answer is known independently.

> *Without it:* Column indices are stored in descending order within a row, so a lookup assuming otherwise returns zero for every gene without ever erroring. Only a known answer catches a silent zero.

**2. verify_pairing.py:41** - The highest-ranked target is deliberately excluded from that check.

> *Without it:* It carries two molecules across all 64,538 malignant cells - a capture failure the ranking stage already documents, not a broken derivation. Including it would stop a correct implementation on the pool's top target.

**3. verify_pairing.py:47** - The control gene is requested as an extra column rather than drawn from the pool, and is one the surface filter can never admit.

> *Without it:* That is what makes it a control: it is the loudest malignant signal in the atlas and it is independent of anything the ranking decided. A control drawn from the pool is entangled with the thing being checked.

**4. verify_pairing.py:155** - The span annotation is allowed to fail without losing the run, and is not registered as a connector.

> *Without it:* It gates nothing and the preflight therefore cannot have checked it. A hard failure this late would discard everything computed above for a source no decision depends on.

**5. verify_pairing.py:177** - A non-finite median is stored as not measured, never as a low value.

> *Without it:* Stored as a number it would make the gene ineligible while still counting as having a measurement, which reads as measured-and-found-absent rather than never-measured.

**6. verify_pairing.py:204** - The agreement tolerance absorbs only the rounding boundary, and the worst excess is printed alongside.

> *Without it:* Both numbers are carried at four decimals. A tolerance wide enough to be safe would let a real disagreement hide underneath it; printing the worst excess means it cannot.

**7. verify_pairing.py:275** - The artifact records which criteria were checked, not only which failed.

> *Without it:* A criterion silently dropped from the run is then visible in the manifest instead of reading as a criterion that passed.

**8. verify_pairing.py:313** - Clearance must be decided on the conservative arm; a pair clearing only on the optimistic arm is unresolved by definition and must not be marked cleared.

> *Without it:* Deciding on the optimistic arm converts an assumption about unmeasured organs into a safety claim.

**9. verify_pairing.py:335** - The substantive check counts cleared pairs that would stop clearing if the unmeasured antigen saturated the organ nobody looked at.

> *Without it:* Having an unresolved organ is not itself the failure - the conservative arm already charges the measured member's score there. The failure is clearance surviving only because the missing antigen was assumed no more prevalent than the measured one.

**10. verify_pairing.py:351** - Two criteria were withdrawn and kept as reported numbers rather than deleted or left gating.

> *Without it:* They asserted properties enforceable only while the patient floor took part in selection. That floor inherits the span confound wholesale, and selection no longer uses it. Gating on it here would reinstate through the criteria the thing that was removed from the stage - so the number stays visible and stops deciding.

**11. verify_pairing.py:434** - The artifact is written whether or not a criterion tripped, and carries which ones did.

> *Without it:* A stage that persists only on a clean run leaves nothing on disk in exactly the state the project is actually in, forcing the next stage to re-derive the numbers under question. The manifest says plainly whether the payload may be read as a result, so it cannot be mistaken for one.

**12. verify_pairing.py:449** - Measurements are printed whether or not a criterion tripped.

> *Without it:* They describe what the atlas contains rather than what should be built. Withholding them behind a passing run hides the strongest thing this stage has to say about the architecture.

**13. verify_pairing.py:475** - The ordering question is measured rather than asserted: how often would the recommendation change if partners were ordered by the composite instead of the risk margin?

> *Without it:* Asserting that the ordering key is right proves nothing. Measuring the disagreement shows whether the choice is load-bearing.

**14. verify_pairing.py:560** - The illustrative pair is the most favourable one on the watch list, not the worst.

> *Without it:* An argument that holds on the best case does not depend on which pair was picked; one built on the worst case invites the reply that the pair was chosen to make it.


## `verify_routing.py`

**1. verify_routing.py:21** - The ceiling in the config is compared against the value recorded in the spec.

> *Without it:* Moving the ceiling without moving the spec is then a tripped criterion rather than a silent retune - and this is the parameter the pipeline is least able to justify from its own data.

**2. verify_routing.py:25** - The admitted count is reported across the whole ceiling range.

> *Without it:* The choice of ceiling becomes visible rather than argued, and a reader can see what the setting buys instead of taking the single number on trust.


## `verify_binders.py`

**1. verify_binders.py:16** - The known targets are imported from the ranking stage's own list rather than restated here.

> *Without it:* Two copies of the same list drift, and the check silently stops covering what the stage actually produces.

**2. verify_binders.py:20** - Every pin was verified present by accession-anchored search before the run, and every negative verified absent.

> *Without it:* Without the negatives the check is one-sided: a stage that returned everything for every target passes it completely.

**3. verify_binders.py:31** - The structure route carries must-appear pins, not only must-not-appear ones.

> *Without it:* This half is not optional. The first run of this stage returned zero structure candidates for all 200 targets because the two sources identify entries differently - and every structure-route check in place at the time was a negative, so all of them passed against a dead route. This is the origin of the positive-pin rule the rest of the criteria follow.

**4. verify_binders.py:42** - The negative pins separate two different absences: entries that exist but contain no antibody, and no entries at all.

> *Without it:* A stage emitting a candidate for either is not filtering, and collapsing the two hides which half of the filter failed.

**5. verify_binders.py:47** - One pin is set where the correct answer is a single specific count.

> *Without it:* It is the sharpest check in the file: six accession-anchored entries and exactly one containing an antibody. A stage echoing its entry count returns six; one that lost the join returns none. Only the correct implementation sits at one, so the check cannot be passed by either failure mode.

**6. verify_binders.py:137** - The filters are required to subtract - some targets partial, some echoing, and the correct answer is a mix.

> *Without it:* A stage returning a candidate per entry shows zero partials; one that lost the join shows every target partial and none echoing. Requiring the mix rules out both at once, where either half alone admits one of them.

**7. verify_binders.py:154** - The final tally is recomputed rather than read from the failures list.

> *Without it:* The early return above guarantees that list is empty by the time this line runs, so the criterion could never fail. A criterion that cannot fail is the thing this project keeps deleting.

**8. verify_binders.py:245** - The therapeutic line is built as one string, not a conditional expression over a concatenation.

> *Without it:* Written the latter way the condition binds to the whole expression, so a therapeutic with only one variable region printed a blank line and lost its name, stage and status with it.


## `verify_construct.py`

**1. verify_construct.py:26** - K2's pin set is derived from the run through `two_armed_duals`, and where it is empty K2 trips and says the two-arm join is untested.

> *Without it:* The pins were originally two named genes, chosen because they were then the only duals carrying a binder on both arms. The property moved off them: both now pair to partners with no binder. A criterion holding a result rather than a property passes on a regression, which is the Stage 4a end-state pin with the sign reversed. Deriving the subjects from the run keeps the criterion pointed at the property; tripping on an empty set is what stops it going green the day one dual assembles while the join is still uncovered.

**2. verify_construct.py:32** - The check reads the cache the binder stage wrote rather than retrieving again, and gates reuse on the recorded hash.

> *Without it:* That cache is the path the service takes, so the join under test is the real one. Passing the hash also persists the records when the cache is cold, so the next verifier does not spend the retrieval again.

**3. verify_construct.py:67** - The pins make a vacuous pass impossible.

> *Without it:* A stage that assembled nothing passes any check phrased only over what it assembled.

**4. verify_construct.py:124** - Both terms are recomputed from the segments rather than read from the construct's own properties.

> *Without it:* Comparing a stored headroom against its own definition is a tautology and cannot fail.

**5. verify_construct.py:141** - A binder is necessary and not sufficient: a construct is owed only where the target also carries a recommendation, and for a dual, where the partner has a binder too.

> *Without it:* Owing a construct on the binder alone makes the criterion trip on designs the pipeline was never going to build, burying real failures in noise.

**6. verify_construct.py:208** - Counts are checked against the manifest's recorded pool size, not this run's own input.

> *Without it:* Comparing the output to the input it was built from cannot fail.

**7. verify_construct.py:86** - K0 asserts the decision set is routed and yields at least one construct, and is reported before every criterion that reads the assembled set.

> *Without it:* The artifact was written with routing disabled, so nothing assembled, and K1, K3, K4, K5 and K6 each cleared over an empty list - K4 reporting "0 parts in the first construct". A criterion phrased over a population is satisfied by the population being empty. K0 states the precondition once where it is legible; each of the five also trips on an empty population in its own right, so removing K0 would not restore the silence.

**8. verify_construct.py:98** - The binder check is per route: a Stage 5 candidate for a single or a dual, the retrieved anti-tag part for an adaptor.

> *Without it:* Both K2 and K7 assumed every construct's binder is a Stage 5 sequence candidate. An adaptor's is not - it is retrieved from a deposited structure - so five correctly assembled constructs produced five failures under each, reading "chosen binder not in Stage 5" and "construct without a usable binder". Checking a construct against a route it does not use reports a defect in the criterion as a defect in the construct.

**9. verify_construct.py:66** - The anti-tag part is retrieved once and the adaptor clause compares the construct against that object's sequence, name and accession rather than against a literal.

> *Without it:* A literal copied into the criterion agrees with itself if the retrieval changes underneath it, which is how a dead retrieval route passed every check once already.


## `verify_safety.py`

**1. verify_safety.py:24** - The trial pins are antigens with many registered trials, so zero is distinguishable from untried.

> *Without it:* Without a pin that must be non-zero, a dead route returns zero for everything and reads as an antigen nobody has tried.

**2. verify_safety.py:34** - Binder records are read from what the binder stage blessed, gated on a hash covering the configuration and not just the gene set.

> *Without it:* Re-deriving costs five minutes and 200 requests. Gating on the gene set alone reuses records built under a different configuration.

**3. verify_safety.py:122** - Counts are checked against the manifest's recorded pool size rather than this run's own input.

> *Without it:* Comparing output to input cannot fail; the specification pins the number independently.

**4. verify_safety.py:190** - The gate reports which of its own questions actually ran.

> *Without it:* A criterion passing on a code path nothing reached has not been tested. The binder stage has already shown what that looks like.

**5. verify_safety.py:108** - S4 is stated against the ceiling applied to each row, with a second clause requiring that anything admitted above the persistent ceiling sits on a route declaring the exposure terminable.

> *Without it:* S4 compared every risk to the persistent ceiling in a design that holds two apart on purpose. It had never seen a target admitted against the terminable one, because none existed in the set it read, and it would have called the first such target a contradiction. Reading the applied ceiling alone would be the opposite error - the record would then justify itself - so the second clause takes the exposure from the routing decision and the ceiling value from the project spec.

**6. verify_safety.py:65** - The gate is given the constructs, so the origin check can see a binder that came from a deposited structure.

> *Without it:* `structure_binders` returns nothing without them, and an adaptor row falls through to NO_GATE with the reason "no binder, so there is nothing to gate" - on a receptor carrying a murine binder. That is a recorded defect, fixed in the service and never in this verifier, and it stayed invisible only because no adaptor row existed in the set this verifier read.



**7. verify_safety.py:214** - S11 re-encodes each construct under a genuinely synonymous codon table and requires the invariant findings to be identical and the map-specific ones to move.

> *Without it:* The basis label is the whole design, and a label nobody can falsify is decoration. Requiring the invariant half to survive re-encoding catches a finding marked invariant but computed from the DNA; requiring the map-specific half to change catches one marked encoding-dependent that is not. The alternate table is synonymous rather than an arbitrary permutation, so the second sequence encodes the same protein and the comparison means what it says.

**8. verify_safety.py:178** - Each detector is tested against a planted known answer and against a control containing none of it.

> *Without it:* A detector that finds everything and one that finds nothing both pass a single-sided check. There is no labelled sequence data here to score against, so a planted answer in both directions is the only test available.

**9. verify_safety.py:262** - S12's positive half is a synthetic domain map, because no shipping design repeats a part.

> *Without it:* All five surviving designs are adaptors and none duplicates a domain, so a criterion phrased only over them would clear on an empty population - the defect recorded as instance 8. The synthetic map supplies the case the pool does not.
## `verify_developability.py`

**0. verify_developability.py:16** - The binder read is gated on the Stage 4 hash, like its siblings, rather than taking whatever the cache slot holds.

> *Without it:* `read_binders` checks the payload against its own digest and nothing else. D5 then compares the rows scored against the records they were scored from, which is a self-comparison that any binder set satisfies, so a set belonging to another configuration would be scored and reported without a single criterion noticing. Recorded in full against `verify_ranking_final.py`, where the same read had the same gap.

**1. verify_developability.py:18** - Binder records are read from the persisted artifact rather than re-queried.

> *Without it:* Re-querying two hundred accessions for every downstream stage is slow and loses the whole run to one dropped connection.


## `verify_package.py`

**1. verify_package.py:65** - Q1 fails on an empty package set rather than clearing on it, and a run with no survivors must report a status instead.

> *Without it:* Every criterion below Q1 reads the package list. On an empty list they all clear, which is the defect this suite has already found in the construct stage, where five criteria reported success over nothing. An empty package set is a real state - no candidate reached the end - and it is reported the way Stage 11 reports it, as a named status.

**2. verify_package.py:163** - Q6 trips when no probe executed, not only when a probe fails.

> *Without it:* A gaps table whose entries all lack probes would pass a criterion that only checks the probes it finds. The positive pin makes a table of unverified assertions fail, which is the difference between a criterion and a formality.

**3. verify_package.py:185** - Q7 requires the conservative-backup refusal to carry counts, not merely to exist.

> *Without it:* "No conservative design exists" with nothing behind it is the same sentence a blank section would produce if someone wrote one. The counts are what make it a measurement: three single-antigen targets recommended and none assembling, no dual assembling because every partner retrieves no binder.

## `verify_api.py`

**1. verify_api.py:94** - The end-state assertion was re-specified: it no longer pins one terminal status.

> *Without it:* It asserted a single no-buildable-construct end state, which encoded the old outcome as an expectation. Once routing sent eight designs to an architecture that fits, the criterion was testing yesterday's answer and would have failed on a correct improvement. What it was protecting is unchanged and is what is checked now: buildable designs are counted, and one carrying a sequence is told apart from one that fits without residues.

**2. verify_api.py:117** - What must hold is the partition, not the value: the chain still partitions the pool exactly, and survivors decompose into complete plus awaiting with nothing unaccounted for.

> *Without it:* Pinning the end state repeats the mistake above - the terminal status is a result, and a criterion that pins a result fails whenever the pipeline gets better.

**3. verify_api.py:133** - The empty-collection guard is written so an empty list cannot raise.

> *Without it:* Indexing a one-element fallback defends only a missing key. An empty list raises IndexError and crashes the check instead of tripping it, so the criterion dies rather than reporting.


## `verify_indications.py`

**1. verify_indications.py:28** - The shared sources are asserted NOT to have gained a per-indication copy.

> *Without it:* Duplicating them would assert that normal-tissue biology changes with the diagnosis. The criterion catches a namespacing fix applied too broadly.

**2. verify_indications.py:33** - The validation-mode pin is a target the platform would not have chosen itself, and the criterion requires it to rank poorly.

> *Without it:* It is the canonical marker for a blood cancer, not a solid tumour, and ranks around 1,300 of 3,400 here. A pin the platform already likes makes a passing verdict indistinguishable from the platform agreeing with itself.

**3. verify_indications.py:115** - The atlas-less criteria call the pipeline end to end instead of inspecting a constructed object or grepping the source.

> *Without it:* Both earlier forms passed while the real behaviour was an AttributeError fifty-five lines before the refusal was ever computed. A criterion that greps for the words it wants is a criterion testing that someone wrote them.

**4. verify_indications.py:158** - An absent rank trips this criterion rather than satisfying it.

> *Without it:* Written the other way round, the criterion whose job is to prove the verdict is not self-agreement passes hardest exactly when the verdict means least.

**5. verify_indications.py:192** - The missing-source criterion is exercised against an indication that really is degraded.

> *Without it:* Otherwise it passes merely because nothing happened to be unavailable, which says nothing about whether degradation is reported.

