# TODO

Work that is known about but not done. Kept in the repository so it survives a
handover: anyone picking up `feature/add_nanopore_support` or starting the
ILLUMINA follow-up should find the open threads here rather than in a chat log.

Sections are ordered by where the work belongs, not by priority. Within each
section, items that could affect analysis results come first.

---

## 1. Needs a Linux machine with Singularity

Cannot be settled on macOS/Docker.

- [x] **Verify `bcftools` works inside the built SIF.** Done on Linux
      (2026-09-18). The SIF is **not** affected: a fresh `apptainer build` from
      `Nanopore_baseContainer.sing` produces a working `bcftools 1.21 / htslib
      1.21` with no `ldconfig` step and no `libhts.so.3` failure, so the `.sing`
      needs no change. Checked beyond `--version`, since that alone would not
      exercise the dynamic link under load: the four subcommands the pipeline
      actually uses — `norm`, `filter -i "FORMAT/AF >= …"`, `index --tbi` and
      `consensus --mask` — were run against `tests/data/bcftools/` and
      reproduced the fixture's expected consensus byte for byte, and a
      malformed VCF still failed non-zero. Those four are the only bcftools
      calls anywhere in the pipeline, so the missing plugin directory in the
      image is immaterial. The whole container suite (14 tests at the time) then
      passed against the rebuilt image.

      Why the Docker build hit it and the SIF does not is worth keeping in
      mind rather than treating as settled luck: both install htslib to
      `/usr/local/lib`, but the `%post` shell in the SIF build leaves a linker
      cache state the Docker layer does not. Any change to the htslib install
      step should re-run the check above.

- [x] **Confirm the truth test passes under Singularity**, not only under
      `-profile docker`. Done on Linux, and repeatedly since — most recently
      2026-09-22 against the camelCase refactor, where it passed in 38.6s
      alongside the multi-sample fan-out test. The Clair3 directive's
      `docker://` handling is fine under Singularity: the digest-pinned image is
      pulled and converted to
      `hkubal-clair3@sha256-1430f7b5….img` in `NXF_SINGULARITY_CACHEDIR`, and
      Clair3 runs from it natively on amd64.

      Note the truth values survived five changes to the nanopore path made
      here — the `--trimLen` normalization, the GENPLOTS staged-path fix, the
      switch to `bam trimBam --clip`, the bamUtil tool spec, and the rename —
      which is the main thing this test is for. The soft-clip change did move
      depth (mean 176.36 → 166.47 on the 5072-read fixture), but the designed
      variants and masked consensus still match.

- [x] **Confirm Singularity/Apptainer pull and run `docker://` images cleanly**,
      including by digest. Done on Linux (apptainer 1.5.3, Singularity 3.8.7
      also installed). Every nanopore run of this review pulled Clair3 from its
      pinned digest and ran it, so this has been exercised repeatedly rather
      than once. Three reference forms are proven, all sitting in the local
      caches:

      | Form | Cached as |
      |---|---|
      | `docker://` by digest | `hkubal-clair3@sha256-1430f7b5….img` (1.4G) |
      | `docker://` by tag | `hkubal-clair3-v1.2.0.img` (1.4G) |
      | third-party registry | `community.wave.seqera.io-library-pip_bio_numpy_pandas-….img` (217M) |

      The digest form is the one that matters: Nextflow converts the OCI image
      to a `.img` in `NXF_SINGULARITY_CACHEDIR` on first use and reuses it after,
      the digest survives into the filename, and `container_manifest.tsv`
      records it as `remote_uri` with checksum `NA`.

      One practical caveat for anyone repeating this: the conversion needs
      several GB of temporary space and Singularity/Apptainer default to `/tmp`.
      On a machine where `/` is tight that fails as a confusing "no space left
      on device" partway through the pull. Set `APPTAINER_TMPDIR`,
      `SINGULARITY_TMPDIR` and `NXF_SINGULARITY_CACHEDIR` somewhere with room —
      both scripts in `viralflow_box/` do this.

      This unblocks the modular-container evaluation below: registry pulls are
      reliable here, so that decision can be made on its merits rather than on
      whether the mechanism works.

- [ ] **Then evaluate breaking the NANOPORE base container into modular pulled
      images.** Today `baseContainer.sif` is one monolith built locally from
      `Nanopore_baseContainer.sing`, carrying Porechop_ABI, minimap2, samtools,
      bcftools and bamUtil. If registry pulls are reliable, each process could
      instead name a small published image the way `run_clair3` already does.
      What it would buy: no local build step, so no `ldconfig`-class surprises
      and no unreproducible SIF (the truth fixture records the base container's
      checksum as `NA` precisely because two people building the same `.sing`
      get different bytes); per-tool version pinning by digest; and a much
      smaller download for anyone who only needs part of the pipeline.
      What to weigh against it: more images to track, registry availability
      becomes a run-time dependency for every process rather than one, and
      airgapped sites would need all of them mirrored. Worth deciding before
      the alpha ships, since it changes what gets published.

- [x] **Confirm `intrahost_analysis:1.1.0.sif` really carries what its recipe
      says.** Done on Linux (2026-09-21). The image matches the expected pins
      exactly, so `ruff.toml` and `INTRAHOST_CONTAINER_PYTHON` need no change:

      | | expected | in the image |
      |---|---|---|
      | Python | 3.8.x | **3.8.0** |
      | pandas | 1.5.3 | **1.5.3** |
      | numpy | 1.23 | **1.23.0** |
      | biopython | 1.81 | **1.81** |

      `bam-readcount` is present too, at
      `/usr/local/bam-readcount/build/bin/bam-readcount`, which matters because
      `runReadCounts` shares this image. `intrahost.py` compiles under the
      container's own interpreter (`python3 -m py_compile`), and the
      `uv run --python 3.8` path the `intrahost-py38` pre-commit hook uses works
      here as well.

      The image provides both `python` and `python3`
      (`/usr/local/bin/mm/bin/`), which is directly relevant to the `fixWGS`
      item in section 4: giving that process this container would satisfy its
      `#!/usr/bin/env python` shebang as it stands, though switching the
      shebang to `python3` is still worth doing rather than relying on a
      `python` alias.

      Confirmed end to end afterwards: `runIntraHostScript` completed for all
      three samples of `viralflow_box/run_illumina_check.sh`.

- [ ] **Decide where the snpEff writable-filesystem workaround belongs.**
      Deleting the commented arch block changed nothing, but it made a
      pre-existing gap visible: `--writable-tmpfs` is now set *only* by
      `-profile singularity`. A plain `nextflow run` (no profile) and
      `-profile apptainer` both get Singularity/Apptainer with no writable
      option at all, and the wrapper's `--profile` defaults to `None`, so no
      profile is the common path. On `develop` the default did supply one via
      the arch switch. This does not affect NANOPORE — nothing in the nanopore
      base container writes inside itself — so it is an ILLUMINA question, but
      it is a real behaviour difference from `develop` and someone should
      confirm it on a box that can actually run snpEff. Three possible answers:
      restore `--writable-tmpfs` as the default in `containers.config`, add it
      to the `apptainer` profile as well, or finally fix snpEff to write
      outside its container, which the GAMBIARRA note has wanted all along.
      Check first whether the target Singularity build even supports
      `--writable-tmpfs`; some HPC builds without overlay support reject it,
      which would rule out making it a global default.

---

## 2. Open review feedback on PR #47

Checked against the live PR via the GraphQL `reviewThreads` API, which reports
resolution state rather than inferring it from reply counts: **28 threads, 20
answered, 0 resolved.**

Nothing nanopore-scoped is left to write. Six threads still need a reply, two
are compliments that need none, and no thread has been marked resolved — worth
doing as each is agreed, so the next reviewer sees what is actually left.

### Still needs a reply on GitHub

Done in code; the thread is just waiting for a note.

- [ ] `docs/parameters.md` and `docs-pt/parameters.md` — the `outDir` `--`
      prefix. The `docs-es/` thread was answered with "sorted"; these two are
      the same fix in the other two languages and were left unanswered.
- [ ] `nanopore_summary.py:56` — the `fileinput.hook_compressed` tip. Taken;
      commit `69213ea`.
- [ ] `getUnmappedReads.nf` — the `> output.gz` question. Answered in full on
      the `getMappedReads.nf` thread; this one needs a pointer to it, since both
      modules were fixed together.
- [ ] `runAmpliconClip.nf` — the formatting suggestion. The module was rewritten
      and is now multi-line, though not in exactly the suggested order.
- [ ] `vfnext/README.md` — the command-formatting suggestion. Applied.

### No reply needed

- Two "nice!" comments, on `.pre-commit-config.yaml` and on reading the version
  from one place in `main.nf`.

### Answered, but the answer has since gone stale

- [ ] **`ILLUMINA.nf` channel naming** was answered with "gonna tidy up that on
      the ILLUMINA work branch" — but it was then done *in this PR*, as part of
      the standardisation agreed on the `NANOPORE.nf` thread. `ILLUMINA.nf` no
      longer mixes `bam_Out_ch`, `bam_output_ch` and `alignCon_Out_ch`. Worth a
      follow-up so the reviewer does not go looking for it on a later branch.
- [ ] **`intrahost.py` auto-formatting** was answered with the `except:` ->
      `except Exception:` change, which is one of *two* behaviour changes ruff
      made. The other is arguably the bigger one: it rewrote the multi-context
      `with` into the parenthesized 3.10+ form, silently raising the script's
      minimum Python and making it unloadable in the very container another
      thread asked us to move to. Both are fixed; the second is worth mentioning
      because the `intrahost_analysis` switch depended on it. That thread also
      asked to drop the `v2` from the filename — already done, the file is
      `vfnext/bin/intrahost.py`, and the reply did not say so.

### Still open — nanopore scope

- [x] **`nanopore_summary.py`: are `masked_bases` and `consensus_n_bases`
      redundant?** No — settled by experiment, not reasoning. They agree only
      because both SARS-CoV-2 references in the repo are N-free; a reference
      carrying `N` makes them diverge by exactly that count. Kept both,
      documented the distinction in the script docstring, `NANOPORE.md` and the
      three parameter tables, and pinned it with an N-containing fixture.
      Commits `50f8701`, `c657d1f`.
- [x] `nanopore_summary.py`: replace the `.gz` if/else with
      `fileinput.hook_compressed`. Commit `69213ea`.
- [x] **`runClair3.nf`: pass the reference index in as an input.** Done via a
      new `run_faidx` process. Commit `d6ee699`.
- [x] `getMappedReads.nf` / `getUnmappedReads.nf`: answer the questions about
      the `-s` removal and whether `> out.gz` already gzips. **Both settled
      experimentally with samtools 1.21.** No, redirection does not compress:
      `> out.fq.gz` produced plain text under a `.gz` name, first bytes `40 73`
      rather than the gzip magic `1f 8b`, and `gzip -t` rejected it. And the
      `-s` removal was a bug fix: on single-end data `-s` receives nothing,
      because those reads are not flagged as paired and so are never
      singletons — `samtools fastq -F 4 -s out.fq.gz` on 3 mapped single-end
      reads wrote a 28-byte empty archive and sent all 3 to stdout, where
      Nextflow files them in `.command.out`. `develop` therefore publishes an
      empty single-end FASTQ. Both branches now use `-0`, which is the flag for
      reads carrying neither READ1 nor READ2, and which compresses from the
      file extension so the separate `gzip` step is gone; output payload
      verified byte-identical to the two-step form. **One follow-up in section
      4.**
- [x] `containers.config`: delete the commented-out architecture-specific
      `runOptions` block rather than leaving it commented, now that engine
      selection is profile-driven. Deleted. Verified as a pure no-op:
      `nextflow config` resolves byte-identically before and after for the
      default and for all five profiles. Left prose in its place recording why
      the arch switch was wrong in the first place — it asked for `--writable`
      on amd64, which needs a sandbox directory, and only `pangolin` and
      `snpeff` are built that way (`build_containers.py` uses
      `--sandbox` despite the `.sif` names); every other container is a pulled
      `.sif`, which Singularity can only open read-only.
- [x] `containers.config`: `runIntraHostScript` uses a remote Wave container
      while every other process uses a local SIF — is that deliberate?
      **It was necessary but wrong, and is now reverted to a local SIF.** Before
      commit `90957f6` (on this branch) the process had no `container` entry at
      all, so with Singularity enabled it could not run. The reviewer's
      suggestion was right: `intrahost_analysis:1.1.0.sif`, already used by
      `runReadCounts`, was built for exactly this pair of processes — its def
      file installs `biopython 1.81`, `pandas 1.5.3` and `numpy 1.23`, which
      bam-readcount has no use for. Switching to it drops the run-time registry
      dependency, the amd64-only limitation and the digest pin in one go.
      One catch, described below, had to be fixed first.
- [x] **`intrahost.py` could no longer run on Python 3.8.** Found while
      checking the above. `ruff format`, with the repo-wide
      `target-version = "py312"`, had rewritten the multi-context `with` at
      line 107 into the parenthesized form, which is **3.10+**; on `develop` it
      is a single line that any Python accepts. Since
      `intrahost_analysis:1.1.0.sif` pins Python 3.8, the reviewer's suggestion
      was impossible until this was undone. Fixed by pinning
      `per-file-target-version` for the file in `ruff.toml` and reverting the
      construct; verified with `uv run --python 3.8 python -m py_compile`, now
      a pre-commit hook. Note this was a latent hazard either way — nothing
      else declares what Python that script needs.
- [x] `vfnext/README.md`: apply the command-formatting suggestion. Applied —
      three stray blank lines inside the NANOPORE code fence, plus the blank
      line before the following `---` that the rest of the file uses.

### Still open — decisions, not code

These need an answer from the team before anything is written.

- [x] **Release strategy** (`NANOPORE.md`) — **decided on the PR.** Modular
      containers are the right direction, but this PR keeps the single base
      container so real-world testing can start. The modular work is the
      section 1 item gated on confirming `docker://` pulls under Singularity.
- [x] **Move `dev.md` into the readthedocs folders.** Moved to
      `docs/development.md` and added to the `toctree` in `docs/index.md`, which
      is what actually puts a page in the sidebar — a Markdown file that no
      toctree lists is built but unreachable, and Sphinx warns about it. The
      docs build is Sphinx with MyST, configured by `.readthedocs.yaml` at the
      repository root. Verified locally with the same Python and requirements
      Read the Docs uses:
      ```bash
      uv run --no-project --python 3.13 --with-requirements docs/requirements.in \
        sphinx-build -b html -W docs /tmp/rtdbuild
      ```
      Two stale claims in the page were corrected on the way: it said the
      pre-push stage runs 14 tests, when it runs 18, and it named three specific
      test files as the container-requiring ones, which has not been true for a
      while. Both are now phrased so they cannot go stale again.
      **Two follow-ups below.**
- [ ] **Decide whether the development page should be translated.** It is in
      `docs/` only; `docs-es/` and `docs-pt/` have their own `index.md`
      toctrees, which still list five pages each. Nothing is broken — their
      builds pass — but a Spanish or Portuguese reader will not find the
      developer docs. Either translate it, or link to the English page from the
      other two indexes.
- [ ] **Work out how `docs-es/` and `docs-pt/` are actually built.**
      `.readthedocs.yaml` names `docs/conf.py` and nothing else, and Read the
      Docs reads one config from the repository root, so nothing in this
      repository explains how the other two trees reach
      `viralflow.readthedocs.io/es/` and `/pt-br/`. Most likely they are
      separate Read the Docs projects configured through its web dashboard,
      which is invisible from here. Worth confirming in the dashboard before
      anyone assumes a change to `docs-es/` will appear online.
- [x] **Drop the `vfnext/` directory** — **agreed on the PR**, deferred to the
      next round of work. It exists only because the Nextflow code was kept
      apart from the wrapper early on. Moved to section 4.
- [x] **An `annotations` subworkflow** — **agreed on the PR**, deferred until
      NANOPORE merges, on the grounds that this PR already carries a lot of
      non-nanopore change. Moved to section 4.
- [x] Should the `coveragePlot` import live in `GENPLOTS.nf`? — **agreed on the
      PR**, deferred to the ILLUMINA branch. Moved to section 4.
- [ ] **The recipes for the pulled containers are no longer in the repo.**
      `def_files/` held one `.def` per image until `c2be157` deleted them in
      favour of pulling prebuilt images from the Sylabs library; only the
      pangolin and snpEff recipes survive. So the only way to learn what is
      inside `intrahost_analysis:1.1.0.sif` — the question this branch just had
      to answer — is `git show 5780713:...`, and nothing guarantees the
      published image still matches that recipe. Worth restoring the def files,
      or recording the build inputs somewhere the pull step can check.
      Relevant beyond bookkeeping: the images pin Python 3.8, 3.9.13, 3.10 and
      one unpinned, and that spread is invisible from the repository.
- [x] **Naming convention**, two comments making related but distinct points —
      both done, standardising on the ILLUMINA style. Worth noting the review
      calls that style PascalCase, but the existing names (`runFastp`,
      `getMappedReads`, `alignConsensus2Ref`) are camelCase; camelCase is what
      the branch now follows throughout.
      - 16 processes renamed from snake_case: the nine nanopore ones, the three
        metadata ones and four test-fixture helpers. No snake_case process name
        is left in the repository.
      - 62 channel and emit names renamed across `ILLUMINA.nf`, `NANOPORE.nf`,
        `GENPLOTS.nf`, `step0-input-handling.nf`, `main.nf` and the test
        fixtures, which also settles the second comment: `ILLUMINA.nf` had been
        mixing `bam_Out_ch`, `bam_output_ch`, `bwaidx_Output_ch` and
        `alignCon_Out_ch` in one file.
      - Process input declarations went with them where they shared the names
        (`ref_fa` -> `refFa`, `ref_gff` -> `refGff`).

---

## 3. NANOPORE bug audit, before the PR merges

From a read-through of the nanopore path on 2026-09-23. Items marked as
confirmed were checked against real output — Clair3 VCFs from the
integration tests, runs of the base image, or upstream documentation; the rest
come from reading the code and want a run before anyone relies on them.

### Can change analysis results

- [ ] **Clair3 `LowQual` calls reach the consensus, so `clair3_qual` does
      nothing to it.** Confirmed from Clair3's own documentation: `--qual` does
      not drop variants, it labels them `PASS` or `LowQual` and keeps both.
      `runBcftools` filters on `FORMAT/AF` alone (`modules/runBcftools.nf`,
      the `bcftools filter` line) and never looks at `FILTER`, so a `LowQual`
      call with AF >= `af_threshold` is written into the consensus. The value
      of `clair3_qual` ends up only in the summary TSV. `NANOPORE.md` already
      says no `FILTER=PASS` condition is applied, but not that this makes the
      parameter inert. Fix: add `-f PASS` (or `-i 'FILTER="PASS" && …'`) and
      a fixture VCF with a `LowQual` row.
      An earlier draft of this audit also claimed `RefCall` rows would be
      applied, because Clair3 defines their `AF` as the *reference* allele
      frequency and `bcftools consensus` is run without `-s`. That part is
      **wrong for bcftools 1.21**: on a single-sample VCF it uses the sample's
      GT anyway (it prints `applying IUPAC codes based on FORMAT/GT in sample
      sample`), so a `0/0` row is not applied. RefCall rows are also off by
      default in Clair3 (`--print_ref_calls`). The same message does raise a
      question worth one fixture: with `--haploid_sensitive`, Clair3 can call
      `0/1`, and bcftools may then write an IUPAC ambiguity code rather than
      the ALT. Nothing tests that today.
- [x] **A sample with no aligned reads published the reference as its
      consensus, reported 100% callable.** Confirmed end to end, then fixed.
      `runBcftoolsConsensus` built its mask from `samtools depth -J -a`, and a
      single `-a` prints zero-depth positions only on contigs with at least one
      read: a contig no read reached is not in the output at all. For a
      negative control (20 random 1 kb reads, none aligned) the run succeeded,
      `low_cov.bed` was empty, the consensus was byte-identical to
      `NC_045512.2`, and the summary said `zero_depth_bases 0`,
      `callable_percent 100.000000`. The same gap left any uncovered contig of
      a multi-contig reference unmasked. Now `-aa`, which reports every
      reference position. Covered by two new cases in
      `tests/workflows/bcftools-fixture.nf.test` (no reads; a two-contig
      reference with one contig unreached) and by
      `integration_tests/nanopore-no-reads.nf.test`, which runs the whole
      NANOPORE workflow, Clair3 included, on a negative control. All three
      fail against `-a`. The truth and multi-sample tests are unchanged by it.
- [x] **An empty FASTQ aborted the whole batch.** Found while testing the item
      above, then fixed by rejecting it in step0. A barcode that demultiplexed
      nothing gives a 20-byte gzip with no reads, which passed input
      validation (it rejected only zero-byte files). Porechop_ABI's ab initio
      adapter inference then failed with `ERROR - Unable to build graph`,
      exit 1, and with the default error strategy every other sample in the
      run was abandoned with it. `validateFastqPath` now reads each input up to
      its first non-blank line and reports `FASTQ contains no reads` for one
      that has none, gzipped or plain, alongside every other input problem
      and before any task is submitted. A `.gz` that is not gzip is reported
      as `FASTQ cannot be read as gzip` instead of failing later inside
      `prepareSampleReads`. It applies to both modes and to `--inDir` as well
      as `--samplesheet`, and it also catches a `concat-fastq` output whose
      reads were all dropped by the length filter. Covered by
      `Rejects FASTQ files that hold no reads` in
      `tests/workflows/input-fixture.nf.test`.

      The check is per file, like the zero-byte check it sits beside, so a
      multi-chunk sample with one read-less chunk is rejected even though its
      other chunks hold reads. The alternative considered was to let a
      read-less sample through to an all-`N` consensus, which the fix above
      makes correct, by skipping Porechop when there is nothing to trim.
      Rejecting was chosen as the simpler, explicit behaviour; revisit if
      empty barcodes turn out to be routine in real runs. For ILLUMINA, where
      empty samples used to be reported rather than rejected, this is a
      behaviour change. See the section 4 item "A read-less sample now stops
      an ILLUMINA run".
- [ ] **Primer clipping trims only one end of each read.**
      `modules/runAmpliconClip.nf` runs `ampliconclip --strand` without
      `--both-ends`, so only the 5′ primer is clipped. Nanopore amplicon reads
      span the whole amplicon, so the 3′ primer stays in the consensus.
      ILLUMINA's `ampliconclip.nf` does pass `--both-ends`. The fixture has one
      left primer on one forward read, so it cannot see this; it needs a read
      carrying both primers. Two smaller points on the same line: `--strand`
      needs a 6th strand column that nothing validates, and there is no
      `--filter-len`, so reads clipped to nothing are kept.
- [ ] **The mapped-reads FASTQ carries duplicate and truncated reads.**
      `samtools fastq -F 4` in `modules/getMappedReads.nf` *replaces* the
      default exclusion mask `0x900`, so secondary and supplementary alignments
      are written as reads too. minimap2 on ONT data produces many
      supplementary records, hard-clipped fragments sharing the primary's name.
      Fix: `-F 0x904`. ILLUMINA shares the module and the problem.
- [ ] **A multi-contig reference gives duplicate FASTA headers.**
      `runBcftoolsConsensus` renames every header to `>${meta.id}` with `sed`,
      and NANOPORE never checks the reference is a single sequence;
      `coveragePlot` likewise draws only `references[0]`. Either reject
      multi-contig references in step0 or name headers `${id}|${contig}`. The
      new two-contig bcftools fixture asserts sequence lines only, so it will
      not pin the current behaviour.

### Documented runs that fail, or record the wrong thing

- [x] **A full NANOPORE run needs ILLUMINA images nobody is told to get.**
      GENPLOTS ran in NANOPORE mode in two ILLUMINA images: `coveragePlot` in
      `generate_plots:2.0.0.sif`, `getMappedReads`/`getUnmappedReads` in
      `generate_consensus:2.0.0.sif`. `runFaidx.nf` explains exactly why
      sharing that image breaks nanopore, and GENPLOTS did it anyway. Under
      `-profile docker` there was no way to run them at all: those images are
      SIFs pulled from the Sylabs library (`library://wallaulabs2/…`), with no
      Docker equivalent. In NANOPORE mode all three now use the base image
      (the mode check is in their `configs/containers.config` closures), and
      the base image gained what they need. Both recipes:
      - pin `bamdash==0.4.4` (the version `generate_plots` carries), with
        `kaleido==0.2.1`, `plotly==5.24.1` and `pysam==0.23.3`;
      - add `python-is-python3`, since `coveragePlot` runs under
        `#!/usr/bin/env python`.
      samtools was already there. The images grow: the Docker image from
      1.3 GB to 2.48 GB and the SIF from 460 MB to 778 MB, mostly kaleido's
      bundled Chromium. A NANOPORE run now needs only the base image and
      Clair3, so the manifest lists just those two again, and bamdash is
      recorded in `software_versions.tsv`. Checked on 2026-09-23 under both
      engines: `main-nanopore.nf.test` passes under Singularity and asserts the
      HTML, PNG and SVG plots. Under Docker, GENPLOTS completes in
      `viralflow/nanopore-base:2.0.0a1` and publishes all five GENPLOTS files.
      The PNG and SVG are new: no run of either mode had produced them before
      (see the section 4 item on ILLUMINA's coverage plots).
- [ ] **`-profile docker` through `main.nf` fails in METADATA.** Confirmed on
      the Linux box on 2026-09-23 with the base image built from
      `nanopore_base.Dockerfile` and Clair3 pulled by digest. Run on the truth
      fixture via `tests/data/input/nanopore-truth.csv`. The analysis itself is
      correct under Docker: variants equal `truth.vcf`, every metric matches
      `expected_metrics.tsv`, and the consensus is byte-identical to the
      Singularity run. Rerun with `process.errorStrategy = 'ignore'` to see
      every failure rather than the first, it showed exactly the three
      predicted causes and no others:

      | Task | Exit | Cause |
      |---|---|---|
      | `METADATA:captureToolVersion (clair3)` | 125 | raw `docker://hkubal/clair3@…` reference; `docker: invalid reference format`. `runClair3` strips the prefix for Docker; this does not |
      | `METADATA:captureContainerMetadata (nanopore_base)` | 1 | recorded as `local_sif`, so the tag becomes the path `<launchDir>/viralflow/nanopore-base:2.0.0a1`: `Configured SIF file does not exist` |
      | `GENPLOTS:coveragePlot`, `getMappedReads`, `getUnmappedReads` | 125 | Docker asked to run the ILLUMINA `.sif` paths. **Fixed**: see "A full NANOPORE run needs ILLUMINA images" above |

      The first two remain, and with them the next item, which the same run
      shows as `profile docker`, `container_engine singularity`. Under
      `--profile docker`, as CI runs it, the unit suite passes (46/46), as do
      the truth, multi-sample and no-reads integration tests.
      `main-nanopore.nf.test` fails, because the first METADATA failure stops
      the run. `NANOPORE.md` documents this exact command.
- [x] **`container_manifest.tsv` omits images a NANOPORE run used.**
      `containerSpecs()` listed only `nanopore_base` and `clair3`, not
      `generate_plots` or `generate_consensus`, the drift `containers.config`
      says the shared map prevents. It now declares `generate_plots` for every
      NANOPORE run and `generate_consensus` when mapped reads are written, both
      read from `params.illumina_containers`. The second is gated by
      `writeMappedReadsEnabled()` in `modules/param_helpers.nf`, the same rule
      GENPLOTS uses to decide whether the two processes run, so the two cannot
      disagree. Caught by `integration_tests/main-nanopore.nf.test`, which
      checks every image the trace says a task ran against the manifest.
      Superseded the same day: GENPLOTS now runs in the base image in NANOPORE
      mode (see "A full NANOPORE run needs ILLUMINA images" above), so a
      NANOPORE manifest lists `nanopore_base` and `clair3` only, which is again
      the truth. The same test guards it either way.
- [ ] **`run_manifest.json` records Docker runs as `singularity`.**
      `container_engine` is guessed from the profile name
      (`metadata_helpers.nf`); `workflow.containerEngine` has the real answer.
- [ ] **The wrapper cannot set any NANOPORE parameter.** Neither the
      `parse_params` allow-list nor `viralflow run` knows `clair3_model`,
      `np_min_depth`, `af_threshold`, `clair3_qual`, `clair3_chunk_size`,
      `base_container` or the per-tool cpus/memory, and a params file naming
      one is rejected. The model matters most: the default
      `r941_prom_sup_g5014` is for R9.4.1 flowcells, so an R10.4.1 user of the
      wrapper gets the wrong model with no way out.
- [ ] **The Docker image's smoke test cannot fail the build.** The last `RUN`
      of `nanopore_base.Dockerfile` ends `… && bam help > /dev/null 2>&1 ||
      true`; the `|| true` binds to the whole `&&` chain, so a broken minimap2,
      samtools or bcftools still builds. Only `bam help` needs the exemption.

### Smaller

- [ ] **`np_min_depth` and ILLUMINA's `--depth` mean different things.** Depth
      <= 20 is masked, so nanopore needs 21x; ILLUMINA's `--depth 25` (ivar
      `-m`) needs 25x. Documented, but two similarly named thresholds with
      opposite edge semantics invite mistakes.
- [ ] **The coverage plot draws the wrong threshold for nanopore** —
      `params.depth` (25), while masking uses `np_min_depth` (20).
- [ ] **Silent failures in the GENPLOTS steps.** `getMappedReads` and
      `getUnmappedReads` have no `set -o pipefail`, so a failed
      `samtools sort` still publishes an empty FASTQ; `coveragePlot` calls
      bamdash through `subprocess.run(..., shell=True)` without checking the
      exit status. That is how ILLUMINA's PNG and SVG coverage plots have gone
      missing on every run without anyone noticing; see the section 4 item.
- [ ] Minor: `runPorechop` publishes an uncompressed `*.chopped.fastq`,
      roughly doubling storage; `runNanoporeSummary` calls
      `${projectDir}/bin/…` rather than relying on `bin/` being on `PATH`,
      which breaks on cloud executors; `concat-fastq` defaults to
      `--max-len 500`, silently dropping nearly every read of 1200 bp or
      whole-genome protocols; legacy `--inDir` discovery fails a nanopore file
      named `*_R1.fastq` as an "orphan Illumina mate".

### Test gaps behind these

- Nothing runs `main.nf --mode NANOPORE`, which would have caught the METADATA,
  GENPLOTS and manifest items above.
- The truth fixture's error-free synthetic reads never produce a `LowQual`
  call, a right-hand primer, or a supplementary alignment. The zero-coverage
  case now has its own integration test, and read-less inputs are rejected
  and tested in step0.
- `integration_tests/nanopore-multisample.nf.test` repeats the Clair3 digest in
  its `params` block, and `tests/test_container_recipes.py` checks only the
  truth test's copy, so that one can drift unnoticed. It could simply inherit
  the pin from `tests/nextflow.config`, as `nanopore-no-reads.nf.test` does.

---

## 4. ILLUMINA follow-up branch

Deliberately kept out of the nanopore PR.

- [ ] **Move `main.nf`, `workflows/` and `modules/` to the repository root**,
      dropping `vfnext/`. Agreed on PR #47. Touches every `includeConfig` and
      `$projectDir` path, the wrapper's `root_path`, `.readthedocs.yaml` and the
      CI workflow, so it wants its own PR with nothing else in it.
- [ ] **ILLUMINA coverage plots have never included the PNG or SVG.**
      `coveragePlot` asks bamdash for HTML, PNG and SVG, but only the HTML is
      ever published, in every run checked, including the 2026-08-12 baseline
      under `viralflow_box/test_box/output/`. The `generate_plots:2.0.0.sif`
      image has kaleido 1.0.0, while bamdash 0.4.4 requires `kaleido==0.2.*`.
      Kaleido 1.x needs a separately installed Chrome for static export, the
      image has none, and the task log says so (`plotly_get_chrome`, then
      `mv: cannot stat 'NC_045512.2_plot.png'`). The task still succeeds,
      because `coveragePlot` ignores bamdash's exit status (the "Silent
      failures in the GENPLOTS steps" item in section 3). The nanopore base
      image pins `kaleido==0.2.1` with `plotly==5.24.1` and produces all three
      formats (checked on the truth fixture), so rebuilding `generate_plots`
      with the same pins is the likely fix. Installing Chrome is the
      alternative. Inherited: the image predates this branch.
- [ ] **Publish through workflow outputs, and retire `--outDir` in favour of
      Nextflow's `outputDir`.** `outputDir` (`-output-dir`) is the root for the
      workflow `output {}` block, which Nextflow documents as "intended to
      replace the publishDir directive" (stable since 25.10; we require 26.04).
      It does not affect `publishDir`, and every module here publishes through
      `publishDir` under `params.outDir` (28 files), so switching the parameter
      alone gains nothing. Today `nextflow.config` sets
      `outputDir = params.outDir`, and step0 stops the run when the two
      disagree. What migrating buys: one output-directory setting instead of
      two, with no conflict check; index files, a structured per-sample
      catalog of what was published, which could replace hand-written files
      like `resolved_sample_inputs.tsv`; and data lineage (next item), which
      records only outputs published through `outputDir`. It does not move the
      trace, report or timeline: those default to the launch directory, not
      `outputDir`, so their paths stay explicit in `nextflow.config`, and
      `run_manifest.json` records where they actually went. Keep `--outDir` as
      a deprecated alias that warns and sets `outputDir`, since the wrapper,
      the three language docs and users' parameter files all use it, and drop
      it in v3 along with `--inDir`. Touches every module, so pair it with the
      repository-root move above.
- [ ] **Capture provenance in the tasks that ran, instead of predicting it.**
      Nearly every metadata bug on the nanopore branch has had the same cause:
      `containerSpecs()` and `toolSpecs()` in `modules/metadata_helpers.nf`
      predict which images and tools a run will use, and the prediction drifts.
      Examples: the GENPLOTS images missing from the manifest, bamUtil needing
      its own gate, version capture using the raw `docker://` Clair3 reference
      under `-profile docker`, ILLUMINA's `toolSpecs()` still hardcoding `.sif`
      paths rather than reading `params.illumina_containers`, and bamdash
      missing from `software_versions.tsv` until GENPLOTS moved into the
      nanopore base image. Instead, have each process
      report what it used, in the same task:
      ```groovy
      tuple val(task.process), val(task.container), eval('samtools --version | head -n 1'), topic: provenance
      ```
      METADATA then reads `channel.topic('provenance')`, de-duplicates, and
      checksums each distinct local image once. Checked on 26.04.6 in an
      ordinary (untyped) process: it emitted the process name, the real
      container path and `samtools 1.21` from inside that container. Topic
      channels are stable since 25.04, `eval` since 24.04. This is the pattern
      nf-core moved to for versions. It removes the separate
      `captureToolVersion` tasks and the Docker METADATA failure by
      construction. Costs: an output line in every process, ILLUMINA's
      included, which is why it is not on the nanopore branch; and a failing
      version command fails the task, so each needs care.

      Not a replacement for our files: the trace, report and timeline cannot
      be read to build the manifest, because `workflow.onComplete` runs before
      Nextflow writes the report and timeline (`Session.shutdown0()` runs the
      shutdown callbacks, then notifies the observers). The trace is written
      asynchronously, nf-test replaces it, and users can move or disable all
      three. Keep `run_manifest.json` for the analysis settings Nextflow knows
      nothing about. Data lineage (`lineage.enabled = true`) records
      parameters, commit, per-task container and input/output checksums, but
      it is experimental ("may change in future releases"), records only
      outputs published through `outputDir`, and writes no new records for
      cached tasks. Worth trying alongside our files once publishing moves to
      workflow outputs.
- [ ] **Add an `annotations` subworkflow** grouping snpEff, pangolin, nextclade
      and compileOutput. Agreed on PR #47. The reviewer's point is that all
      three could serve nanopore output too, so this is what would let NANOPORE
      reuse the ILLUMINA annotation stack rather than reimplement it.
- [ ] **Move the `coveragePlot` import into `GENPLOTS.nf`.** Agreed on PR #47.
- [ ] **`fixWGS` fails for every sample and always has.** It runs its script
      under `#!/usr/bin/env python` (`modules/fixWGS.nf:17`) and is the only
      ILLUMINA process with no entry in `configs/containers.config`, so it runs
      on the host — where Ubuntu 24.04 provides `python3` and no `python`:
      ```
      ILLUMINA:fixWGS (ART1)  exit: 127
      /usr/bin/env: 'python': No such file or directory
      ```
      The run does not stop. `compileOutputs` still writes a batch summary and
      the only trace is the task status, which is how this went unnoticed.
      Whatever the step contributes has therefore never been produced, and it
      is worth establishing what that is before fixing it. Seen on all three
      samples of `viralflow_box/run_illumina_check.sh`, and identically in the
      2026-08-12 baseline under `viralflow_box/test_box/output/`, so nothing
      about it is new. Inherited from `develop`, where the same shebang and the
      same absent container directive are already present. The fix is a
      container plus `python3`: the script imports `pandas` and `Bio`, which
      `intrahost_analysis:1.1.0.sif` already carries for `runReadCounts` and
      `runIntraHostScript`.
- [ ] **Boolean parameters given on the command line are ignored, and the
      metadata layer disagrees with the workflow about them.** Nextflow hands
      over command line parameters as Strings, and the two idioms that read
      them are each wrong for a String in a different direction:

      | Idiom | Where | `--runSnpEff true` | `--runSnpEff false` |
      |---|---|---|---|
      | `params.X == true` | `workflows/ILLUMINA.nf:112`, `modules/param_helpers.nf:35` | **false** | false |
      | `if (params.X)` | `modules/metadata_helpers.nf:269`, `:315` | true | **true** |

      `param_helpers.nf:35` is `writeMappedReadsEnabled()`, which GENPLOTS and
      `containerSpecs()` both call, so for `--writeMappedReads` the workflow and
      the metadata already agree. Both still reject the String `"true"`, so
      normalizing there fixes both at once.

      So `--runSnpEff true` **silently does not run snpEff** — `"true" == true`
      is false in Groovy, and only the `nextflow.config` default, a real
      boolean, ever enables it. Same for `--writeMappedReads true` and
      `--dedup true`. Meanwhile `container_manifest.tsv` and
      `software_versions.tsv` declare `snpeff` and `generate_report` for a
      `--runSnpEff false` run, because the metadata layer reads the same
      parameter as truthy; confirmed in
      `viralflow_box/test_box/output_illumina_check/RUN_METADATA/`.

      The `== true` half is inherited from `develop` (`main.nf:164` and `:181`
      there; this branch only moved the gates into `GENPLOTS.nf` and
      `ILLUMINA.nf`). The disagreement is not: the metadata layer is new here,
      so a provenance record that contradicts the run is this branch's to own.
      `--writeMappedReads` also gates NANOPORE through `GENPLOTS`, so this is
      not purely an ILLUMINA concern — fixing that one instance before the PR
      merges is defensible.

      Same class as the `--trimLen` String/Integer bug fixed in `074caf3`. The
      fix is the same shape: a `normalizeFlag(value)` beside `normalizeTrimLen`
      in `modules/param_helpers.nf`, rejecting anything that is not a
      recognised boolean, read by the workflow gate and the metadata gate
      alike. Auditing for the remaining `params.X == true` and bare
      `if (params.X)` sites is part of the job.
- [ ] **A read-less sample now stops an ILLUMINA run instead of being
      reported.** A behaviour change this branch introduces, not a bug, but
      nobody has decided it for ILLUMINA. Before this branch, ILLUMINA let an
      empty sample through and reported it: the 2023 fixture in
      `viralflow_box/empty_fastq_set/` (ART1 plus a sample with zero-byte R1
      and R2, run by its `run_test.sh`) finished ART1 and wrote
      `sample,ERROR,No consensus sequence obtained` to
      `COMPILED_OUTPUT/errors_detected.csv`. On this branch the same input
      stops at input validation, before any task is submitted, and ART1 is
      never processed. Confirmed on the Linux box on 2026-09-23, with ART1
      next to a zero-byte sample and next to a 20-byte gzip sample:
      ```
      Input validation failed:
       - legacy sample sample fastq_1: FASTQ is empty: .../sample_R1.fastq.gz
       - legacy sample sample fastq_1: FASTQ contains no reads: .../sample_R1.fastq.gz
      ```
      (first line from the zero-byte run, second from the 20-byte run;
      `completed=0` both times).

      Two commits, both on this branch: the zero-byte check came with
      `4a80659` (samplesheet handling) and is absent from `develop` and
      `main`. `28d2fd4` extended it to read-less files, which is the section 3
      item "An empty FASTQ aborted the whole batch". That extension is what
      makes this likely to happen: a sample that demultiplexed to nothing
      usually arrives as a 20-byte gzip, not a zero-byte file. Low-read
      samples still pass. `test_box`'s Cneg negative control has 63 pairs and
      runs as before.

      For NANOPORE, rejecting is clearly better than what it replaced:
      Porechop_ABI crashed and took every other sample down with it. For
      ILLUMINA, a plate with one failed library no longer gives partial
      results, and the operator has to remove the sample and rerun. Options:
      - Keep rejecting in both modes, and say so in the release notes and
        `docs/parameters.md` so users are not surprised.
      - In ILLUMINA, warn and let read-less samples through to the existing
        `No consensus sequence obtained` path. This needs the pre-branch path
        re-checked with a 20-byte gzip: the 2023 fixture only covers zero-byte
        files.
      - Warn and drop read-less samples from the batch in both modes, keeping
        a record of each in the compiled output.

      Worth deciding before the PR merges, even if the code change waits for
      this branch: the current behaviour ships with the merge.
- [ ] **`getMappedReads.nf` / `getUnmappedReads.nf`: the paired branch
      desynchronizes R1 and R2.** Neither branch passes `-s`, so a read whose
      mate was removed by the `-F 4` / `-f 4` filter is written to the R1 file
      rather than set aside as a singleton, and the two files end up with
      different numbers of records in different orders. Demonstrated with a
      four-record BAM: without `-s`, R1 held `orphan_b` and `pair_a` while R2
      held only `pair_a`; adding `-s` put `pair_a` in both and `orphan_b` in a
      singleton file. Anything that reads the published FASTQs as a matched
      pair gets mismatched mates. This is inherited from `develop`, not caused
      by this branch, which is the only reason it is deferred — it produces
      wrong data. The fix is `-s ${meta.id}.mapped.singleton.fq.gz`; the
      existing `*.mapped.*.fq.gz` glob already picks the new file up, but it
      changes what a run publishes, so it wants a test alongside it.
- [ ] **`compileOutput.py`: `mepf_reads_aligned` is a typo** for
      `pf_reads_aligned`, and it becomes a **column header in the published
      `reads_count.csv`**. Confirmed against picard 2.27.2's
      `AlignmentSummaryMetrics` source that index 5 is `PF_READS_ALIGNED`. The
      cause is implicit string concatenation — `"me"` and `"pf_reads_aligned"`
      on adjacent lines, which Python joins silently. Renaming changes a
      published column, so it belongs in a release note.
- [ ] **`compileOutput.py` has no tests at all.** It is how the wrong column
      name shipped unnoticed. Cover at least `__parse_metrics`,
      `__parse_wgs` and the two `virus_tag` branches.
- [ ] **`compileOutput.py`: a missing `<cod>.depth<N>.fa.bc` crashes with a raw
      `FileNotFoundError`.** That file is not in the checked-files list, so
      unlike every other per-sample input its absence is an unhandled traceback
      instead of the usual "missing output" warning plus skip.
- [ ] **`compileOutput.py`: `get_lineages_summary` reads `./wgs.csv` from the
      working directory** while `compile_output_fls` writes it to
      `--outputDir`. These agree only because the Nextflow process passes
      `-oD ./`. Any other output directory silently skips the lineage summary.
- [ ] **`intrahost_scriptv2.py` has no tests either**, and the reviewer suggests
      renaming it to `intrahost.py` (dropping the `v2`).
- [ ] **Move `intrahost.py` onto a current Python.** It is pinned to 3.8 today
      only because `intrahost_analysis:1.1.0.sif` is, and 3.8 has been
      end-of-life since October 2024. Doing this properly means rebuilding that
      container on a supported Python — it also holds bam-readcount for
      `runReadCounts`, so both processes move together — then dropping the
      `per-file-target-version` pin from `ruff.toml`, the `intrahost-py38`
      pre-commit hook, the `INTRAHOST_CONTAINER_PYTHON` constant and the
      parenthesized-`with` guard in `tests/test_container_recipes.py`, and
      letting ruff reformat the file. Worth pairing with the rest of this
      section, since the script needs tests before anyone changes it with
      confidence. The pins exist to stop the formatter silently outrunning the
      container again; they are scaffolding, not the goal.
- [ ] **`intrahost.py` has two invalid escape sequences**, at lines 63-64:
      `re.sub(".*\/", ...)` and `re.sub("\..*", ...)`. Running the script
      inside its container surfaced them as `SyntaxWarning` under Python 3.14,
      and they become a `SyntaxError` in 3.15. Ruff is not catching them, so
      check whether `W605` is enabled. Raw strings (`r".*/"`, `r"\..*"`) fix
      both. Worth doing sooner than the rest of this section: the container the
      script runs in is unpinned upstream, so the Python it gets is whatever
      Wave last built against.
- [ ] **Document the six remaining parameters**: `snpEffDBCatalog`,
      `databaseDir`, `nxtclade_jobs`, `pangolin_threads`, `queue_size`,
      `minBamSize` — absent from all three parameter tables.
- [ ] `compileOutput.py`: delete the commented-out debug blocks at roughly
      lines 210, 227 and 320 (reviewer's suggestion).

---

## 5. Repository maintenance, after the PR merges

- [ ] **History cleanup.** Two 4.38 MB blobs live only on this branch —
      `test_files/nanopore/test.fastq.gz` and the `test.fastq.tar.gz` it
      replaced. `test.fastq.gz` was deleted from the working tree, but the
      blobs remain and enter `develop` on merge, because the project merges with
      merge commits rather than squashing.

      Removing them was deliberately **not** done during the PR: it needs a
      force-push, and PR #47 carries 29 live inline review comments that a
      rewrite would risk detaching. It also strips the GPG signature from
      `8f27ef0`, a GitHub-generated merge commit authored by @dezordi.
      Authorship itself is preserved — only the "Verified" badge is lost.

      If this is ever done, do it once for the whole repository: the ~50 MB of
      ILLUMINA test FASTQs already in `main`
      (`test_files/sars-cov-2/input/ART*.fq.gz`) dominate the 58 MB history far
      more than these two. Use `--prune-empty` with care — it also removes
      commits that were already empty, such as
      `Fix exec command in Singularity_snpEff`.

- [ ] **Revisit the pinned runner image.** CI runs on `ubuntu-24.04` rather
      than `ubuntu-latest`, pinned on 2026-09-22 because `ubuntu-latest`
      migrates to Ubuntu 26 on 2026-10-19 and an unannounced base-image change
      during review is not worth the surprise. Someone should move it
      deliberately once Ubuntu 26 has settled. `.readthedocs.yaml` pins the same
      way, so the two are consistent.
- [ ] **`astral-sh/setup-uv` is pinned to an exact version** (`v10.2.0`) while
      every other action uses a moving major tag, because astral-sh stopped
      publishing bare major tags after `v7` — `@v10` does not resolve and the
      run fails. This one needs a manual bump; the others do not.
- [x] **Build the docs in CI.** Added as a `docs` job with a matrix over
      `docs`, `docs-es` and `docs-pt`, `fail-fast: false` so one tree failing
      does not mask the others. It runs `sphinx-build -W --keep-going`, and
      reads the Python version out of `.readthedocs.yaml` rather than repeating
      it, the same way the container job reads its image tag from
      `profiles.config`. Mutation-tested: dropping `development` from the
      toctree and restoring the deprecated `display_version` option each fail
      the build. The command is documented in `docs/development.md` so the
      local and CI invocations are the same one.
- [ ] **Add a CI status badge to `README.md`** once the workflow has run on
      `develop`.

---

## Notes for whoever picks this up

- `tests/test_container_recipes.py` guards the version pins that appear in more
  than one file. If you change a pinned tool version, it will tell you which
  other files need the same change.
- The container tests need an image. With Singularity, build
  `containers/baseContainer.sif` from `Nanopore_baseContainer.sing`. Without it,
  build `nanopore_base.Dockerfile` and pass `--profile docker`; see
  `docs/development.md`.
- Clair3's published image is amd64 only, so it runs under emulation on Apple
  Silicon. The workflow still completes, the truth test included — just slowly.
