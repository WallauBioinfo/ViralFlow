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

- [ ] **Build the SIF with bamdash's pinned dependencies.** The Docker image was
      rebuilt and tested with the full closure pinned, `--no-deps` and
      `pip check` (section 3), but a SIF cannot be built on macOS, so
      `Nanopore_baseContainer.sing` has only had its pins checked by
      `tests/test_container_recipes.py`. Two things only a real build shows:
      that the pinned set installs on amd64 as it does on arm64, and that
      `pip check` passes there. The second is the one to watch — it checks the
      whole environment, and the `.sing` installs apt packages with their
      recommendations where the Dockerfile uses `--no-install-recommends`, so
      its system Python packages differ. Then run `main-nanopore.nf.test` under
      Singularity, which asserts all three coverage plots.

      While building, check that `%test` can actually fail — the Docker smoke
      test could not (section 3). Its commands have no pipes, but every tool
      check is followed by an `echo`, so if Apptainer runs `%test` without
      `-e` only the last `echo` decides. Apptainer's documentation says a
      failing command halts the build without naming the shell flags. The quick
      test: add `rm -f /usr/local/bin/bcftools` at the end of `%post` and
      confirm the build fails at `%test`.

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
- [x] **bamdash's own dependencies were unpinned.** Follow-up to the item
      above, found when the Docker image was rebuilt on Apple Silicon
      (2026-09-24). bamdash 0.4.4 asks only for `pandas>=1.4.4` and
      `biopython>=1.79`, so pip took the newest of everything each time the
      image was built: that build resolved pandas 3.0.6, numpy 2.5.3 and
      biopython 1.88 — pandas two majors past anything bamdash 0.4.4 was
      written against. It worked, but the next build could have pulled a
      release that did not, on any platform. Both recipes now pin the whole
      closure at those resolved versions — `pandas`, `biopython`, `numpy`,
      `python-dateutil`, `six`, `tenacity`, `packaging` — beside the four
      already pinned, all eleven in `PINNED` in `tests/test_container_recipes.py`
      so the recipes cannot drift. They are installed with `--no-deps` and
      followed by `pip check`, so nothing outside the list can enter the
      image, and a bump that brings in a new dependency fails the build
      instead of fetching it unpinned (checked: with `six` removed,
      `pip check` exits 1). Verified on arm64 under `-profile docker`: the
      pinned build installs exactly the same eleven packages, all with native
      aarch64 wheels, kaleido's bundled Chromium exports PNG and SVG, the
      coverage plot shows the truth fixture's gap and deletion correctly, and
      the unit and integration suites pass. The `.sing` side still needs a
      build on Linux; see section 1.
- [x] **`-profile docker` through `main.nf` fails in METADATA.** Confirmed on
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
      | `METADATA:captureToolVersion (clair3)` | 125 | raw `docker://hkubal/clair3@…` reference; `docker: invalid reference format`. `runClair3` strips the prefix for Docker; this does not. **Fixed** |
      | `METADATA:captureContainerMetadata (nanopore_base)` | 1 | recorded as `local_sif`, so the tag becomes the path `<launchDir>/viralflow/nanopore-base:2.0.0a1`: `Configured SIF file does not exist`. **Fixed** |
      | `GENPLOTS:coveragePlot`, `getMappedReads`, `getUnmappedReads` | 125 | Docker asked to run the ILLUMINA `.sif` paths. **Fixed**: see "A full NANOPORE run needs ILLUMINA images" above |

      Both METADATA causes came from the metadata layer ignoring the engine.
      `containerSpecs()` and `toolSpecs()` in `modules/metadata_helpers.nf`
      now take it as an argument; `main.nf` passes `workflow.containerEngine`.
      Under Docker:
      - both NANOPORE images are recorded as a new kind, `docker_image`, by
        reference;
      - version commands run in the reference Docker accepts
        (`engineImageReference()`, the rule `runClair3` applies);
      - `captureContainerMetadata` records the image ID and size from
        `docker image inspect`, pulling first when the image is absent,
        because nothing orders it after the tasks that make Docker fetch it.

      Singularity and Apptainer runs are unchanged. The same `main.nf` run
      now completes with 22 tasks and no failures and no override, and the
      manifest's IDs and sizes match `docker image inspect`. Pinned in
      `tests/workflows/metadata-fixture.nf.test`, whose fixtures now fix the
      engine per test so they assert the same thing under either profile.
      With the next item fixed too, `main-nanopore.nf.test` passes under
      `--profile docker`, as it does under Singularity.
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
- [x] **`run_manifest.json` records Docker runs as `singularity`.**
      `container_engine` was guessed from the profile name
      (`metadata_helpers.nf`): `apptainer` if the name said so, otherwise
      `singularity`. It now records `workflow.containerEngine`, the engine
      Nextflow actually used, or null when none is enabled.
      `main-nanopore.nf.test` checks it against the engine a task's
      `.command.run` launched, and passes under both `--profile singularity`
      and `--profile docker`. `tests/main.metadata.nf.test` also checks it on
      every push and in CI, and fails under Docker with the old guess.
      `executor` was guessed the same way (`profile.contains('pbs') ? 'pbs' :
      'local'`), so any other executor was recorded as `local`. It is now
      resolved from the session config as Nextflow resolves it:
      `process.executor`, then `executor.name`, then `NXF_EXECUTOR`, then
      `local`. Checked with `-process.executor=slurm`, a `-c` config setting
      `executor.name`, `-profile fiocruz_pbs` and `NXF_EXECUTOR`, each against
      the old guess. `tests/main.metadata.nf.test` pins it with
      `-process.executor=slurm`, and fails with the old guess. A `withName` or
      `withLabel` selector can still move a single process to another
      executor; no ViralFlow configuration does.
- [x] **The wrapper could not set any NANOPORE parameter.** Neither the
      `parse_params` allow-list nor `viralflow run` knew `clair3_model`,
      `np_min_depth`, `af_threshold`, `clair3_qual`, `clair3_chunk_size`,
      `base_container` or the per-tool cpus/memory, and a params file naming
      one was rejected. The model mattered most: the default
      `r941_prom_sup_g5014` is for R9.4.1 flowcells, so an R10.4.1 user of the
      wrapper had the wrong model with no way out. Fixed 2026-09-24:
      `wrapper.NANOPORE_PARAMS` names all thirteen, the params file accepts
      them, and `viralflow run` has one option each (`--clair3-model`,
      `--np-min-depth`, `--clair3-memory`, …), built from one table in
      `wrapper/cli.py`. Unlike the older options (section 5), they default to
      `None` and are forwarded only when given, so `nextflow.config` stays the
      one source of their defaults; `--help` still shows them, and a test fails
      if those shown values stop matching `nextflow.config`. Using one outside
      NANOPORE mode is an error, in a params file or on the command line, and
      so is combining one with `--params-file`, rather than being silently
      dropped like the other CLI options there. Values are checked before
      Nextflow starts (ranges, and memory in any form the directive accepts,
      `8.GB` or `8GB`), and a container given as a local file is made absolute
      while a registry reference passes through untouched. A CLI `--x` on a
      plain `cpus = params.x` directive was checked first and does take
      effect. Covered by `tests/test_wrapper_nanopore.py`, and run for real
      through `viralflow run` under `-profile docker` with
      `--clair3-model r1041_e82_400bps_sup_v500 --np-min-depth 30
      --af-threshold 0.6 --clair3-cpus 2 --clair3-memory 3.GB`: Clair3 ran
      with that model and `--threads=2`, its container got 3 GB, and the
      summary and `run_manifest.json` record the thresholds (`masked_bases`
      588 → 738 on the truth fixture).
- [x] **The Docker image's smoke test could not fail the build.** Confirmed
      with deliberately broken builds on 2026-09-24, then fixed. The last
      `RUN` of `nanopore_base.Dockerfile` had two independent ways of hiding a
      failure, and this item originally named only the first:
      - `… && bam help > /dev/null 2>&1 || true` — the `|| true` binds to the
        whole `&&` chain, so nothing in it could fail. With minimap2 deleted,
        the image built.
      - `samtools --version | head -n 1` and the same for bcftools — Docker
        runs `RUN` under `/bin/sh` without `pipefail`, so the step saw
        `head`'s exit status. With `libhts.so*` deleted, recreating the
        historical `libhts.so.3: cannot open shared object file` failure,
        the image built **even with the `|| true` removed**. Dropping only the
        `|| true`, as first proposed here, would not have fixed it.
      `bam help` never needed the exemption either: it exits 0. The step now
      has no pipes and no fallback, and also runs `bamdash --help`, which the
      `.sing` `%test` already checked. Re-run against the fix, all three
      broken builds (libhts removed, minimap2 removed, bamdash removed) fail
      and name the tool, and the unmodified recipe still builds.
      `tests/test_container_recipes.py` now asserts the step contains no `|`
      at all, which rules out both forms, and that it runs every tool the
      pipeline takes from the image. The `.sing` has no pipes in `%test`;
      whether its shell stops at the first failure is on the section 1 SIF
      item.

### Smaller

- [ ] **`np_min_depth` and ILLUMINA's `--depth` mean different things.** Depth
      <= 20 is masked, so nanopore needs 21x; ILLUMINA's `--depth 25` (ivar
      `-m`) needs 25x. Documented, but two similarly named thresholds with
      opposite edge semantics invite mistakes.
- [x] **The coverage plot used the wrong threshold for nanopore** —
      `params.depth` (25), while masking uses `np_min_depth` (20). Fixed
      2026-09-24. bamdash draws no threshold line: `-c` feeds only the plot's
      stats, the title's "% recovery >= Nx" figure. So every NANOPORE plot
      reported recovery at ILLUMINA's threshold, a number that described no
      step of the run. On the truth fixture it read 97.77% at 25x; at
      `np_min_depth` it reads 98.02%, against the consensus's
      `callable_percent` of 98.03%. The remaining gap is expected: bamdash
      counts only bases at quality 15 or more and no deletions, while
      `samtools depth -J` counts both. `coveragePlot` now takes the threshold
      as an input instead of reading `params.depth`, `GENPLOTS` takes it from
      its caller, and `main.nf` passes `params.depth` for ILLUMINA and
      `params.np_min_depth` for NANOPORE. No offset is needed: bamdash's label
      says `>=`, but its code counts a position as recovered only when coverage
      is strictly above the threshold, the exact complement of the
      `depth <= np_min_depth` mask. Pinned by two genplots fixture tests
      (threshold 0 gives 50%; threshold 1 over depth-1 reads gives 0%, as the
      mask would), and `main-nanopore.nf.test` asserts the published plot's
      threshold is 20, not 25, and its figure is within 0.1 of
      `callable_percent`. That assertion failed before the fix, reading 25x
      and 97.77. ILLUMINA's plot has the same class of mismatch, off by one;
      see section 4.
- [x] **Silent failures in the GENPLOTS steps.** Confirmed and fixed
      2026-09-24. Two separate holes:
      - **`getMappedReads` / `getUnmappedReads` had no `pipefail`.** This item
        first said a failed `samtools sort` would publish an empty FASTQ; it
        does not — a sort that fails outright writes nothing, and
        `samtools fastq` then fails on the empty stream too. The real case is
        a sort that dies *partway through* its output, as an OOM kill does.
        Reproduced by cutting a name-sorted BAM at a BGZF block boundary (BGZF
        is written a whole block per `write()`, so that is where a killed
        process stops): under Nextflow's default `bash -ue` the task exited 0
        and published 1720 of 3382 reads, the only trace an "EOF marker is
        absent" line in `.command.err`. With `set -euo pipefail` it fails
        with the sort's status. Both modules now set it. The two ILLUMINA
        modules with the same gap are in section 4.
      - **`coveragePlot` ignored bamdash entirely**, running it through
        `subprocess.run(..., shell=True)` without looking at the exit status
        or the output — how ILLUMINA's PNG and SVG went missing unnoticed. The
        inline Python moved to `vfnext/bin/coverage_plot.py`: a missing HTML
        plot fails the task; a missing PNG or SVG is written to
        `coveragePlot_result.txt`, which GENPLOTS already logs as a warning
        for the sample. Static formats are deliberately not fatal: ILLUMINA's
        `generate_plots` image cannot draw them at all (section 4), so failing
        on them would fail every ILLUMINA run. Make them fatal once that image
        is rebuilt. A format counts as drawn only if its file exists, not on
        bamdash's exit status alone. bamdash now gets an argument list, not a
        shell line: the reference name comes from the BAM header, and SAM
        allows `;|&$` in it. The output glob was narrowed from
        `*coveragePlot*`, which also published the result file as a plot.
        The script is held to Python 3.8 in `ruff.toml`, since the
        `generate_plots` image's Python is recorded nowhere. Covered by
        `tests/test_coverage_plot.py` (six cases, fake bamdash), by the
        genplots fixture and `main-nanopore.nf.test` with the real one, and
        checked in the nanopore image with kaleido's Chromium removed: exit 0,
        HTML published, and "not produced: PNG (bamdash exited 1); SVG
        (bamdash exited 1). The HTML plot is complete."
- [x] **`runPorechop` published its trimmed reads uncompressed.** Fixed
      2026-09-24. The cost was far more than the "roughly doubling" this item
      first said: on the truth fixture the published `.chopped.fastq` was
      6.9 MB, against 188 KB compressed, about 35 times (the gzipped input
      is 87 KB). The task now
      runs `bgzip -@ ${task.cpus}` on Porechop's output and publishes only
      `${id}.chopped.fastq.gz`, which minimap2 reads as is. Streaming
      Porechop's stdout into bgzip was tried and rejected: `-abi` runs a
      helper that prints 43 lines of progress to stdout ahead of the reads,
      which shifts every FASTQ record. Porechop's own `.gz` output was not
      used either: it writes the same uncompressed temporary file, then
      compresses it with single-threaded `gzip` through a shell. The read
      content is unchanged: the decompressed output matches the old file
      byte for byte on the fixture. `main-nanopore.nf.test` asserts the `.gz`
      is published with every input read, and that no uncompressed copy is.
- [x] **`runNanoporeSummary` called `${projectDir}/bin/nanopore_summary.py`.**
      Fixed 2026-09-24. Nextflow puts `bin/` on every task's `PATH`
      (`.command.run` exports it), uploading it first on executors that do
      not share the launch host's filesystem; `${projectDir}` there still
      names the launch host's copy, which the task cannot see. It worked
      locally only because the path is the same inside and outside the
      container. Not reproduced on a cloud executor, since none is set up. The
      module now calls the script by name, and the script is now executable
      (it was `100644`; its shebang was already `#!/usr/bin/env python3`).
      `tests/test_bin_scripts.py` fails if a module outside a two-entry
      ILLUMINA allowlist uses `projectDir/bin`, or if a script a module calls
      by name loses its executable bit or shebang. The allowlisted modules
      are in section 4.
- [x] **`concat-fastq` defaulted to `--max-len 500`.** Fixed 2026-09-24.
      Reproduced with the envs' pinned seqkit 2.12 on three barcodes of 200
      reads: ~400 bp amplicons kept 200, ~1200 bp amplicons (Midnight) kept
      **0**, and a whole-genome run of 300 bp to 20 kb kept **3**, and the
      command reported success for each. `--max-len` now defaults to no
      maximum and is passed to seqkit only when given; with the fix all
      three keep 200. `--min-len` stays at 200. A `--max-len` below
      `--min-len`, which would drop every read, is rejected before anything
      is written (exit 1). Covered in `tests/test_wrapper_helpers.py`: the
      seqkit arguments with and without a maximum, the rejection, and the
      command's default. What the filters drop is still not reported; see
      section 5.
- [x] **Legacy `--inDir` discovery failed a Nanopore file named
      `*_R1.fastq` as an "orphan Illumina mate".** Fixed 2026-09-24.
      Reproduced with a new `input-fixture.nf.test` case: `--mode NANOPORE`
      over a folder holding only `orphan_R1.fastq` stopped with "Legacy input
      sample 'orphan' has an orphan Illumina mate". `_R1`/`_R2` marks a mate
      only in ILLUMINA, so in NANOPORE mode a lone mate is now an ordinary
      single-end sample, named after the whole file (`orphan_R1`) as any
      other single file is. A complete `_R1`/`_R2` pair in NANOPORE mode is
      still refused ("fastq_2 must be empty in NANOPORE mode"), since it is
      most likely Illumina data given the wrong `--mode`; a second new test
      pins that. ILLUMINA is unchanged. Two `_R1` chunk files of one name in
      NANOPORE mode (`x_R1_001`, `x_R1_002`) still stop with the existing
      "duplicate sample ID … use --samplesheet" error, which says what to do.

### Test gaps behind these

- ~~Nothing runs `main.nf --mode NANOPORE`~~. Closed by
  `integration_tests/main-nanopore.nf.test`, which found the report-path bug
  as well as the METADATA, GENPLOTS and manifest items above. It passes under
  both engines and runs in CI's integration step.
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
      `mv: cannot stat 'NC_045512.2_plot.png'`). The task still succeeds, now
      by design rather than by accident: since the "Silent failures in the
      GENPLOTS steps" fix in section 3, each ILLUMINA sample logs a warning
      naming the missing PNG and SVG instead of saying nothing. Once the image
      is fixed, make those two formats fatal in `vfnext/bin/coverage_plot.py`
      (`STATIC_FORMATS` is the list reported rather than failed). The nanopore base
      image pins `kaleido==0.2.1` with `plotly==5.24.1` and produces all three
      formats (checked on the truth fixture), so rebuilding `generate_plots`
      with the same pins is the likely fix. Installing Chrome is the
      alternative. Inherited: the image predates this branch.
- [ ] **Two ILLUMINA scripts pipe without `pipefail`, one of them into the
      consensus.** The GENPLOTS read-extraction fix in section 3 showed what
      that hides: under Nextflow's default `bash -ue`, an upstream command
      that dies partway through its output leaves the downstream one exiting
      0 on the part it got, and the task succeeds. Reproduced there on a
      name-sorted BAM cut at a block boundary: 1720 of 3382 reads published,
      no error. Still exposed:
      - `runIvar.nf`: `samtools mpileup … | ivar consensus` (and the
        `ivar variants` calls) — a pileup killed mid-stream would publish a
        truncated consensus.
      - `prepareDatabase.nf`: `esearch … | efetch -format fasta > ref.fa` — a
        failed search can leave a partial or empty reference.

      Either add `set -euo pipefail` to each script, as the NANOPORE modules
      do, or set it once for every process with
      `process.shell = ['/bin/bash', '-euo', 'pipefail']` in `nextflow.config`
      (the nf-core default). The global option is the better end state, but
      it changes how every ILLUMINA script fails, so it wants a full ILLUMINA
      run before merging, and a look for commands that exit non-zero by
      design in a pipe (a `head` closing it early gives SIGPIPE, exit 141).
- [ ] **Two ILLUMINA modules still call scripts through `$projectDir/bin`.**
      Same problem as NANOPORE's `runNanoporeSummary` (section 3): the path
      names the launch host's copy, which a cloud executor's task never sees.
      - `runIntraHostScript.nf`: `python $projectDir/bin/intrahost.py`
      - `runIvar.nf`: `python $projectDir/bin/tsv_to_vcf.py`

      Neither script is executable (`100644`), so calling them by name needs
      `chmod +x` as well. The shebangs also need checking against the
      images: `intrahost.py` has `#!/usr/bin/python3`, and the call uses
      `python`, so check whether `intrahost_analysis:1.1.0.sif` has a
      `/usr/bin/python3` and whether it is the 3.8.0 with pandas 1.5.3 that
      section 1 confirmed; `#!/usr/bin/env python3` is the safer choice.
      `tsv_to_vcf.py` already uses `env python3`, but the ivar image must
      then provide `python3`, not only `python`. Once a module is fixed,
      remove it from `PROJECT_DIR_BIN_ALLOWED` in
      `tests/test_bin_scripts.py`; the test fails until you do.
- [ ] **ILLUMINA's coverage plot counts recovery one read short of its
      consensus.** Follow-up to the section 3 NANOPORE plot-threshold fix.
      bamdash counts a position as recovered only when coverage is strictly
      above `-c`, and ILLUMINA passes `params.depth`. But `depth` is ivar's
      minimum to call a base — `docs/parameters.md`: "Positions with lower
      coverage depth will not be called" — so a position at exactly `depth`
      is called in the consensus but not recovered in the plot. Passing
      `params.depth - 1` in `main.nf`'s ILLUMINA `GENPLOTS` call would line the
      two up, once ivar's `-m` rule has been confirmed on a real run (not done:
      nothing here runs ILLUMINA end to end under Docker). Small in practice,
      but the plot is the figure people quote. The quality filters differ
      too: bamdash counts bases at quality 15 or more, ivar at `base_quality`
      (30).
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
      its own gate, the Docker METADATA failures (the specs ignored the
      engine), ILLUMINA's `toolSpecs()` still hardcoding `.sif`
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
      |:---|:---|:---|:---|
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

## 5. Wrapper follow-up branch

From a review of `wrapper/` (`cli.py`, `__init__.py`) and the container
scripts it drives, on 2026-09-24, while adding the NANOPORE options (section
3). Kept for a branch of its own: most of these change what `viralflow run`
sends to Nextflow, so each wants a test that pins the command line it builds.
Items confirmed on a real run through the wrapper say so.

### Changes what a run does

- [ ] **The wrapper switches off snpEff, mapped-read output and deduplication
      on every run, and no option can switch them back on.** Confirmed.
      `--run-snpeff`, `--write-mapped-reads` and `--dedup/--no-dedup` are
      flags that default to `False`, and every option is forwarded, so every
      run sends `--runSnpEff false --writeMappedReads false --dedup false`,
      overriding `nextflow.config`'s `true` for the first two. A NANOPORE run
      through the wrapper publishes no mapped-read FASTQs for that reason.
      Turning a flag on does not help either: it sends the String `"true"`,
      which the `params.X == true` gates reject (the boolean-parameter item in
      section 4), and a params file's `runSnpEff true` arrives the same way. So
      the documented quick start, `viralflow run --params-file
      test_files/sars-cov-2.params`, which asks for snpEff, has never run it
      through the wrapper, while the metadata layer, reading the same String as
      truthy, records snpEff as used. Needs both halves: `normalizeFlag` in
      `param_helpers.nf` (section 4) and flags here that default to `None` and
      are forwarded only when given.
- [ ] **Every other option is always forwarded at the wrapper's own copy of
      its default.** `run`'s docstring says the defaults come from
      `nextflow.config`, but thirteen options carry their own copies
      (`--min-len 75`, `--depth 25`, `--out-dir ./output/`, …) and all of them
      are sent on every run, so `nextflow.config`'s values never apply
      through the wrapper and the two can drift. It also puts ILLUMINA
      settings into NANOPORE runs: `--virus sars-cov2` is sent, and
      `run_manifest.json` records `virus: sars-cov2` for a NANOPORE run
      (confirmed). The NANOPORE options now show the pattern to follow:
      default `None`, forward only what was given, show the default in
      `--help` from a table a test checks against `nextflow.config`.
- [ ] **With `--params-file`, the other options are silently ignored.**
      `viralflow run --params-file p.txt --trim-len 5` runs with the file's
      `trimLen`. Only `--mode` (and now the NANOPORE options) raise an error.
      The wrapper cannot tell a given option from a default today — the item
      above fixes that, or click's `ctx.get_parameter_source()` can — after
      which the choice is to reject the mix or let the CLI override the file.
- [ ] **Relative paths in a params file resolve against the current directory,
      not the file.** `parse_params` calls `os.path.abspath(value)`, while
      sample-sheet paths resolve against the CSV. The shipped
      `test_files/*.params` therefore work only from the repository root, and
      `docs/quickstart.md` tells users to always write absolute paths, which
      the wrapper already makes unnecessary.
- [ ] **No way past the fixed Nextflow command.** `-resume` is always added,
      `work/` lands in whatever directory the wrapper is started from, and
      there is no way to pass `-c extra.config`, `-work-dir`, `-with-tower` or
      any other Nextflow option. A `--` passthrough would cover them.
- [ ] **CLI-given numbers are recorded as strings.** Everything the wrapper
      sends arrives in Nextflow as a String, so `run_manifest.json` records
      `mapping_quality: "30"` and `af_threshold: "0.6"` beside
      `clair3_qual: 10` from `nextflow.config` (confirmed). Harmless to the
      run, awkward for anyone reading the provenance record by type. The same
      root as the boolean item; normalizing typed parameters once on the
      Nextflow side fixes both.

### Setup and containers

- [ ] **`build-containers` does nothing for NANOPORE.** It pulls the ILLUMINA
      images from the Sylabs library and builds the pangolin and snpEff
      sandboxes. The nanopore base image is never built, so a NANOPORE user
      still follows `NANOPORE.md` by hand. A `--mode NANOPORE` (or `--all`)
      that runs `singularity build` on `Nanopore_baseContainer.sing`, or
      `docker build` on the Dockerfile for `-profile docker`, would close it.
- [ ] **Singularity only.** `pull_containers.py`, `build_containers.py` and
      `update-pangolin*` call `singularity` by name, so an Apptainer-only host
      fails, although `NANOPORE.md` documents Apptainer's `library://` setup
      and `-profile apptainer` exists. `--profile`'s help lists neither
      `docker` nor `singularity`.
- [ ] **`build_containers.py` insists on `/usr/local/bin/unsquashfs`**, and its
      advice hard-codes `$HOME/miniconda3/envs/viralflow/bin/unsquashfs`, a
      layout nothing else in the repository sets up.
- [ ] **The wrapper works only as an editable install from a clone.**
      `VF_ROOT_PATH` is two directories up from `wrapper/cli.py`, and
      `setup.py` packages `wrapper` alone. After `pip install .` or from a
      wheel, `viralflow run` points Nextflow at a `vfnext/main.nf` that does
      not exist. Either ship `vfnext/` as package data or detect the missing
      tree and say so.
- [ ] **A third copy of the Nextflow version.** `run_vfnext` defaults
      `NXF_VER` to `26.04.6`, beside the CI workflow's `NXF_VER` and the docs;
      `nextflow.config` asks only for `>=26.04.0`. Nothing keeps the three in
      step.

### `concat-fastq`

- [ ] **Its output needs the deprecated input path.** It writes
      `filtered/<barcode>.concat.fastq.gz` and no sample sheet, so the only way
      to run the result is `--inDir`, deprecated for v3, which names the
      samples `barcode01.concat` and so on. Writing a `samplesheet.csv` beside
      the FASTQs would let it feed `--samplesheet` directly.
- [ ] Smaller: it writes into the input tree; it needs `seqkit` on the host,
      which nothing declares or checks up front, so a missing `seqkit` shows
      up as one failure per barcode; the command has no docstring, so its
      `--help` is empty; and it never says how many reads its length filters
      dropped, so a `--min-len` or `--max-len` that does not suit the
      protocol still loses reads without a word. Printing "kept N of M
      reads" per barcode would fix that; counting the input costs a second
      pass over it, or a counting stage in the pipe. The `--max-len 500`
      default that made this matter most is fixed (section 3).

### Documentation

- [ ] **The quick start uses command names that do not exist**, in all three
      languages: `viralflow add_entry_to_snpeff --org_name … --genome_code …`,
      `update_pangolin` and `update_pangolin_data`. The commands are
      `add-entry-to-snpeff --org-name … --genome-code …`, `update-pangolin`
      and `update-pangolin-data`; the underscore forms fail with "No such
      command" (confirmed).
- [ ] **No NANOPORE quick start.** The docs cover only ILLUMINA, and
      `test_files/` has no NANOPORE params file (`test_files/nanopore/` holds
      just a reference). The truth fixture and
      `vfnext/tests/data/input/nanopore-truth.csv` would make a working
      example.
- [ ] **The CLI option names are documented nowhere.** The parameter tables
      use the params-file names (`np_min_depth`); only `--help` shows
      `--np-min-depth`.
- [ ] **`docs/parameters.md` has gone stale in two NANOPORE rows**, in all
      three languages: `trimLen` says NANOPORE "masks the bases … with
      bamUtil", but it has soft-clipped since `--clip`; and `clair3_qual` is
      "Minimum variant quality for Clair3 to report a call", when Clair3
      reports every call and only labels those at or below it `LowQual` —
      which the consensus then ignores (the section 3 `LowQual` item).

### Code health

- [ ] `parse_csv` in `wrapper/__init__.py` is dead code.
- [ ] `pull_containers.py` and `build_containers.py` are scripts run at import
      through `sys.executable`, with module-level state (`success`,
      `failed_containers`), so nothing can test them without running
      Singularity. Moving the logic into functions the wrapper calls would
      let the existing mocked-subprocess tests cover them.

---

## 6. Repository maintenance, after the PR merges

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
