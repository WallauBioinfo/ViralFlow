# TODO

Work that is known about but not done. Kept in the repository so it survives a
handover: anyone picking up `feature/add_nanopore_support` or starting the
ILLUMINA follow-up should find the open threads here rather than in a chat log.

Sections are ordered by where the work belongs, not by priority. Within each
section, items that could affect analysis results come first.

---

## 1. Needs a Linux machine with Singularity

Cannot be settled on macOS/Docker.

- [ ] **Verify `bcftools` works inside the built SIF.**
      ```bash
      singularity exec vfnext/containers/baseContainer.sif bcftools --version
      ```
      `bcftools` is built from the GitHub *source archive*, which — unlike the
      release tarball `samtools` uses — does not bundle htslib, so it links
      dynamically against `/usr/local/lib/libhts.so`. Nothing in
      `Nanopore_baseContainer.sing` runs `ldconfig`. Building the Docker
      equivalent surfaced this as a hard runtime failure
      (`libhts.so.3: cannot open shared object file`) and
      `nanopore_base.Dockerfile` fixes it with an `ldconfig` after the htslib
      install. Whether the SIF is affected depends on incidental linker-cache
      state, so it needs checking rather than assuming. If it fails, apply the
      same one-line fix to the `.sing`.

- [ ] **Confirm the truth test passes under Singularity**, not only under
      `-profile docker`. The container directive for Clair3 strips the
      `docker://` prefix only when the engine is Docker; the Singularity path is
      unchanged but has not been run since that change.

- [ ] **Confirm Singularity/Apptainer pull and run `docker://` images cleanly**,
      including by digest. Expected to work — `clair3_container` already relies
      on it for NANOPORE, so the truth test exercises it — but it has only been
      run under Docker on this Mac, and the answer gates the item below.

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

- [ ] **Confirm `intrahost_analysis:1.1.0.sif` really carries what its recipe
      says.** `runIntraHostScript` now uses it instead of a registry image, on
      the strength of a def file deleted eight months ago (`c2be157`); the
      published library image could have moved since. One command settles it:
      ```bash
      singularity exec vfnext/containers/intrahost_analysis:1.1.0.sif \
        python -c "import sys, pandas, numpy, Bio; print(sys.version, pandas.__version__, numpy.__version__, Bio.__version__)"
      ```
      Expect Python 3.8.x, pandas 1.5.3, numpy 1.23, biopython 1.81. The
      equivalents were verified here against a Python 3.8 environment built to
      those pins — `intrahost.py` compiles and runs `--help` clean — so the only
      untested link is the image itself. If the Python differs, update
      `per-file-target-version` in `ruff.toml` and
      `INTRAHOST_CONTAINER_PYTHON` in `tests/test_container_recipes.py` to
      match; a newer Python needs no other change.

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
      3.**
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
      apart from the wrapper early on. Moved to section 3.
- [x] **An `annotations` subworkflow** — **agreed on the PR**, deferred until
      NANOPORE merges, on the grounds that this PR already carries a lot of
      non-nanopore change. Moved to section 3.
- [x] Should the `coveragePlot` import live in `GENPLOTS.nf`? — **agreed on the
      PR**, deferred to the ILLUMINA branch. Moved to section 3.
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

## 3. ILLUMINA follow-up branch

Deliberately kept out of the nanopore PR.

- [ ] **Move `main.nf`, `workflows/` and `modules/` to the repository root**,
      dropping `vfnext/`. Agreed on PR #47. Touches every `includeConfig` and
      `$projectDir` path, the wrapper's `root_path`, `.readthedocs.yaml` and the
      CI workflow, so it wants its own PR with nothing else in it.
- [ ] **Add an `annotations` subworkflow** grouping snpEff, pangolin, nextclade
      and compileOutput. Agreed on PR #47. The reviewer's point is that all
      three could serve nanopore output too, so this is what would let NANOPORE
      reuse the ILLUMINA annotation stack rather than reimplement it.
- [ ] **Move the `coveragePlot` import into `GENPLOTS.nf`.** Agreed on PR #47.
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

## 4. Repository maintenance, after the PR merges

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
