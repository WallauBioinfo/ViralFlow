# Nanopore truth fixture

This fixture is generated deterministically from the repository SARS-CoV-2
reference with five designed variants and a 500-base read gap at reference
positions 15001-15500. Minimap2 may soft-clip reads at the gap boundaries, so
the exact masked interval is asserted through the committed completeness
metrics rather than assumed to equal the raw read gap.

The reads are high-quality, ONT-length synthetic reads. They validate workflow
wiring and exact truth recovery; they are not intended to model empirical ONT
error profiles.

Regenerate the committed fixture from the `vfnext` directory with:

```bash
python3 tests/integration/data/nanopore_truth/generate_fixture.py
```

That script rewrites `SHA256SUMS` for the five data files it produces. It does
not touch `expected_containers.tsv`, which is maintained by hand.

## Container pinning

`expected_containers.tsv` records the images the truth values were generated
with, but the two are pinned differently:

- **clair3** is pulled, so it is pinned by digest. The same digest is declared
  in the `params` block of `integration_tests/nanopore-truth.nf.test`, and the
  test asserts the two agree. Moving to a new Clair3 image means regenerating
  the truth values and updating both places.
- **nanopore_base** is built locally from `containers/Nanopore_baseContainer.sing`
  and its `sha256` is recorded as `NA`. That build is not reproducible byte for
  byte: `%post` runs an unpinned `apt-get install` and clones Porechop_ABI and
  bamUtil from their default branches, and SIF images embed build metadata. Two
  people building from the same definition get different checksums, so the test
  only asserts the file is present.

The checksum of the base container that actually ran is captured at runtime in
`RUN_METADATA/container_manifest.tsv`, which is where that provenance belongs.
