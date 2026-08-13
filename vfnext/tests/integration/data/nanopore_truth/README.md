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
