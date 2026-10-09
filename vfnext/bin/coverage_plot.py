#!/usr/bin/env python
"""Draw one sample's coverage plot with bamdash, as HTML, PNG and SVG.

Every way this can go wrong is either fatal or said. The inline script it
replaced ran bamdash through shell=True and never looked at the result, so a
plot that was not drawn left no trace but a missing file - which is how
ILLUMINA's PNG and SVG plots went unpublished in every run.

- HTML is required: it is the plot people open, and both images can draw it.
  If it is not produced, the task fails.
- PNG and SVG are reported, not fatal. They need kaleido 0.2's bundled
  Chromium; the nanopore base image has it, but ILLUMINA's generate_plots image
  has kaleido 1.x and no Chrome, so neither is ever drawn there. Failing on
  that would fail every ILLUMINA run until the image is rebuilt (TODO.md,
  section 4). Instead the result file names what is missing, and GENPLOTS logs
  it as a warning for the sample.
- A sample with no mapped reads gets no plot, and the result file says so.

bamdash's exit status is not taken on trust: a format counts as drawn only when
its file exists. It is also run from an argument list, never a shell line,
because the reference name comes from the BAM header and SAM allows ;|&$ in it.

Runs in ILLUMINA's generate_plots image, whose Python version is not recorded
anywhere, so ruff.toml holds this file to 3.8 syntax.
"""

import argparse
import os
import subprocess
import sys

STATIC_FORMATS = ("png", "svg")


def mapped_reads(bam):
    """Mapped read count, and the first reference the plot is drawn for."""
    import pysam

    with pysam.AlignmentFile(bam, "rb") as alignments:
        if not alignments.references:
            raise ValueError(f"No references found in the header of {bam}")
        return alignments.count(), alignments.references[0]


def draw(bam, reference, threshold, fmt, target, run):
    """Draw one format; return None on success, or what went wrong."""
    command = ["bamdash", "-b", bam, "-r", reference, "-c", str(threshold)]
    if fmt != "html":
        command += ["-e", fmt]
    status = run(command).returncode
    # bamdash names its output after the reference, in the working directory.
    produced = f"{reference}_plot.{fmt}"
    if status != 0:
        return f"bamdash exited {status}"
    if not os.path.isfile(produced):
        return f"bamdash exited 0 but wrote no {produced}"
    os.replace(produced, target)
    return None


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--threshold",
        required=True,
        help="bamdash -c: depth a position must exceed to count as recovered",
    )
    parser.add_argument(
        "--result-file",
        required=True,
        help="written only when there is something to report; GENPLOTS logs it",
    )
    return parser.parse_args(argv)


def main(argv=None, run=subprocess.run):
    args = parse_args(argv)

    count, reference = mapped_reads(args.bam)
    if count == 0:
        with open(args.result_file, "w") as handle:
            handle.write(
                f"No mapped reads were found in the sorted BAM file for sample "
                f"{args.sample}. The coverage plot will not be generated for it."
            )
        return 0

    def target(fmt):
        return f"{args.sample}_coveragePlot.{fmt}"

    problem = draw(args.bam, reference, args.threshold, "html", target("html"), run)
    if problem:
        print(
            f"Coverage plot for sample {args.sample}: the HTML plot was not "
            f"produced ({problem}).",
            file=sys.stderr,
        )
        return 1

    missing = []
    for fmt in STATIC_FORMATS:
        problem = draw(args.bam, reference, args.threshold, fmt, target(fmt), run)
        if problem:
            missing.append(f"{fmt.upper()} ({problem})")
    if missing:
        with open(args.result_file, "w") as handle:
            handle.write(
                f"Coverage plot for sample {args.sample}: not produced: "
                f"{'; '.join(missing)}. The HTML plot is complete."
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
