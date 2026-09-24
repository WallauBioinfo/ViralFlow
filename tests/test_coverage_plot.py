"""vfnext/bin/coverage_plot.py: every way bamdash can fail is either fatal or said.

The script replaced inline Python that ran bamdash through shell=True and never
looked at the result, so a plot that was not drawn left no trace but a missing
file. That is how ILLUMINA's PNG and SVG plots went unpublished in every run.
"""

import importlib.util
import tempfile
import unittest
from contextlib import chdir
from pathlib import Path
from types import SimpleNamespace

SCRIPT = Path(__file__).resolve().parents[1] / "vfnext" / "bin" / "coverage_plot.py"


def load_script():
    spec = importlib.util.spec_from_file_location("coverage_plot", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


coverage_plot = load_script()


class FakeBamdash:
    """Stands in for subprocess.run: records calls, writes or withholds output."""

    def __init__(self, fail=(), exit_zero_without_output=()):
        self.fail = set(fail)
        self.silent = set(exit_zero_without_output)
        self.calls = []

    def __call__(self, command, **_kwargs):
        self.calls.append(command)
        reference = command[command.index("-r") + 1]
        fmt = command[command.index("-e") + 1] if "-e" in command else "html"
        if fmt in self.fail:
            return SimpleNamespace(returncode=1)
        if fmt not in self.silent:
            Path(f"{reference}_plot.{fmt}").write_text(fmt)
        return SimpleNamespace(returncode=0)


class CoveragePlotTests(unittest.TestCase):
    def run_script(self, bamdash, mapped=(3, "NC_045512.2")):
        self.addCleanup(
            setattr, coverage_plot, "mapped_reads", coverage_plot.mapped_reads
        )
        coverage_plot.mapped_reads = lambda _bam: mapped
        with tempfile.TemporaryDirectory() as temp_dir, chdir(temp_dir):
            status = coverage_plot.main(
                [
                    "--bam",
                    "sample.bam",
                    "--sample",
                    "sample",
                    "--threshold",
                    "20",
                    "--result-file",
                    "result.txt",
                ],
                run=bamdash,
            )
            files = sorted(path.name for path in Path(".").iterdir())
            result = (
                Path("result.txt").read_text() if Path("result.txt").exists() else None
            )
        return status, files, result

    def test_all_three_plots_are_published(self):
        status, files, result = self.run_script(FakeBamdash())

        self.assertEqual(status, 0)
        self.assertEqual(
            files,
            [
                "sample_coveragePlot.html",
                "sample_coveragePlot.png",
                "sample_coveragePlot.svg",
            ],
        )
        self.assertIsNone(result)

    def test_a_failed_html_plot_fails_the_task(self):
        status, files, _ = self.run_script(FakeBamdash(fail={"html"}))

        self.assertEqual(status, 1)
        self.assertNotIn("sample_coveragePlot.html", files)

    def test_an_exit_status_of_zero_is_not_taken_on_trust(self):
        """bamdash's status alone is not enough: the file must exist too."""
        status, _, _ = self.run_script(FakeBamdash(exit_zero_without_output={"html"}))

        self.assertEqual(status, 1)

    def test_a_failed_static_plot_is_reported_not_fatal(self):
        """ILLUMINA's generate_plots image cannot export PNG or SVG today.

        Making that fatal would fail every ILLUMINA run; hiding it is how it
        went unnoticed. It is reported through the result file, which GENPLOTS
        logs as a warning for the sample.
        """
        status, files, result = self.run_script(FakeBamdash(fail={"png", "svg"}))

        self.assertEqual(status, 0)
        self.assertEqual(files, ["result.txt", "sample_coveragePlot.html"])
        self.assertIn("sample", result)
        self.assertIn("PNG", result)
        self.assertIn("SVG", result)
        self.assertIn("HTML", result)

    def test_no_mapped_reads_skips_bamdash_and_says_so(self):
        bamdash = FakeBamdash()
        status, files, result = self.run_script(bamdash, mapped=(0, "NC_045512.2"))

        self.assertEqual(status, 0)
        self.assertEqual(bamdash.calls, [])
        self.assertEqual(files, ["result.txt"])
        self.assertIn("No mapped reads", result)

    def test_bamdash_gets_an_argument_list_not_a_shell_line(self):
        """Reference names come from the BAM header, and SAM allows ;|&$ in them."""
        reference = "chr1;touch pwned"
        bamdash = FakeBamdash()
        status, files, _ = self.run_script(bamdash, mapped=(3, reference))

        self.assertEqual(status, 0)
        self.assertNotIn("pwned", files)
        for command in bamdash.calls:
            self.assertIsInstance(command, list)
            self.assertEqual(command[command.index("-r") + 1], reference)
            self.assertEqual(command[command.index("-c") + 1], "20")


if __name__ == "__main__":
    unittest.main()
