import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import call, patch

from click.testing import CliRunner

from wrapper import (
    add_entries_to_DB,
    build_containers,
    concat_fastqs,
    update_pangolin,
    update_pangolin_data,
)
from wrapper.cli import cli


FASTQ = b"@read1\nACGT\n+\nIIII\n"


class HelperCommandTests(unittest.TestCase):
    def test_snpeff_arguments_are_literal_and_checked(self):
        with patch("wrapper.subprocess.run") as run:
            add_entries_to_DB(
                "/viral flow", "Dengue; touch /tmp/pwned", "NC_001474.2", "amd64"
            )

        run.assert_called_once_with(
            [
                "bash",
                "/viral flow/vfnext/containers/add_entries_SnpeffDB.sh",
                "Dengue; touch /tmp/pwned",
                "NC_001474.2",
                "amd64",
            ],
            cwd=Path("/viral flow/vfnext/containers"),
            check=True,
        )

    def test_snpeff_rejects_control_characters_and_option_like_codes(self):
        with self.assertRaisesRegex(ValueError, "control characters"):
            add_entries_to_DB("/viralflow", "Dengue\nInjected", "NC_001474.2", "amd64")
        with self.assertRaisesRegex(ValueError, "must start"):
            add_entries_to_DB("/viralflow", "Dengue", "--help", "amd64")

    def test_container_build_steps_stop_and_propagate_failure(self):
        failure = subprocess.CalledProcessError(
            9, [sys.executable, "pull_containers.py", "amd64"]
        )
        with patch("wrapper.subprocess.run", side_effect=failure) as run:
            with self.assertRaises(subprocess.CalledProcessError):
                build_containers("/viralflow", "amd64")
        self.assertEqual(run.call_count, 1)

    def test_container_build_and_updates_use_checked_argument_lists(self):
        containers = Path("/viralflow/vfnext/containers")
        with patch("wrapper.subprocess.run") as run:
            build_containers("/viralflow", "arm64")
            update_pangolin("/viralflow")
            update_pangolin_data("/viralflow")

        self.assertEqual(
            run.call_args_list,
            [
                call(
                    [sys.executable, "pull_containers.py", "arm64"],
                    cwd=containers,
                    check=True,
                ),
                call(
                    [sys.executable, "build_containers.py", "arm64"],
                    cwd=containers,
                    check=True,
                ),
                call(
                    [
                        "singularity",
                        "exec",
                        "--writable",
                        "./pangolin:4.4.sif",
                        "pangolin",
                        "--update",
                    ],
                    cwd=containers,
                    check=True,
                ),
                call(
                    [
                        "singularity",
                        "exec",
                        "--writable",
                        "./pangolin:4.4.sif",
                        "pangolin",
                        "--update-data",
                    ],
                    cwd=containers,
                    check=True,
                ),
            ],
        )

    def test_click_reports_helper_failure(self):
        runner = CliRunner()
        failure = subprocess.CalledProcessError(7, ["singularity", "exec"])
        with patch("wrapper.cli._update_pangolin", side_effect=failure):
            result = runner.invoke(cli, ["update-pangolin"])
        self.assertEqual(result.exit_code, 1)
        self.assertIn("Command failed with exit status 7", result.output)


class ConcatFastqTests(unittest.TestCase):
    def make_seqkit(self, root):
        bin_dir = root / "bin"
        bin_dir.mkdir()
        seqkit = bin_dir / "seqkit"
        seqkit.write_text("#!/bin/sh\ncat\n", encoding="utf-8")
        seqkit.chmod(0o755)
        return bin_dir

    def test_literal_names_and_atomic_success_for_uncompressed_input(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            bin_dir = self.make_seqkit(root)
            barcode = root / "barcode;touch injected"
            barcode.mkdir()
            (barcode / "reads $one.fastq").write_bytes(FASTQ)
            environment = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
            with patch.dict(os.environ, environment, clear=True):
                concat_fastqs(root, "barcode;", ".fastq", 1, 100)

            output = root / "filtered" / "barcode;touch injected.concat.fastq.gz"
            with gzip.open(output, "rb") as result:
                self.assertEqual(result.read(), FASTQ)
            self.assertFalse((root / "injected").exists())
            self.assertEqual(list((root / "filtered").glob("*.tmp")), [])

    def test_compressed_batch_continues_then_reports_all_failures(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            bin_dir = self.make_seqkit(root)
            good = root / "barcode01"
            bad = root / "barcode02"
            good.mkdir()
            bad.mkdir()
            with gzip.open(good / "reads.fastq.gz", "wb") as output:
                output.write(FASTQ)
            (bad / "reads.fastq.gz").write_bytes(b"not gzip")
            environment = {**os.environ, "PATH": f"{bin_dir}:{os.environ['PATH']}"}
            with patch.dict(os.environ, environment, clear=True):
                with self.assertRaisesRegex(
                    RuntimeError, "barcode02: failed pipeline stages: reader"
                ):
                    concat_fastqs(root, "barcode", ".fastq.gz", 1, 100)

            self.assertTrue((root / "filtered" / "barcode01.concat.fastq.gz").is_file())
            self.assertFalse((root / "filtered" / "barcode02.concat.fastq.gz").exists())
            self.assertEqual(list((root / "filtered").glob("*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
