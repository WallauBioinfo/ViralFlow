import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from wrapper import run_vfnext
from wrapper.cli import cli


class WrapperModeTests(unittest.TestCase):
    def test_parameter_file_nanopore_mode_is_not_overridden(self):
        with patch("wrapper.subprocess.run") as run:
            with patch("wrapper.parse_params", return_value="--mode NANOPORE -resume"):
                run_vfnext("/viralflow/", "params.txt", None)

        command = run.call_args.args[0]
        self.assertEqual(command.count("--mode"), 1)
        self.assertEqual(command[command.index("--mode") + 1], "NANOPORE")

    def test_parameter_file_without_mode_adds_no_mode(self):
        with patch("wrapper.subprocess.run") as run:
            with patch("wrapper.parse_params", return_value="--virus custom -resume"):
                run_vfnext("/viralflow/", "params.txt", None)

        self.assertNotIn("--mode", run.call_args.args[0])

    def test_parameter_file_rejects_explicit_cli_mode(self):
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as temp_dir:
            params_file = Path(temp_dir) / "params.txt"
            params_file.write_text("mode NANOPORE\n", encoding="utf-8")
            with patch("wrapper.cli._run_vfnext") as run:
                result = runner.invoke(
                    cli,
                    ["run", "--params-file", str(params_file), "--mode", "ILLUMINA"],
                )

        self.assertEqual(result.exit_code, 2)
        self.assertIn("--mode cannot be used with --params-file", result.output)
        run.assert_not_called()

    def test_cli_only_run_defaults_to_illumina(self):
        with patch("wrapper.subprocess.run") as run:
            run_vfnext("/viralflow/", None, None, {"virus": "custom"})

        command = run.call_args.args[0]
        self.assertEqual(command[command.index("--mode") + 1], "ILLUMINA")

    def test_cli_only_run_preserves_explicit_nanopore_mode(self):
        with patch("wrapper.subprocess.run") as run:
            run_vfnext("/viralflow/", None, "NANOPORE", {"virus": "custom"})

        command = run.call_args.args[0]
        self.assertEqual(command[command.index("--mode") + 1], "NANOPORE")


if __name__ == "__main__":
    unittest.main()
