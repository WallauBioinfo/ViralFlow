import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from wrapper import parse_params, run_vfnext
from wrapper.cli import cli


class WrapperModeTests(unittest.TestCase):
    def test_parameter_file_nanopore_mode_is_not_overridden(self):
        with patch("wrapper.subprocess.run") as run:
            with patch("wrapper.parse_params", return_value=["--mode", "NANOPORE", "-resume"]):
                run_vfnext("/viralflow/", "params.txt", None)

        command = run.call_args.args[0]
        self.assertEqual(command.count("--mode"), 1)
        self.assertEqual(command[command.index("--mode") + 1], "NANOPORE")

    def test_parameter_file_without_mode_adds_no_mode(self):
        with patch("wrapper.subprocess.run") as run:
            with patch("wrapper.parse_params", return_value=["--virus", "custom", "-resume"]):
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

    def test_cli_paths_with_spaces_remain_one_argument(self):
        cases = {
            "samplesheet": "--samplesheet",
            "inDir": "--inDir",
            "outDir": "--outDir",
            "referenceGFF": "--referenceGFF",
            "referenceGenome": "--referenceGenome",
            "primersBED": "--primersBED",
        }
        for key, option in cases.items():
            with self.subTest(key=key), patch("wrapper.subprocess.run") as run:
                path = f"/data/run 1/{key}"
                run_vfnext("/viralflow/", None, None, {key: path})
                command = run.call_args.args[0]
                self.assertEqual(command[command.index(option) + 1], path)

    def test_parameter_file_path_with_spaces_remains_one_argument(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            params_file = Path(temp_dir) / "params.txt"
            relative_path = "data/run 1/samples.csv"
            params_file.write_text(f"samplesheet {relative_path}\n", encoding="utf-8")
            args = parse_params(params_file)

        self.assertEqual(args[args.index("--samplesheet") + 1], str(Path(relative_path).absolute()))

    def test_empty_and_null_parameter_values_are_omitted(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            params_file = Path(temp_dir) / "params.txt"
            params_file.write_text("primersBED\nreferenceGFF null\nmode NANOPORE\n", encoding="utf-8")
            args = parse_params(params_file)

        self.assertEqual(args, ["--mode", "NANOPORE", "-resume"])

    def test_non_path_parameter_rejects_multiple_values(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            params_file = Path(temp_dir) / "params.txt"
            params_file.write_text("virus custom extra\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "virus accepts a single value"):
                parse_params(params_file)

    def test_profile_and_resume_are_added_once(self):
        with patch("wrapper.subprocess.run") as run:
            run_vfnext("/viralflow/", None, None, {"virus": "custom"}, "apptainer")

        command = run.call_args.args[0]
        self.assertEqual(command.count("-resume"), 1)
        self.assertEqual(command.count("-profile"), 1)
        self.assertEqual(command[command.index("-profile") + 1], "apptainer")

    def test_printed_command_quotes_paths_with_spaces(self):
        path = "/data/run 1/samples.csv"
        with patch("wrapper.subprocess.run"), patch("builtins.print") as output:
            run_vfnext("/viralflow/", None, None, {"samplesheet": path})

        self.assertIn("'/data/run 1/samples.csv'", output.call_args.args[0])


if __name__ == "__main__":
    unittest.main()
