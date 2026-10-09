"""The wrapper's NANOPORE parameters, on the command line and in a params file.

The wrapper used to know none of them: a params file naming one was rejected,
and `viralflow run` had no option for any. The one that matters most is the
Clair3 model, whose default is for R9.4.1 flowcells - an R10.4.1 user of the
wrapper had no way to change it.
"""

import re
import tempfile
from contextlib import chdir
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from wrapper import NANOPORE_PARAMS, parse_params
from wrapper.cli import NANOPORE_OPTIONS, cli

PROJECT_ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_CONFIG = PROJECT_ROOT / "vfnext" / "nextflow.config"


def write_params(directory, text):
    params_file = Path(directory) / "params.txt"
    params_file.write_text(text, encoding="utf-8")
    return params_file


def forwarded(args, option):
    return args[args.index(option) + 1]


def nextflow_config_default(name):
    """The literal a params entry in nextflow.config assigns, quotes removed."""
    match = re.search(
        rf"^\s*{name}\s*=\s*(.+?)\s*(//.*)?$",
        NEXTFLOW_CONFIG.read_text(),
        flags=re.MULTILINE,
    )
    return match.group(1).strip("\"'") if match else None


class NanoporeParamsFileTests(unittest.TestCase):
    def test_every_nanopore_parameter_is_accepted_in_nanopore_mode(self):
        values = {
            "clair3_model": "r1041_e82_400bps_sup_v500",
            "clair3_qual": "12",
            "clair3_chunk_size": "5000",
            "af_threshold": "0.7",
            "np_min_depth": "30",
            "clair3_container": "docker://hkubal/clair3@sha256:abc",
            "base_container": "viralflow/nanopore-base:2.0.0a1",
            "porechop_cpus": "2",
            "porechop_memory": "8.GB",
            "minimap_cpus": "2",
            "minimap_memory": "8GB",
            "clair3_cpus": "8",
            "clair3_memory": "16.GB",
        }
        self.assertEqual(sorted(values), sorted(NANOPORE_PARAMS))
        body = "mode NANOPORE\n" + "".join(f"{k} {v}\n" for k, v in values.items())
        with tempfile.TemporaryDirectory() as temp_dir:
            args = parse_params(write_params(temp_dir, body))

        for key, value in values.items():
            with self.subTest(key=key):
                self.assertEqual(forwarded(args, f"--{key}"), value)

    def test_nanopore_parameters_are_rejected_outside_nanopore_mode(self):
        for mode_line in ("mode ILLUMINA\n", ""):
            with self.subTest(mode=mode_line.strip() or "unset"):
                with tempfile.TemporaryDirectory() as temp_dir:
                    params_file = write_params(
                        temp_dir, f"{mode_line}clair3_model r1041_e82_400bps_sup_v500\n"
                    )
                    with self.assertRaisesRegex(
                        ValueError, "clair3_model.*only apply to mode NANOPORE"
                    ):
                        parse_params(params_file)

    def test_a_local_container_file_is_made_absolute(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            image = Path(temp_dir) / "baseContainer.sif"
            image.write_bytes(b"")
            params_file = write_params(
                temp_dir, f"mode NANOPORE\nbase_container {image.name}\n"
            )
            with chdir(temp_dir):
                args = parse_params(params_file)

        self.assertEqual(forwarded(args, "--base_container"), str(image.resolve()))

    def test_a_registry_image_reference_is_left_alone(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            args = parse_params(
                write_params(
                    temp_dir,
                    "mode NANOPORE\n"
                    "base_container viralflow/nanopore-base:2.0.0a1\n"
                    "clair3_container docker://hkubal/clair3@sha256:abc\n",
                )
            )

        self.assertEqual(
            forwarded(args, "--base_container"), "viralflow/nanopore-base:2.0.0a1"
        )
        self.assertEqual(
            forwarded(args, "--clair3_container"), "docker://hkubal/clair3@sha256:abc"
        )


class NanoporeCliTests(unittest.TestCase):
    def invoke(self, *arguments):
        with patch("wrapper.cli._run_vfnext") as run:
            result = CliRunner().invoke(cli, ["run", *arguments])
        return result, run

    def test_options_are_forwarded_under_their_nextflow_names(self):
        result, run = self.invoke(
            "--mode",
            "NANOPORE",
            "--clair3-model",
            "r1041_e82_400bps_sup_v500",
            "--np-min-depth",
            "30",
            "--af-threshold",
            "0.7",
            "--clair3-memory",
            "16.GB",
        )

        self.assertEqual(result.exit_code, 0, result.output)
        cli_params = run.call_args.args[3]
        self.assertEqual(cli_params["clair3_model"], "r1041_e82_400bps_sup_v500")
        self.assertEqual(cli_params["np_min_depth"], 30)
        self.assertEqual(cli_params["af_threshold"], 0.7)
        self.assertEqual(cli_params["clair3_memory"], "16.GB")

    def test_unset_options_are_not_forwarded(self):
        """nextflow.config stays the one source of the NANOPORE defaults."""
        result, run = self.invoke("--mode", "NANOPORE")

        self.assertEqual(result.exit_code, 0, result.output)
        cli_params = run.call_args.args[3]
        self.assertEqual([key for key in NANOPORE_PARAMS if key in cli_params], [])

    def test_options_are_rejected_outside_nanopore_mode(self):
        for mode in (["--mode", "ILLUMINA"], []):
            with self.subTest(mode=mode or "default"):
                result, run = self.invoke(*mode, "--clair3-model", "r1041_x")
                self.assertEqual(result.exit_code, 2)
                self.assertIn("only apply to --mode NANOPORE", result.output)
                run.assert_not_called()

    def test_options_are_rejected_with_a_params_file(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            params_file = write_params(temp_dir, "mode NANOPORE\n")
            result, run = self.invoke(
                "--params-file", str(params_file), "--np-min-depth", "30"
            )

        self.assertEqual(result.exit_code, 2)
        self.assertIn("set them in the parameter file", result.output)
        run.assert_not_called()

    def test_values_are_checked_before_nextflow_starts(self):
        bad = {
            "--af-threshold": "1.5",
            "--np-min-depth": "-1",
            "--clair3-cpus": "0",
            "--clair3-memory": "lots",
        }
        for option, value in bad.items():
            with self.subTest(option=option):
                result, run = self.invoke("--mode", "NANOPORE", option, value)
                self.assertEqual(result.exit_code, 2, result.output)
                run.assert_not_called()

    def test_a_local_container_file_is_made_absolute(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            image = Path(temp_dir) / "baseContainer.sif"
            image.write_bytes(b"")
            with chdir(temp_dir):
                result, run = self.invoke(
                    "--mode", "NANOPORE", "--base-container", image.name
                )

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(run.call_args.args[3]["base_container"], str(image.resolve()))


class NanoporeOptionTableTests(unittest.TestCase):
    def test_the_cli_covers_every_nanopore_parameter(self):
        self.assertEqual(
            sorted(option.param for option in NANOPORE_OPTIONS),
            sorted(NANOPORE_PARAMS),
        )

    def test_documented_defaults_match_nextflow_config(self):
        """Help shows each default; nextflow.config is what applies it."""
        for option in NANOPORE_OPTIONS:
            if option.default is None:
                continue
            with self.subTest(param=option.param):
                self.assertEqual(
                    str(option.default), nextflow_config_default(option.param)
                )

    def test_help_lists_the_options_with_their_defaults(self):
        result = CliRunner().invoke(cli, ["run", "--help"])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("--clair3-model", result.output)
        self.assertIn("r941_prom_sup_g5014", result.output)


if __name__ == "__main__":
    unittest.main()
