import importlib.util
import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch


MODULE_PATH = (
    Path(__file__).parents[1] / "vfnext" / "containers" / "spython_functions.py"
)
SPEC = importlib.util.spec_from_file_location(
    "spython_functions_under_test", MODULE_PATH
)
SPYTHON_FUNCTIONS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SPYTHON_FUNCTIONS)


class ContainerPullTests(unittest.TestCase):
    def test_pull_uses_literal_checked_arguments(self):
        container = (
            "viralflow-amd64",
            "name;touch pwned",
            "org/project/name;touch pwned",
        )
        with patch.object(SPYTHON_FUNCTIONS.subprocess, "run") as run:
            SPYTHON_FUNCTIONS.container_pull("/containers with spaces", [container])

        run.assert_called_once_with(
            [
                "singularity",
                "pull",
                "-F",
                "name;touch pwned.sif",
                "library://org/project/name;touch pwned",
            ],
            cwd="/containers with spaces",
            check=True,
        )

    def test_retry_exhaustion_raises_after_three_attempts(self):
        container = ("project", "missing:1", "org/project/missing:1")
        failure = subprocess.CalledProcessError(1, ["singularity", "pull"])
        with patch.object(
            SPYTHON_FUNCTIONS, "container_pull", side_effect=failure
        ) as pull:
            with patch.object(
                SPYTHON_FUNCTIONS, "check_containers", return_value=[container]
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "Failed to download containers after retries"
                ):
                    SPYTHON_FUNCTIONS.containers_routine_pull(
                        [container], "/containers", [container]
                    )
        self.assertEqual(pull.call_count, 3)

    def test_retry_stops_after_container_appears(self):
        container = ("project", "missing:1", "org/project/missing:1")
        failure = subprocess.CalledProcessError(1, ["singularity", "pull"])
        with patch.object(
            SPYTHON_FUNCTIONS, "container_pull", side_effect=[failure, None]
        ) as pull:
            with patch.object(
                SPYTHON_FUNCTIONS, "check_containers", side_effect=[[container], []]
            ):
                SPYTHON_FUNCTIONS.containers_routine_pull(
                    [container], "/containers", [container]
                )
        self.assertEqual(pull.call_count, 2)


if __name__ == "__main__":
    unittest.main()
