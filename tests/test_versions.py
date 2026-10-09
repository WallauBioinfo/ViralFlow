import ast
import re
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def python_package_version():
    setup_tree = ast.parse((PROJECT_ROOT / "setup.py").read_text())
    setup_call = next(
        node
        for node in ast.walk(setup_tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "setup"
    )
    version_keyword = next(
        keyword for keyword in setup_call.keywords if keyword.arg == "version"
    )
    return ast.literal_eval(version_keyword.value)


def nextflow_manifest_version():
    config = (PROJECT_ROOT / "vfnext" / "nextflow.config").read_text()
    match = re.search(
        r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]",
        config,
        flags=re.MULTILINE,
    )
    if match is None:
        raise AssertionError("Nextflow manifest version was not found")
    return match.group(1)


class VersionConsistencyTests(unittest.TestCase):
    def test_python_and_nextflow_versions_match(self):
        self.assertEqual(python_package_version(), nextflow_manifest_version())


if __name__ == "__main__":
    unittest.main()
