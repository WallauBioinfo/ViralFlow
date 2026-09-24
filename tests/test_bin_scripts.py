"""vfnext/bin scripts are called by name, the way Nextflow expects.

Nextflow puts bin/ on every task's PATH, uploading it first on executors that
do not share the launch host's filesystem. A module that calls
${projectDir}/bin/script.py names the launch host's copy instead, which works
locally and fails on a cloud executor. Calling a script by name in turn needs
it executable with a shebang, or the task fails with "Permission denied".
"""

import os
import re
import stat
import unittest
from pathlib import Path

VFNEXT = Path(__file__).resolve().parents[1] / "vfnext"
BIN = VFNEXT / "bin"
MODULES = VFNEXT / "modules"

# ILLUMINA modules still on the old form, left for the ILLUMINA branch
# (TODO.md, section 4). Remove each entry as it is fixed; a new one fails.
PROJECT_DIR_BIN_ALLOWED = {"runIntraHostScript.nf", "runIvar.nf"}

PROJECT_DIR_BIN = re.compile(r"\$\{?projectDir\}?/bin/")


def module_texts():
    return {path.name: path.read_text() for path in sorted(MODULES.glob("*.nf"))}


class BinScriptTests(unittest.TestCase):
    def test_modules_do_not_reach_into_project_dir_bin(self):
        offenders = {
            name
            for name, text in module_texts().items()
            if PROJECT_DIR_BIN.search(text)
        }
        self.assertEqual(offenders - PROJECT_DIR_BIN_ALLOWED, set())
        # Keeps the allowlist honest: a fixed module must leave it too.
        self.assertEqual(PROJECT_DIR_BIN_ALLOWED - offenders, set())

    def test_scripts_called_by_name_are_executable_with_a_shebang(self):
        texts = "\n".join(module_texts().values())
        called = [
            script
            for script in sorted(BIN.glob("*.py"))
            # Not preceded by a path: bin/x.py in a comment is not a call.
            if re.search(rf"(?<![/\w]){re.escape(script.name)}\b", texts)
        ]
        self.assertIn(BIN / "coverage_plot.py", called)
        self.assertIn(BIN / "nanopore_summary.py", called)
        for script in called:
            with self.subTest(script=script.name):
                self.assertTrue(
                    os.stat(script).st_mode & stat.S_IXUSR, "not executable"
                )
                self.assertTrue(script.read_text().startswith("#!"), "no shebang line")


if __name__ == "__main__":
    unittest.main()
