"""Guard the places a container version is written down more than once.

Nanopore_baseContainer.sing and nanopore_base.Dockerfile describe one
environment two ways, and the pinned Clair3 digest appears in three files.
Nothing at build or run time forces any of them to stay in step, so a version
bumped in one place and not the others would silently produce two different
containers - the exact drift the pins were added to prevent.
"""

import re
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONTAINERS = PROJECT_ROOT / "vfnext" / "containers"
SINGULARITY_RECIPE = CONTAINERS / "Nanopore_baseContainer.sing"
DOCKER_RECIPE = CONTAINERS / "nanopore_base.Dockerfile"
PROFILES = PROJECT_ROOT / "vfnext" / "configs" / "profiles.config"
NEXTFLOW_CONFIG = PROJECT_ROOT / "vfnext" / "nextflow.config"
TEST_CONFIG = PROJECT_ROOT / "vfnext" / "tests" / "nextflow.config"
TRUTH_TEST = PROJECT_ROOT / "vfnext" / "integration_tests" / "nanopore-truth.nf.test"
TRUTH_FIXTURE = (
    PROJECT_ROOT
    / "vfnext"
    / "tests"
    / "integration"
    / "data"
    / "nanopore_truth"
    / "expected_containers.tsv"
)

PINNED = (
    "HTSLIB_VERSION",
    "MINIMAP2_TAG",
    "NETWORKX_VERSION",
    "PORECHOP_ABI_COMMIT",
    "BAMUTIL_COMMIT",
    "LIBSTATGEN_COMMIT",
)


def singularity_pins():
    text = SINGULARITY_RECIPE.read_text()
    pins = {}
    for name in PINNED:
        match = re.search(rf"^\s*{name}=(\S+)\s*$", text, flags=re.MULTILINE)
        if match:
            pins[name] = match.group(1)
    return pins


def docker_pins():
    text = DOCKER_RECIPE.read_text()
    pins = {}
    for name in PINNED:
        match = re.search(rf"^ARG\s+{name}=(\S+)\s*$", text, flags=re.MULTILINE)
        if match:
            pins[name] = match.group(1)
    return pins


def docker_profile_image():
    text = PROFILES.read_text()
    match = re.search(r"base_container\s*=\s*['\"]([^'\"]+)['\"]", text)
    return match.group(1) if match else None


def manifest_version():
    text = NEXTFLOW_CONFIG.read_text()
    match = re.search(
        r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]", text, flags=re.MULTILINE
    )
    return match.group(1) if match else None


def clair3_pins():
    """Every place the pinned Clair3 image is written down.

    nextflow.config is what production runs; the other three exist for testing.
    The truth test declares it in its own params block because nf-test exposes
    only that block to a `then` block, so the fixture cross-check has nothing
    else to compare against. That block also outranks tests/nextflow.config, so
    without this guard a digest bumped in the config would leave the truth test
    quietly running the previous image while still passing.
    """
    pins = {}

    match = re.search(
        r"^\s*clair3_container\s*=\s*['\"]([^'\"]+)['\"]",
        NEXTFLOW_CONFIG.read_text(),
        flags=re.MULTILINE,
    )
    if match:
        pins["nextflow.config"] = match.group(1)

    match = re.search(
        r"params\.clair3_container\s*=\s*['\"]([^'\"]+)['\"]",
        TEST_CONFIG.read_text(),
    )
    if match:
        pins["tests/nextflow.config"] = match.group(1)

    match = re.search(
        r"clair3_container\s*=\s*['\"]([^'\"]+)['\"]", TRUTH_TEST.read_text()
    )
    if match:
        pins["nanopore-truth.nf.test"] = match.group(1)

    for line in TRUTH_FIXTURE.read_text().splitlines()[1:]:
        fields = line.split("\t")
        if fields and fields[0] == "clair3":
            pins["expected_containers.tsv"] = fields[1]

    return pins


class Clair3PinTests(unittest.TestCase):
    def test_every_location_declares_the_pin(self):
        self.assertEqual(
            sorted(clair3_pins()),
            [
                "expected_containers.tsv",
                "nanopore-truth.nf.test",
                "nextflow.config",
                "tests/nextflow.config",
            ],
        )

    def test_all_locations_agree(self):
        pins = clair3_pins()
        self.assertEqual(
            len(set(pins.values())),
            1,
            f"Clair3 is pinned inconsistently across files: {pins}",
        )

    def test_the_pin_is_a_digest_not_a_mutable_tag(self):
        for where, pin in clair3_pins().items():
            self.assertIn("@sha256:", pin, f"{where} pins a mutable tag: {pin}")


class ContainerRecipeTests(unittest.TestCase):
    def test_both_recipes_exist(self):
        self.assertTrue(SINGULARITY_RECIPE.is_file(), SINGULARITY_RECIPE)
        self.assertTrue(DOCKER_RECIPE.is_file(), DOCKER_RECIPE)

    def test_every_pin_is_declared_in_both_recipes(self):
        self.assertEqual(
            sorted(singularity_pins()), sorted(PINNED), "Singularity recipe"
        )
        self.assertEqual(sorted(docker_pins()), sorted(PINNED), "Dockerfile")

    def test_pins_match_between_recipes(self):
        self.assertEqual(
            singularity_pins(),
            docker_pins(),
            "Nanopore_baseContainer.sing and nanopore_base.Dockerfile pin "
            "different versions; update both together",
        )

    def test_docker_profile_image_tag_tracks_the_pipeline_version(self):
        image = docker_profile_image()
        self.assertIsNotNone(image, "docker profile does not set base_container")
        self.assertIn(":", image, f"image reference is not tagged: {image}")
        self.assertEqual(
            image.rsplit(":", 1)[1],
            manifest_version(),
            "the docker profile image tag and the pipeline manifest version "
            "have diverged",
        )


if __name__ == "__main__":
    unittest.main()
