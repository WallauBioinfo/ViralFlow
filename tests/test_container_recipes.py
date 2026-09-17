"""Keep the Singularity and Docker nanopore recipes in agreement.

Nanopore_baseContainer.sing and nanopore_base.Dockerfile describe one
environment two ways. Nothing at build time forces them to stay in step, so a
pin bumped in one and not the other would silently produce two different
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

PINNED = (
    "HTSLIB_VERSION",
    "MINIMAP2_TAG",
    "NETWORKX_VERSION",
    "PORECHOP_ABI_COMMIT",
    "BAMUTIL_COMMIT",
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
