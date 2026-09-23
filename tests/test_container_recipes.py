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
CONTAINERS_CONFIG = PROJECT_ROOT / "vfnext" / "configs" / "containers.config"
METADATA_TEST = (
    PROJECT_ROOT / "vfnext" / "tests" / "workflows" / "metadata-fixture.nf.test"
)
NEXTFLOW_CONFIG = PROJECT_ROOT / "vfnext" / "nextflow.config"
METADATA_HELPERS = PROJECT_ROOT / "vfnext" / "modules" / "metadata_helpers.nf"
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
    "BAMDASH_VERSION",
    "KALEIDO_VERSION",
    "PLOTLY_VERSION",
    "PYSAM_VERSION",
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


INTRAHOST_SCRIPT = PROJECT_ROOT / "vfnext" / "bin" / "intrahost.py"
RUFF_CONFIG = PROJECT_ROOT / "ruff.toml"

# The Python that intrahost_analysis:1.1.0.sif ships, from the def file that
# built it (vfnext/containers/def_files/intrahost_analysis.def at 5780713,
# deleted in c2be157 when the images moved to the Sylabs library).
INTRAHOST_CONTAINER_PYTHON = "py38"


class IntrahostContainerTests(unittest.TestCase):
    """intrahost.py has no container of its own; it borrows runReadCounts'.

    That container pins Python 3.8, which is older than everything else in the
    fleet, so the formatter can silently put syntax in the file that the
    container cannot parse. It already did: left at the repo-wide py312 target,
    ruff rewrote a multi-context `with` into the 3.10+ parenthesized form.
    """

    def test_ruff_targets_the_container_python(self):
        match = re.search(
            r"\[per-file-target-version\][^\[]*?"
            r'"vfnext/bin/intrahost\.py"\s*=\s*"([^"]+)"',
            RUFF_CONFIG.read_text(),
            flags=re.DOTALL,
        )
        self.assertIsNotNone(
            match, "ruff.toml no longer pins a target version for intrahost.py"
        )
        self.assertEqual(match.group(1), INTRAHOST_CONTAINER_PYTHON)

    def test_no_parenthesized_context_managers(self):
        """The one construct the formatter reaches for that 3.8 rejects.

        The pre-commit hook compiles the file under a real 3.8 and is the
        authority; this runs everywhere and names the likely cause.
        """
        offenders = [
            number
            for number, line in enumerate(
                INTRAHOST_SCRIPT.read_text().splitlines(), start=1
            )
            if re.match(r"\s*with\s*\($", line)
        ]
        self.assertEqual(
            offenders,
            [],
            f"parenthesized `with` is Python 3.10+, at line(s) {offenders}",
        )


def illumina_container_images():
    """The name -> path entries of params.illumina_containers."""
    text = CONTAINERS_CONFIG.read_text()
    block = re.search(r"illumina_containers\s*=\s*\[(.*?)\]", text, re.DOTALL)
    assert block, "params.illumina_containers not found in containers.config"
    return dict(re.findall(r"(\w+)\s*:\s*[\"'](.+?)[\"']", block.group(1)))


def illumina_container_references():
    """The image names the process directives resolve, in file order.

    Read from inside each container closure rather than matching the whole
    closure, because the GENPLOTS directives choose between the nanopore base
    image and an ILLUMINA one depending on the mode.
    """
    return [
        name
        for body in container_closure_bodies(CONTAINERS_CONFIG)
        for name in re.findall(
            r"params\.getOrDefault\(\s*'illumina_containers'\s*,\s*\[:\]\s*\)\.(\w+)",
            body,
        )
    ]


def container_closure_bodies(path):
    """The source of every `container = { ... }` closure in a config file."""
    text = path.read_text()
    bodies = []
    for match in re.finditer(r"\bcontainer\s*=\s*\{", text):
        depth, start = 1, match.end()
        for index in range(start, len(text)):
            depth += {"{": 1, "}": -1}.get(text[index], 0)
            if depth == 0:
                bodies.append(text[start:index])
                break
    return bodies


class ContainerConfigTests(unittest.TestCase):
    def test_illumina_processes_use_local_images(self):
        """Every ILLUMINA image is a local .sif under containers/.

        runIntraHostScript was briefly an exception, pointing at a Wave image
        pulled at run time. That image is published for linux/amd64 only, while
        the project ships arm64 SIFs too, and it is absent from
        containers/repositories/, so `viralflow build-containers` never fetched
        it. Clair3 is the one deliberate registry dependency and lives in
        nextflow.config, not here.
        """
        images = illumina_container_images()
        self.assertTrue(images, "no images declared in params.illumina_containers")
        for name, image in images.items():
            self.assertRegex(
                image,
                r"^\$projectDir/containers/",
                f"non-local container declared for {name}: {image}",
            )

    def test_no_process_hardcodes_a_container_path(self):
        """Directives resolve through the map rather than naming an image.

        A literal path here would be a second place to bump a version, which is
        what let containers.config and metadata_helpers.containerSpecs() drift
        apart: the run used one image while container_manifest.tsv named
        another. It would also slip past the local-image check above, which now
        reads the map.
        """
        literals = [
            line.strip()
            for line in CONTAINERS_CONFIG.read_text().splitlines()
            if re.search(r"^\s*container\s*=\s*[\"']", line)
        ]
        self.assertEqual(
            literals,
            [],
            "container directives must read params.illumina_containers",
        )

    def test_container_closures_do_not_read_params_directly(self):
        """Container closures read params through getOrDefault.

        Nextflow evaluates every container closure once at startup, to fill
        workflow.container, before it binds the config params. A plain
        params.x read there finds nothing and prints "Access to undefined
        parameter" on every run - harmless, since tasks resolve the closure
        again later, but it reads like a misconfiguration to anyone running
        the pipeline. getOrDefault returns the same value once params are bound
        and stays quiet before.
        """
        configs = [
            path
            for path in sorted((PROJECT_ROOT / "vfnext").rglob("*.config"))
            if ".nf-test" not in path.parts
        ]
        offenders = [
            f"{path.relative_to(PROJECT_ROOT)}: {body.strip()}"
            for path in configs
            for body in container_closure_bodies(path)
            if re.search(r"\bparams\.(?!getOrDefault\()", body)
        ]
        closures = sum(len(container_closure_bodies(path)) for path in configs)
        self.assertGreater(closures, 20, "container closures not found")
        self.assertEqual(offenders, [])

    def test_every_referenced_image_is_declared(self):
        """No directive resolves a name the map does not declare.

        Nextflow would hand the task a null container and run it on the host,
        silently, rather than failing.
        """
        images = illumina_container_images()
        references = illumina_container_references()
        self.assertTrue(references, "no process reads params.illumina_containers")
        for name in sorted(set(references)):
            self.assertIn(
                name,
                images,
                f"containers.config resolves undeclared image '{name}'",
            )

    def test_every_declared_image_is_used(self):
        """Nothing is declared that no process runs and no manifest reports.

        An unused entry is either a process that lost its container or an image
        the project stopped shipping; both are worth noticing rather than
        carrying.
        """
        images = illumina_container_images()
        used = set(illumina_container_references()) | set(
            re.findall(r"'(\w+)'", METADATA_HELPERS.read_text())
        )
        for name in sorted(images):
            self.assertIn(
                name,
                used,
                f"'{name}' is declared but no process or manifest entry uses it",
            )


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
