import subprocess
import os
import sys
import shutil

arch = sys.argv[1]

containers = [
    "pangolin:4.4.sif",
    "snpeff:5.0.sif",
]

container_commands = [
    f"singularity build -F --fakeroot --sandbox pangolin:4.4.sif def_files/{arch}/Singularity_pangolin",
    f"singularity build -F --fakeroot --sandbox snpeff:5.0.sif def_files/{arch}/Singularity_snpEff"
]

failed_containers = []
already_built = []
success = True  # Control variable

def container_exists(container):
    return os.path.isdir(container)

def build_container(container, command):
    print(f"@ Building {container}...")

    if container_exists(container):
        print(f"  > Container already exists. If you desire to rebuild the container, use the following command on 'ViralFlow/vfnext/containers/' directory:")
        print(f"    > {command}")
        already_built.append(container)
        return True

    try:
        subprocess.check_call(command, shell=True)
        print(f"  > Done <")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  > Failed <")
        print(f"Error: {e}")
        failed_containers.append((container, command))
        return False

print("Building containers:")

for container, command in zip(containers, container_commands):
    if not build_container(container, command):
        success = False

print("\nSummary:")

if failed_containers:
    print("\nSome containers failed to build. Try to build the containers on 'ViralFlow/vfnext/containers/' directory. Here are the details:")
    for container, command in failed_containers:
        print(f"\nContainer {container} failed to build.")
        print(f"Suggested command to build {container}:")
        print(f"  {command}")
    success = False

if success:
    print("\nAll containers were successfully built.")
else:
    print("\nSome containers failed to build. Please check the error messages above.")

if success:
    print("\nExecuting additional steps:\n")

    print("  > Loading sars-cov2 nextclade dataset...\n")
    nextclade_command = "singularity exec -B nextclade_dataset/sars-cov-2:/tmp nextclade:3.18.sif nextclade dataset get --name 'sars-cov-2' --output-dir '/tmp'"
    try:
        subprocess.check_call(nextclade_command, shell=True)
        print("    > Done <\n")
    except subprocess.CalledProcessError as e:
        print("    > Failed <")
        print(f"Error: {e}")
        success = False

    print("  > Downloading snpeff database catalog...")
    snpeff_command = "singularity exec snpeff:5.0.sif snpEff databases > snpEff_DB.catalog"
    try:
        subprocess.check_call(snpeff_command, shell=True)
        print("    > Done <")
    except subprocess.CalledProcessError as e:
        print("    > Failed <")
        print(f"Error: {e}")
        success = False

    if not shutil.which("unsquashfs"):
        print("Error: unsquashfs was not found in PATH. Install squashfs-tools.")
        success = False

if success:
    print("\nAll steps from '-build_containers' completed successfully. You can test ViralFlow using the following command:")
    print("   > viralflow run --params-file test_files/sars-cov-2.params")

sys.exit(0 if success else 1)
