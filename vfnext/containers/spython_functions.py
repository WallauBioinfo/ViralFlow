from spython.main import Client
import os


def get_repository(repository_list):
    containers_name_list = []
    with open(repository_list, "r") as containers:
        file = containers.readlines()
        for container in file:
            container = container.strip()
            if not container:
                continue
            # Full path: e.g. wallaulabs2/viralflow-amd64/edirect:1.1.0
            parts = container.split("/")
            org = parts[0]  # wallaulabs2
            project = parts[1]  # viralflow-amd64
            container_version = parts[2]  # edirect:1.1.0
            full_repo = f"{org}/{project}/{container_version}"
            containers_name_list.append((project, container_version, full_repo))
    return containers_name_list


def container_pull(containers_dir, containers_name_list):
    for container in containers_name_list:
        container_version = container[1]
        full_repo = container[2]
        print(f"Downloading container {container_version}. This could be take a while. Please Wait ...")
        os.system(f"singularity pull -F {container_version}.sif library://{full_repo}")  



def check_containers(containers_name_list, downloaded_list):
    download_list = os.scandir(downloaded_list)
    container_download_list = [container.name for container in download_list]
    missing_containers = []
    for cont in containers_name_list:
        # cont is (project, container_version, full_repo)
        expected_file = f"{cont[1]}.sif"
        if expected_file not in container_download_list:
            missing_containers.append(cont)
    return missing_containers


def containers_routine_pull(missing_containers_list, containers_dir, containers_names_list):
    lost_containers = missing_containers_list
    attempts = len(lost_containers) * 3
    while len(lost_containers) > 0 or attempts > 0:
        container_pull(containers_dir, missing_containers_list)
        lost_containers = check_containers(containers_names_list, containers_dir)
        attempts -= 1
        if len(lost_containers) == 0:
            break

