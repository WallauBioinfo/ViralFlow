import spython_functions
import sys

arch = sys.argv[1]
repositories =  {
    "arm64": "repositories_arm64.txt",
    "amd64": "repositories_amd64.txt"
}

containers_dir = spython_functions.os.getcwd()
repositories = f"{containers_dir}/repositories/{repositories.get(arch)}"
containers_names_list = spython_functions.get_repository(repositories)
missing_containers = spython_functions.check_containers(containers_names_list, containers_dir)
if missing_containers:
    spython_functions.containers_routine_pull(missing_containers, containers_dir, containers_names_list)
else:
    print(f":: All containers were sucessfully downloaded! :: ")
