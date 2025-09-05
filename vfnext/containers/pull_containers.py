import spython_functions

containers_dir = spython_functions.os.getcwd()
repositories = containers_dir + '/repositories/repositories_viral_flow.txt'
containers_names_list = spython_functions.get_repository(repositories)
missing_containers = spython_functions.check_containers(containers_names_list, containers_dir)
if missing_containers:
    spython_functions.containers_routine_pull(missing_containers, containers_dir, containers_names_list)
else:
    print(f":: All containers were sucessfully downloaded! :: ")
