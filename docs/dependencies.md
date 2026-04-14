# Dependencies

To install and run ViralFlow, you will need at least 4.4GB of disk storage (to build containers), and three system dependencies:

- pip
- Git
- uidmap

The list with tool versions used in ViralFlow can be checked [at this link](https://github.com/WallauBioinfo/ViralFlow/tree/main/versions).

Moreover, for an installation with a controlled environment, we use Conda. ViralFlow uses the Singularity containerization tool and the Nextflow workflow language. To ensure that ViralFlow works correctly, the correct versions of Singularity and Nextflow are essential, both of which can be installed together with ViralFlow.

## Conda

Conda is an open source package manager. It was initially developed to manage Python language package environments, allowing version control.

So that all users can use ViralFlow with the exact versions of Singularity and Nextflow, we created an environment with these versions, which allows a standardization in the ViralFlow installation environments, as well as the resolution of eventual errors.

## Singularity v3.11.4

Singularity is a tool that performs virtualization at the operating system level, a task popularly known as containerization.

For a computational pipeline like ViralFlow to work, it requires several dependencies (other tools and libraries of the operating system) that can be complex to install. In order to make things easier to the user we use the container system which already comes with all dependencies needed to run ViralFlow. Therefore, the user does not need to install everything on his own.

With containerization, the user only needs to install one tool, in this case Singularity, and then build the containers based on recipe files (in the case of ViralFlow, this was also automated with the setupContainers.sh script). In this way, different research groups can run the tool with the same computational environment, without changes in ViralFlow's behavior, which guarantees the reproducibility of the results.

## Nextflow 22.04

Nextflow is a bioinformatics workflow manager that allows the development of portable and reproducible workflows.

Nextflow allows the execution of workflows in different computing environments, such as on-premises, high-performance computing (HPC) environments, cloud computing services, such as AWS and Google Cloud, and also in container orchestration systems, such as Kubernetes.

In addition, Nextflow provides support for managing dependencies with Conda, Docker, Podman or Singularity.
