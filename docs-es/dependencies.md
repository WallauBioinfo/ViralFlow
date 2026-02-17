# Dependencias

Para la instalación y ejecución del ViralFlow, será necesario disponer de al menos 4,4 GB de memoria, para la elaboración de los contenedores, y de tres dependencias:

- pip
- Git
- uidmap

La lista con las versiones de cada herramienta utilizada en el ViralFlow puede ser verificada [en el siguiente link](https://github.com/WallauBioinfo/ViralFlow/tree/main/versions).

Adicionalmente, para la instalación con un ambiente controlado, se usa el Conda. El ViralFlow utiliza la herramienta de creación de contenedores Singularity y del lenguaje de workflow Nextflow. Para garantizar el funcionamiento correcto del ViralFlow es esencial tener las versiones adecuadas del Singularity y del Nextflow, y ambos pueden ser instalados simultáneamente con el ViralFlow.

## Conda

Conda es un administrador de paquetes de código abierto. Fue inicialmente desarrollado para gerenciar ambientes de paquetes de lenguaje Python, permitiendo el control de versiones.

Para que todos los usuarios utilicen el ViralFlow con las versiones exactas del Singularity y del Nextflow, hemos creado un ambiente con estas versiones, permitiendo una estandarización en los ambientes de instalación del ViralFlow, así como la resolución de eventuales errores.

## Singularity v3.11.4

Singularity es una herramienta que ejecuta la virtualización a nivel del sistema operacional, tarea popularmente conocida como creación de contenedores.

Para un pipeline computacional funcionar, como el ViralFlow, se necesita de varias dependencias (otras herramientas y bibliotecas del sistema operacional); y para evitar que el usuario necesite instalar todas las dependencias manualmente en su máquina, utilizamos el sistema de contenedores.

Con la creación de contenedores, el usuario solo necesita instalar una herramienta, en este caso el Singularity, y después construir los contenedores basándose en archivos de receta (en el caso del ViralFlow, también fue automatizado con el script setupContainers.sh). De esta forma, diferentes grupos de investigación pueden rodar la herramienta en el mismo ambiente computacional, sin cambios en el comportamiento del ViralFlow, garantizando la reproductibilidad de los resultados.

## Nextflow 22.04

Nextflow es un administrador de workflows de bioinformática que permite el desarrollo de workflows portátiles y reproducibles.

Nextflow permite la ejecución de los workflows en diversos ambientes computacionales, como ambientes locales, de computación de alta performance (HPC), servicios de computación en la nube, AWS y Google Cloud, y también, en sistemas de orquestación de contenedores como Kubernetes.

Adicionalmente, el Nextflow proporciona soporte para gerenciar las dependencias de Conda, Docker, Podman o Singularity.
