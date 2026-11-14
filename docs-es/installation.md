# Instalación

## Ubuntu 26.04

Esta branch utiliza Micromamba 1.5.7, Nextflow 22.04.0 y Java 17.
ViralFlow llama a `singularity` para construir, descargar, mantener y ejecutar
los contenedores, como en upstream. La GUI instala Apptainer y proporciona su
comando de compatibilidad `singularity`. Pangolin y snpEff
permanecen como **directorios sandbox**, a pesar de su sufijo `.sif`; las demás
imágenes son SIF.

La GUI de ViralFlow administra la instalación en pasos separados: Micromamba,
Nextflow, Apptainer, ViralFlow y contenedores. Las rutas predeterminadas son
`~/ViralFlowGUI/bin`, `~/ViralFlowGUI/micromamba` y `~/ViralFlowGUI/ViralFlow`.
Esta branch no proporciona `install.sh`.

Para instalar manualmente, proporcione Micromamba y un comando `singularity`
funcional en el host (Singularity o el comando de compatibilidad de Apptainer),
colóquelos en `PATH` y cree el entorno con el YAML de la arquitectura:

```bash
git clone -b feat/develop/ubuntu26.04 https://github.com/camilodotto/ViralFlow.git
cd ViralFlow
micromamba env create -n viralflow -f envs/amd64.yml -y
micromamba run -n viralflow python -m pip install -e .
NXF_VER=22.04.0 micromamba run -n viralflow nextflow -version
micromamba run -n viralflow viralflow --version
```

Java 17 y Nextflow 22.04.0 están definidos en el YAML, como en el entorno
original. La GUI también instala su launcher Nextflow en `bin` y selecciona
esa ruta absoluta. El runtime de contenedores está instalado en el host;
el YAML no instala un segundo runtime en el entorno Conda.
Para actualizar un entorno existente, utilice
`micromamba update -n viralflow -f envs/amd64.yml --prune -y`.
La GUI selecciona su propio binario Nextflow y mantiene la caché bajo la raíz
de la instalación. Consulte el [registro de migración](../UBUNTU26.04-CHANGES.md)
para cambios y validación. Las instrucciones siguientes describen las
instalaciones anteriores.

## Instalación en MacOS

Considerando las limitaciones del uso del Singularity en el MacOS, para rodar el ViralFlow en este tipo de sistema se sugiere utilizar un software de virtualización Linux llamado [Lima](https://github.com/lima-vm/lima).

### Instalando el Lima

Instale el Lima utilizando el administrador de paquetes Homebrew:

```bash
brew install lima
```

### Instalando la instancia Ubuntu

Siga el paso a paso de la instalación de la instancia Ubuntu:

```bash
limactl start
```

### Inicie el Ubuntu

Cuando la máquina virtual sea instalada, utilice el comando a continuación para iniciar el Ubuntu y siga los pasos de instalación del ViralFlow normalmente:

```bash
lima
```

## Instalación en Ubuntu

Para realizar la instalación del ViralFlow, deben seguirse cuatro pasos:
1. Instalar las dependencias del sistema
2. Instalar el Conda
3. Instalar el ViralFlow
4. Montar los contenedores para los análisis

Este proceso es realizado una única vez.

La instalación y uso del ViralFlow ya fue probada con los sistemas:
- Ubuntu 20.04 LTS
- Ubuntu 22.04 LTS

### Instalando las Dependencias del Sistema

En el caso de no tener el instalador de dependencias pip, el sistema de control de versiones git y el paquete uidmap, usted debe realizar la instalación con las siguientes líneas de código:

```bash
sudo apt update -y && \
  sudo apt upgrade -y && \
  sudo apt install curl git python3-pip uidmap -y
```

### Instalando y Configurando el Ambiente Micromamba

Recomendamos el administrador de ambientes micromamba por su paralelización al realizar el download e instalación de las dependencias. En el caso de que usted utilice otro administrador de ambientes (conda, miniconda, mamba), puede continuar con la instalación del ViralFlow ignorando esta etapa, sin embargo, pueden ocurrir conflictos durante la instalación.

En el caso de no tener el administrador de ambientes micromamba instalado:

```bash
cd $HOME
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/1.5.7 | tar -xvj bin/micromamba
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
micromamba activate
```

### Instalando el ViralFlow

Si usted ya dispone de las dependencias citadas y el conda instalado, podrá instalar el ViralFlow con 5 líneas de código:

```bash
git clone https://github.com/WallauBioinfo/ViralFlow
cd ViralFlow/
micromamba env create -f envs/env.yml
micromamba activate viralflow
pip install -e .
```

### Construyendo los Contenedores

Todas las etapas del ViralFlow son ejecutadas en ambientes controlados. El ViralFlow posee un método propio para realizar toda esa construcción de ambientes, rodando apenas una línea de código.

Para que la construcción de los contenedores ocurra, el ViralFlow necesita que la herramienta "unsquashfs" se encuentre disponible en el directorio `/usr/local/bin/`. Para garantizar esto, cree un link simbólico:

```bash
sudo ln -s /usr/bin/unsquashfs /usr/local/bin/unsquashfs
```

Después de garantizar que el unsquashfs se encuentra en el lugar adecuado, ejecuta el comando para la construcción de los contenedores:

```bash
viralflow build_containers
```

```{note}
Este proceso descargará aproximadamente 4.4GB de imágenes de contenedores. Asegúrese de tener suficiente espacio en disco y una conexión de internet estable.
```
