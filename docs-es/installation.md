# Instalación

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
viralflow -build_containers
```

```{note}
Este proceso descargará aproximadamente 4.4GB de imágenes de contenedores. Asegúrese de tener suficiente espacio en disco y una conexión de internet estable.
```
