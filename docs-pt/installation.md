# Instalação

## Instalação no MacOS

Devido a limitações do uso do singularity no MacOS, para rodar o ViralFlow neste tipo de sistema, sugerimos a utilização de um software de virtualização Linux chamado [Lima](https://github.com/lima-vm/lima).

### Instalando o Lima

Instale o Lima utilizando o gerenciador de pacotes Homebrew:

```bash
brew install lima
```

### Instalando a instância Ubuntu

Siga o passo a passo da instalação de uma instância Ubuntu:

```bash
limactl start
```

### Inicie o Ubuntu

Quando a máquina virtual for instalada, utilize o comando abaixo para iniciar o Ubuntu e siga os passos de instalação do ViralFlow normalmente:

```bash
lima
```

## Instalação no Ubuntu

Para realizar a instalação do ViralFlow, são necessários quatro passos:
1. Instalar as dependências de sistema
2. Instalar o Conda
3. Instalar o ViralFlow
4. Montar os containers para as análises

Este processo é realizado uma única vez.

A instalação e uso do ViralFlow já foi testada nos sistemas:
- Ubuntu 20.04 LTS
- Ubuntu 22.04 LTS

### Instalando as Dependências do Sistema

Caso você não tenha o instalador de dependências pip, o sistema de controle de versões git, e o pacote uidmap, você deve realizar a instalação com as seguintes linhas:

```bash
sudo apt update -y && \
  sudo apt upgrade -y && \
  sudo apt install curl git python3-pip uidmap -y
```

### Instalando e Configurando o Ambiente Micromamba

Recomendamos o gerenciador de ambientes micromamba devido a sua paralelização ao fazer o download e instalação das dependências. Caso você utilize outro gerenciador de ambientes (conda, miniconda, mamba), você pode pular esta etapa, porém, conflitos durante a instalação podem ocorrer.

Caso você não tenha o gerenciador de ambientes micromamba instalado:

```bash
cd $HOME
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/1.5.7 | tar -xvj bin/micromamba
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
micromamba activate
```

### Instalando o ViralFlow

Caso você já tenha as dependências citadas e o conda instalado, você pode instalar o ViralFlow com 5 linhas de código:

```bash
git clone https://github.com/WallauBioinfo/ViralFlow
cd ViralFlow/
micromamba env create -f envs/env.yml
micromamba activate viralflow
pip install -e .
```

### Construindo os Containers

Todas as etapas do ViralFlow são executadas em ambientes controlados. O ViralFlow possui um método próprio para realizar toda essa construção de ambientes, rodando apenas uma linha de código.

Para que a construção dos containers ocorra, o ViralFlow necessita que a ferramenta "unsquashfs" esteja disponível no diretório `/usr/local/bin/`. Para garantir isto, crie um link simbólico:

```bash
sudo ln -s /usr/bin/unsquashfs /usr/local/bin/unsquashfs
```

Após garantir que o unsquashfs esteja no local apropriado, rode o comando para a construção dos containers:

```bash
viralflow -build_containers
```

```{note}
Este processo irá baixar aproximadamente 4,4GB de imagens de containers. Certifique-se de ter espaço suficiente em disco e uma conexão de internet estável.
```
