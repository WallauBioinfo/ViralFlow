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
  sudo apt install curl git python3-pip uidmap squashfs-tools -y
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

O ViralFlow disponibiliza um arquivo de ambiente por arquitetura: `envs/amd64.yml` (x86_64) e `envs/arm64.yml` (ARM64/aarch64). Verifique a arquitetura da sua máquina com `uname -m` (`x86_64` → AMD64, `aarch64` → ARM64).

#### AMD64

```bash
git clone https://github.com/WallauBioinfo/ViralFlow
cd ViralFlow/
micromamba env create -f envs/amd64.yml
micromamba activate viralflow
pip install -e .
```

#### ARM64

Em ARM64, o Singularity não está disponível via conda e deve ser compilado a partir do código-fonte antes da instalação do ViralFlow (o micromamba também deve ser baixado de `linux-aarch64` em vez de `linux-64` na etapa anterior):

```bash
sudo apt update
sudo apt install -y build-essential git wget pkg-config \
    libseccomp-dev squashfs-tools cryptsetup \
    libglib2.0-dev uuid-dev libssl-dev libgpgme-dev \
    libarchive-dev runc golang

git clone --recursive https://github.com/sylabs/singularity.git
cd singularity
git checkout v3.11.4
git submodule update --init --recursive
./mconfig
make -C builddir
sudo make -C builddir install
cd ..

git clone https://github.com/WallauBioinfo/ViralFlow
cd ViralFlow/
micromamba env create -f envs/arm64.yml
micromamba activate viralflow
pip install -e .
```

### Construindo os Containers

Todas as etapas do ViralFlow são executadas em ambientes controlados. O ViralFlow possui um método próprio para realizar toda essa construção de ambientes, rodando apenas uma linha de código.

Para que a construção dos containers ocorra, o ViralFlow necessita que as ferramentas "unsquashfs" e "mksquashfs" estejam disponíveis no diretório `/usr/local/bin/`. Para garantir isto, crie os links simbólicos:

```bash
sudo ln -sf /usr/bin/unsquashfs /usr/local/bin/unsquashfs
sudo ln -sf /usr/bin/mksquashfs /usr/local/bin/mksquashfs
```

Após garantir que ambas as ferramentas estejam no local apropriado, rode o comando para a construção dos containers:

```bash
# AMD64
viralflow build-containers

# ARM64
viralflow build-containers --arch arm64
```

```{note}
Este processo irá baixar aproximadamente 4,4GB de imagens de containers. Certifique-se de ter espaço suficiente em disco e uma conexão de internet estável.
```

## Solução de Problemas: AppArmor no Ubuntu 24.04+

A partir do Ubuntu 24.04, o AppArmor bloqueia por padrão os user namespaces sem privilégio. O Singularity instalado via conda depende deles, então o `viralflow build-containers` ou o `viralflow run` podem falhar com erros mencionando `permission denied` ou `user namespace`.

Verifique se a restrição está ativa (`1` indica que está):

```bash
sysctl kernel.apparmor_restrict_unprivileged_userns
```

Caso esteja, crie um perfil do AppArmor que libere os user namespaces apenas para o binário do Singularity. Rode os comandos abaixo com o ambiente `viralflow` ativado:

```bash
micromamba activate viralflow
STARTER=$(realpath $(dirname $(which singularity))/../libexec/singularity/bin/starter)
sudo tee /etc/apparmor.d/singularity-ce <<EOT
abi <abi/4.0>,
include <tunables/global>
profile singularity-ce ${STARTER}{,-suid} flags=(unconfined) {
  userns,
}
EOT
sudo systemctl reload apparmor
```

```{note}
O perfil está atrelado ao caminho do binário do Singularity. Caso você recrie o ambiente `viralflow` em outro local, rode estes comandos novamente.
```
