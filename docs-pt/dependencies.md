# Dependências

Para instalação e execução do ViralFlow, você precisará de pelo menos 4,4GB de armazenamento de disco, para a criação dos contêineres, e de três dependências:

- pip
- Git
- uidmap

A lista com as versões de cada ferramenta utilizadas no ViralFlow pode ser conferida [neste link](https://github.com/WallauBioinfo/ViralFlow/tree/main/versions).

Além disso, para uma instalação com um ambiente controlado, utilizamos o Conda. O ViralFlow utiliza a ferramenta de containerização Singularity e da linguagem de workflow Nextflow. Para garantir o funcionamento correto do ViralFlow, é essencial as versões corretas do Singularity e do Nextflow, ambos podem ser instalados juntamente com o ViralFlow.

## Conda

Conda é um gerenciador de pacotes de código aberto. Foi desenvolvido inicialmente para gerenciar ambientes de pacotes da linguagem Python, permitindo o controle de versões.

Para que todos os usuários utilizem o ViralFlow com as versões exatas do Singularity e do Nextflow, nós criamos um ambiente com estas versões, o que permite uma padronização nos ambientes de instalação do ViralFlow, bem como a resolução de eventuais erros.

## Singularity v3.11.4

Singularity é uma ferramenta que executa a virtualização no nível do sistema operacional, tarefa popularmente conhecida como conteinerização.

Para uma pipeline computacional, como o ViralFlow, funcionar, ela precisa de várias dependências (outras ferramentas e bibliotecas do sistema operacional) e para evitar que o usuário precise instalar todas as dependências manualmente em sua máquina, nós usamos o sistema de containers.

Com a conteinerização, o usuário só precisa instalar uma ferramenta, neste caso o Singularity, e depois construir os containers com base em arquivos de receita (no caso do ViralFlow, isto também foi automatizado com o script setupContainers.sh). Dessa forma, diferentes grupos de pesquisa conseguem rodar a ferramenta com o mesmo ambiente computacional, sem mudanças no comportamento do ViralFlow, o que garante a reprodutibilidade dos resultados.

## Nextflow 22.04

Nextflow é um gerenciador de workflows de bioinformática que permite o desenvolvimento de workflows portáteis e reprodutíveis.

Nextflow permite a execução dos workflows em diversos ambientes computacionais, como ambientes locais, de computação de alta performance (HPC), serviços de computação na nuvem, como AWS e Google Cloud, e também em sistemas de orquestração de containers, como Kubernetes.

Além disso, o Nextflow provê suporte para gerenciamento de dependências com Conda, Docker, Podman ou Singularity.
