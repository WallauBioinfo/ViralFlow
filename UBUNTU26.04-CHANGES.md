# Diferenças em relação a WallauBioinfo/develop

## Resumo

A branch `feat/develop/ubuntu26.04` adapta a instalação e o ambiente de
execução do ViralFlow para Ubuntu 26.04, preservando Nextflow 22.04.0,
Java 17 e as chamadas ao comando `singularity` do upstream. Construção,
download, manutenção e execução dos containers usam esse comando; o ViralFlow
não exige um executável chamado `apptainer`. A GUI utiliza Micromamba 1.5.7,
instala o ambiente diretamente e fornece Apptainer no host com o comando de
compatibilidade `singularity`.

Pangolin e snpEff continuam como diretórios sandbox, embora seus nomes terminem
em `.sif`; as demais imagens permanecem em formato SIF. As diferenças de
execução tratam da construção e manutenção desses containers, dependências,
cache gravável, tratamento de erros e seleção explícita do Nextflow da
instalação escolhida pela GUI.

A comparação também inclui quatro gráficos SVG acrescentados à compilação
dos resultados em `vfnext/bin/compileOutput.py`. O fluxo principal, os módulos,
os workflows e os parâmetros científicos padrão permanecem iguais ao upstream.

A referência desta comparação é `WallauBioinfo/develop`, commit
`16cb3f30b1d3`. O documento descreve o conteúdo atual do diretório de trabalho.
São **18 arquivos diferentes:
14 modificados e 4 adicionados**. Nenhum arquivo da referência foi removido.

## Diferenças por arquivo

### `README.md` — modificado

Acrescenta um link para a seção de instalação em Ubuntu 26.04
e para este documento. O restante do README é preservado.

### `UBUNTU26.04-CHANGES.md` — adicionado

Documenta o estado atual da branch e suas diferenças por arquivo em relação
a `WallauBioinfo/develop`.

### `docs/installation.md` — modificado

Acrescenta uma seção em inglês para Ubuntu 26.04. Documenta Micromamba 1.5.7,
Nextflow 22.04.0, Java 17, o comando `singularity`, os sandboxes e a instalação
direta ou pela GUI. Distingue o comando utilizado pelo ViralFlow do runtime
Apptainer instalado pela GUI. Explica os
caminhos próprios da GUI, a criação/atualização do ambiente e a seleção do
launcher Nextflow em `bin`. As instruções preexistentes são preservadas.

### `docs-pt/installation.md` — modificado

Acrescenta a mesma orientação de instalação e execução em português,
preservando as instruções preexistentes.

### `docs-es/installation.md` — modificado

Acrescenta a mesma orientação de instalação e execução em espanhol,
preservando as instruções preexistentes.

### `envs/amd64.yml` — modificado

Remove `singularityce=3.11.4` e `spython=0.3.1`. O comando `singularity` é
fornecido pelo host; instalar SingularityCE também no ambiente Conda criaria
um segundo runtime, que poderia preceder o comando de compatibilidade do
Apptainer instalado pela GUI. `spython` não é utilizado pelos scripts; seus
downloads são feitos pelo CLI `singularity`, como já eram no upstream.
Mantém `nextflow=22.04.0`, `openjdk=17` e as demais dependências do upstream.

Reduz os canais a `bioconda` e `conda-forge`, retirando
`bioconda/label/cf201901`, `conda-forge/label/cf201901`,
`conda-forge/label/main`, `pkgs/main` e `wallaulab`. A lista reduzida evita
consultas a canais antigos, duplicados ou desnecessários para esse ambiente.
A consulta ao índice `noarch` do canal `pkgs/main`, conforme resolvido em
`conda.anaconda.org`, retornou HTTP 404; isso não significa que todos os
canais retirados estejam indisponíveis.

### `envs/arm64.yml` — modificado

Remove somente a dependência não utilizada `spython`. Os downloads utilizam
o CLI `singularity`. Preserva os canais, Nextflow 22.04.0, Java 17 e as demais
dependências do upstream.

### `tests/host-smoke.sh` — adicionado

Inclui um teste de execução da pipeline em modo `custom`, com referência
aleatória não biológica e leituras sintéticas pareadas e simples. Verifica
a criação dos consensos e do diretório de resultados compilados. Permite
selecionar o comando ViralFlow e o diretório de validação, preservando o
código e os limiares científicos da pipeline.

### `tests/pangolin-cache-smoke.sh` — adicionado

Executa duas tarefas mínimas Nextflow usando a configuração real de
`runPangolin` e um sandbox selecionado. Inicializa os caches de fontes e de
runtime do Snakemake e verifica escrita e isolamento entre tarefas. Não
executa a classificação Pangolin nem utiliza sequências biológicas.

### `tests/test_runtime.py` — adicionado

Inclui 12 testes para propagação de falhas de execução/download/manutenção,
seleção explícita do Nextflow e isolamento entre instalações, argumentos
com espaços, limpeza dos arquivos temporários, distinção entre os modos de
atualização Pangolin, preservação do catálogo snpEff diante de falhas,
seleção explícita da biblioteca de containers e término das tentativas de
download. Os testes de download e atualização verificam chamadas a
`singularity`; o executável simulado do teste de manutenção snpEff também
usa esse nome. Os testes usam mocks e executáveis simulados.

### `vfnext/bin/compileOutput.py` — modificado

Acrescenta quatro gráficos SVG aos resultados compilados:

- Contagem de leituras por amostra.
- Relação entre abrangência e profundidade média de cobertura.
- Distribuição da abrangência de cobertura por faixas de classificação.
- Distribuição da abrangência de cobertura por decis.

Inclui funções auxiliares para escala dos eixos, regressão da visualização,
conversão de proporções para porcentagens, escape do texto SVG e destaque dos
controles negativos. Os gráficos utilizam os dados das tabelas de resultados
e são gerados durante sua compilação. Essa diferença deve ser considerada na
comparação completa da branch, além das adaptações do ambiente de execução.

### `vfnext/configs/containers.config` — modificado

Mantém o backend `singularity`, habilitado com montagem automática, como no
upstream. Remove a seleção global de escrita por arquitetura: `--writable`
para AMD64 e `--writable-tmpfs` para ARM64. A escrita temporária passa a ser
restrita ao processo `runSnpEff`, evitando aplicar opções de escrita a todas
as imagens SIF.

Define `XDG_CACHE_HOME` no diretório de trabalho de cada tarefa `runPangolin`,
permitindo ao Snakemake criar seus caches em uma pasta gravável e independente
por tarefa. Preserva os nomes e caminhos dos containers atribuídos aos
processos.

### `vfnext/containers/add_entries_SnpeffDB.sh` — modificado

Preserva as chamadas `singularity` e acrescenta interrupção em caso de erro,
validação dos argumentos e tratamento de caminhos/nomes com
espaços. Obtém o registro GenBank completo em um diretório temporário e
verifica a presença de sequência antes de alterar a configuração do sandbox.

Evita duplicar uma entrada já existente na configuração, mantém a construção
do banco pelo snpEff e verifica a existência do arquivo de banco gerado.
Produz o catálogo em arquivo temporário e o publica após a execução
bem-sucedida, com limpeza dos arquivos temporários ao encerrar.

### `vfnext/containers/build_containers.py` — modificado

Mantém `singularity` na construção dos sandboxes Pangolin/snpEff e na execução
das etapas de preparação dos containers. Preserva as opções de construção
sandbox e fakeroot e as versões declaradas das ferramentas.

Verifica `unsquashfs` pelo `PATH`, em vez de exigir sua presença em
`/usr/local/bin/unsquashfs`, e devolve código de saída diferente de zero quando
alguma etapa falha. Isso permite à GUI reconhecer a falha do build.

### `vfnext/containers/def_files/amd64/Singularity_pangolin` — modificado

Troca a base Debian 11 obtida pela biblioteca pela base `debian:12-slim`
obtida via transporte Docker, corrigindo a construção diante dos
problemas de disponibilidade dos pacotes da base anterior. Esse transporte
não exige Docker Engine para executar a pipeline.

Acrescenta interrupção em caso de erro e instalação/verificação das dependências
Python declaradas por Pangolin e Snakemake, para conciliar essas dependências
com os pacotes Conda. Preserva as versões explícitas de Python, Pangolin,
Snakemake e Usher presentes no upstream e o Micromamba próprio do container.

### `vfnext/containers/def_files/amd64/Singularity_snpEff` — modificado

Preserva a imagem base e as versões das ferramentas. Redireciona os repositórios
APT da distribuição antiga para `archive.debian.org` e ajusta a verificação de
validade dos índices arquivados, permitindo instalar os pacotes necessários.
Acrescenta interrupção em caso de erro na preparação do container.

### `vfnext/containers/spython_functions.py` — modificado

Remove o import não utilizado de `spython` e preserva downloads pelo CLI
`singularity`, informando explicitamente a biblioteca Sylabs. Essa opção
continua necessária quando `singularity` executa Apptainer, pois não depende
da biblioteca padrão desse runtime. O download usa o diretório de containers
selecionado e propaga falhas, sem depender do remote padrão do usuário.

Corrige a condição do laço de tentativas, tenta somente os containers ainda
ausentes e informa erro se faltarem imagens ao final das tentativas.

### `wrapper/__init__.py` — modificado

Remove o import não utilizado de `distutils`, ausente nas versões atuais do
Python. O import de `logging` permanece igual ao upstream. Executa construção
e manutenção por subprocessos com
propagação de falhas e preservação dos argumentos/caminhos com espaços.
A construção utiliza o mesmo interpretador Python que executa o wrapper.

As operações de atualização Pangolin mantêm `singularity` e usam um diretório
temporário montado em `/tmp` e uma restrição de ferramentas de build
`setuptools<81`, necessária a pacotes que ainda utilizam `pkg_resources`.
A atualização completa também instala e verifica as dependências Python
declaradas por Pangolin e Snakemake. A atualização somente de dados conserva
sua operação própria.

Mantém Nextflow 22.04.0 como padrão, permitindo selecionar a versão por
`NXF_VER`. Aceita `VIRALFLOW_NEXTFLOW` como caminho absoluto para o executável
escolhido pela GUI, valida esse caminho e aplica quoting ao comando. Sem essa
variável, utiliza Nextflow pelo `PATH`, como na chamada CLI tradicional.
Propaga o código de falha do Nextflow. O parser dos parâmetros e a opção
`-resume` permanecem iguais ao upstream.

## Arquivos e comportamentos preservados

Os comandos de containers do wrapper e dos scripts auxiliares continuam
sendo `singularity`, assim como o backend do Nextflow. Não há seleção automática
de Apptainer ou exigência desse nome de executável no código de execução.
A GUI fornece o runtime e o link de compatibilidade em seu próprio repositório.

`vfnext/nextflow.config` e `vfnext/configs/profiles.config` são idênticos a
`WallauBioinfo/develop`. Não existe diferença de versão Nextflow nos YAMLs:
ambos mantêm 22.04.0, assim como o padrão do wrapper. Java permanece em 17.
Micromamba 1.5.7 é selecionado pela GUI; esta branch não possui um instalador
que fixe sua versão.

`install.sh` e `tests/install-script.test.sh` estão ausentes tanto nesta branch
quanto no upstream. Portanto, não aparecem como arquivos removidos nessa
comparação. A instalação é feita pela GUI ou pelas instruções manuais.

O fluxo principal, os módulos, os workflows, os parâmetros científicos padrão
e os dados de teste versionados permanecem iguais à referência. As diferenças
dos scripts de análise estão concentradas nos gráficos adicionais de
`compileOutput.py`. As definições de containers ARM64 também são preservadas.
A construção continua sem opção de limpeza prévia de containers.

## Resultado da revisão das diferenças

A revisão cobre todos os arquivos diferentes da referência. A substituição
do comando `singularity` por `apptainer` foi descartada. Também foi descartada
a remoção do import `logging`, que não era necessária à compatibilidade.
A documentação e os testes descrevem o uso de `singularity`.

As diferenças restantes do ambiente tratam de problemas verificados:
instalação de dependências/canais, repositórios APT das imagens, dependências
Python dos containers, cache gravável, manutenção dos sandboxes, propagação
de falhas e seleção do Nextflow próprio da GUI. Nenhum desses ajustes exige
o backend nativo Apptainer ou versões mais recentes de Nextflow/Micromamba.
As versões de Nextflow nos YAMLs e as configurações principal e de perfis
continuam iguais ao upstream.

Os quatro gráficos de `compileOutput.py` são a diferença preexistente que
não é necessária à instalação ou execução em Ubuntu 26.04. Foram preservados
porque pertencem à apresentação dos resultados, fora do escopo desta revisão
do ambiente. Os três arquivos de testes são verificações adicionais e não
são dependências para executar a pipeline.

## Validação disponível

A instalação inicial, a reinstalação e o uso de Micromamba 1.5.7 sobre um
ambiente criado com 2.9.0 foram verificados em diretórios isolados. Nextflow
22.04.0 carregou a configuração principal e os perfis Fiocruz e executou
tarefas mínimas de infraestrutura em SIF e nos sandboxes existentes. A escrita
e o isolamento dos caches reais do Snakemake também foram verificados.
Os 12 testes Python passam com as chamadas `singularity`. A verificação
de sintaxe Bash também passa.

Uma nova execução isolada com os fontes atuais confirmou Nextflow 22.04.0
chamando `singularity exec` em uma imagem SIF e nos dois sandboxes, além da
escrita dos caches reais do Snakemake e da execução do wrapper atual com um
workflow vazio. Construção, downloads e operações de manutenção foram
verificados com mocks/executáveis simulados; os containers existentes não
foram reconstruídos nem atualizados nesta revisão.

Para executar os testes locais de infraestrutura:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -p test_runtime.py -v
bash -n vfnext/containers/add_entries_SnpeffDB.sh tests/pangolin-cache-smoke.sh tests/host-smoke.sh
```

Essas verificações foram realizadas em Ubuntu 26.04.1 AMD64 com Apptainer
1.5.4 e Java 17. A validação do downgrade cobre instalação e infraestrutura;
não equivale à execução completa da pipeline científica. ARM64, PBS,
Ubuntu 20.04 e WSL não foram executados nessa validação. O script
`tests/host-smoke.sh` consta na branch, mas não foi executado nessa etapa.

A documentação da GUI, que pertence a outro repositório, está em
[UBUNTU26.04-COMPATIBILITY.md](https://github.com/camilodotto/ViralFlowGui/blob/v0.53/UBUNTU26.04-COMPATIBILITY.md).

Para reproduzir a comparação dos arquivos versionados com o conteúdo atual
do diretório de trabalho:

```bash
git diff --name-status WallauBioinfo/develop --
git diff WallauBioinfo/develop --
```
