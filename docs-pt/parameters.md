# Parâmetros

Esta página descreve todos os argumentos de linha de comando e parâmetros disponíveis no ViralFlow.

## Arquivo de Parâmetros

O ViralFlow requer um arquivo de parâmetros que contém todas as opções de configuração. Exemplos de arquivos de parâmetros podem ser encontrados no [diretório test_files](https://github.com/WallauBioinfo/ViralFlow/tree/main/test_files).

## Referência de Parâmetros

| Argumento | Valor Padrão | Descrição |
|-----------|--------------|-----------|
| `mode` | ILLUMINA | Tecnologia de sequenciamento dos dados de entrada (ILLUMINA ou NANOPORE) |
| `virus` | sars-cov2 | Tipo de análise (sars-cov2 ou custom) |
| `primersBED` | null | Caminho absoluto para o arquivo bed com informações dos primers usados na amplificação genômica (opcional). Fornecê-lo habilita o recorte de primers nos modos ILLUMINA e NANOPORE; sem ele, nenhum recorte é realizado |
| `outDir` | launchDir/output/ | Diretório de saída do ViralFlow onde os resultados e os metadados da execução serão armazenados |
| `samplesheet` | null | Arquivo CSV com as colunas `sample_id`, `fastq_1` e `fastq_2`. Repita um sample ID para fornecer chunks ou lanes ordenados. Não pode ser combinado com `inDir` |
| `inDir` | null | Obsoleto, será removido na v3; use `samplesheet`. Caminho absoluto para o diretório com os dados de entrada (diretório com os arquivos FASTQ). Se não for definido, o ViralFlow usa `launchDir/input/` |
| `runSnpEff` | true | Necessário para executar a ferramenta snpEff (true ou false) |
| `writeMappedReads` | true | Necessário para gerar os arquivos FASTQ contendo as reads de sequenciamento que mapearam no genoma de referência |
| `minLen` | 75 | Tamanho mínimo que as reads devem ter. Reads abaixo deste limite serão eliminadas pelo FastP |
| `depth` | 25 | Profundidade de cobertura mínima para chamar bases consenso. Posições com profundidade de cobertura menor não serão chamadas e um "-" será adicionado à respectiva posição genômica consenso |
| `mapping_quality` | 30 | Limiar de qualidade de mapeamento usado para variant calling |
| `base_quality` | 30 | Limiar de qualidade de base usado para variant calling |
| `minDpIntrahost` | 100 | Profundidade mínima de cobertura por sítio genômico para ser considerado na análise intrahospedeiro |
| `trimLen` | 0 | Número de bases cortadas em ambas as extremidades das reads; 0 desativa. O ILLUMINA aplica no fastp antes do alinhamento; o NANOPORE mascara as bases no BAM alinhado com o bamUtil |
| `refGenomeCode` | null | Código do genoma a ser usado na análise custom |
| `referenceGFF` | null | Arquivo GFF do genoma a ser usado na análise custom |
| `referenceGenome` | null | Arquivo Fasta do genoma a ser usado na análise custom |
| `nextflowSimCalls` | 6 | Número de chamadas simultâneas que o nextflow pode realizar |
| `fastp_threads` | 1 | Número de threads a serem usadas na etapa de filtragem de reads do fastp |
| `bwa_threads` | 1 | Número de threads a serem usadas na etapa de mapeamento do bwa |
| `mafft_threads` | 1 | Número de threads a serem usadas na etapa de alinhamento do mafft |
| `dedup` | false | Este argumento habilita o modo dedup do fastp. Para ativá-lo, mude o valor para true no arquivo de parâmetros de teste |
| `ndedup` | 3 | Quando o modo dedup está ativo, você pode usar níveis de precisão (1 - 6). Você pode alterar este valor, mas recomendamos o padrão. Quanto maior, mais RAM e tempo são consumidos. Para ativar, mude o valor de 1 a 6 no arquivo de parâmetros de teste |

## Parâmetros do NANOPORE

Aplicam-se apenas quando `mode` é `NANOPORE` e são ignorados caso contrário.
`virus`, `refGenomeCode`, `referenceGFF`, `runSnpEff` e as opções de
fastp/bwa/mafft acima pertencem ao modo ILLUMINA; `base_quality` também é
exclusivo do ILLUMINA, enquanto `mapping_quality` é usado nos dois.

| Argumento | Valor Padrão | Descrição |
|-----------|--------------|-----------|
| `base_container` | projectDir/containers/baseContainer.sif | Contêiner que fornece Porechop_ABI, Minimap2, Samtools, BCFtools e bamUtil. Um caminho local `.sif` com Singularity/Apptainer; o perfil `docker` o substitui por uma referência de imagem |
| `clair3_container` | docker://hkubal/clair3@sha256:1430f7b5… | Imagem do Clair3, fixada por digest para que uma reconstrução não possa trocar o chamador de variantes. Corresponde à versão v1.2.0 |
| `clair3_model` | r941_prom_sup_g5014 | Modelo de basecalling usado pelo Clair3, passado como `--model_path`. Deve corresponder a um diretório presente em `/opt/models` dentro da imagem do Clair3 e deve ser compatível com o basecaller e a química que geraram as reads |
| `clair3_qual` | 10 | Qualidade mínima para que o Clair3 reporte uma variante, passada como `--qual` |
| `clair3_chunk_size` | 10000 | Tamanho em bases dos blocos em que o Clair3 divide a referência para chamada paralela, passado como `--chunk_size`. Afeta o tempo de execução e a memória, não os resultados |
| `mapping_quality` | 30 | Qualidade de mapeamento mínima para que uma read seja usada na chamada de variantes, passada ao Clair3 como `--min_mq` |
| `af_threshold` | 0.51 | Limiar de frequência alélica aplicado à saída do Clair3: o BCFtools mantém uma variante quando `FORMAT/AF >= af_threshold`. Nenhuma condição adicional de profundidade ou de `FILTER=PASS` é aplicada. O padrão, acima de 0.5, mantém o alelo majoritário em cada sítio |
| `np_min_depth` | 20 | Limiar de mascaramento do consenso. A cobertura vem de `samtools depth -J -a`, e toda posição cuja profundidade seja **menor ou igual** a este valor é escrita como `N`. No padrão, uma posição precisa de pelo menos 21 reads para ser chamada |
| `porechop_cpus` | 4 | CPUs para a etapa de remoção de adaptadores com o Porechop_ABI |
| `porechop_memory` | 4.GB | Memória para a etapa do Porechop_ABI |
| `minimap_cpus` | 4 | CPUs para a etapa de alinhamento com o Minimap2 |
| `minimap_memory` | 4.GB | Memória para a etapa do Minimap2 |
| `clair3_cpus` | 4 | CPUs para o Clair3, passados como `--threads` |
| `clair3_memory` | 4.GB | Memória para o Clair3 |

`af_threshold` e `np_min_depth` são aplicados de forma independente, de modo que
uma variante de baixa profundidade pode passar pelo filtro de frequência alélica
enquanto essa mesma posição é mascarada no consenso. Cada diretório de amostra
contém um arquivo `<amostra>.nanopore_summary.tsv` que reporta os limiares
configurados, as contagens de variantes e o total de bases mascaradas, de modo
que isso pode ser inspecionado a cada execução.
