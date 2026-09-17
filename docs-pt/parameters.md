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
