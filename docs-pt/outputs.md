# Arquivos de Saída

A ferramenta ViralFlow gera 2 tipos de saída: **específica** e **compilada**.

## Saídas Específicas

Para cada amostra na análise, um diretório de resultados será criado com o padrão `prefixo_results` onde o prefixo é um código criado a partir do nome do arquivo fastq da amostra. Cada diretório possui os seguintes resultados:

| Arquivo | Descrição |
|---------|-----------|
| `prefixo.all.fa.pango.out.csv` | Arquivo tabular com os resultados da ferramenta pangolin |
| `prefixo.ann.vcf` | Arquivo em formato TSV (Tab separated value) com anotação das variantes feitas pela ferramenta snpEff |
| `prefixo.depth5.all.fa.nextclade.csv` | Arquivo tabular com os resultados da ferramenta nextclade |
| `prefixo.depth5.amb.fa` | Genoma consenso com nucleotídeos ambíguos em posições de múltiplos alelos |
| `prefixo.depth5.fa` | Genoma consenso com nucleotídeos majoritários. Consenso majoritário, normalmente depositado no GISAID |
| `prefixo.depth5.fa.algn` | Genoma consenso alinhado com genoma de referência (seguindo mesmo tamanho, mafft --keeplenght) |
| `prefixo.depth5.fa.algn.minor.fa` | Genoma consenso com nucleotídeos minoritários |
| `prefixo.depth5.fa.bc.intrahost.short.tsv` | Arquivo tabular resumindo as posições genômicas onde variantes intrahospedeiro são suportadas |
| `prefixo.depth5.fa.bc.intrahost.tsv` | Arquivo tabular com todas as informações sobre posições de variantes intrahospedeiro |
| `prefixo.mapped.R1.fq.gz` | Arquivo FASTQ R1 com reads mapeadas |
| `prefixo.mapped.R2.fq.gz` | Arquivo FASTQ R2 com reads mapeadas |
| `prefixo.metrics.genome.tsv` | Arquivo tabular com métricas de profundidade de mapeamento e cobertura |
| `prefixo.fastp.html` | Arquivo HTML resumindo os resultados da ferramenta fastp |
| `prefixo_snpEff_summary.html` | Arquivo HTML resumindo os resultados da ferramenta snpEff |
| `prefixo.sorted.bam` | Arquivo com o mapeamento ordenado das reads da amostra contra o genoma de referência |
| `prefixo.tsv` | Arquivo similar ao VCF gerado pelo iVar |
| `prefixo.unmapped.R1.bam.fq` | Arquivo FASTQ R1 com reads não mapeadas |
| `prefixo.unmapped.R2.bam.fq` | Arquivo FASTQ R2 com reads não mapeadas |
| `prefixo.vcf` | Arquivo VCF |
| `metrics.alignment_summary_metrics` | Arquivo de texto com um resumo de várias métricas do mapeamento |
| `nextclade.errors.csv` | Arquivo relatando erros do nextclade |
| `nextclade_gene_*.translation.fasta` | Arquivo com as proteínas de cada gene |
| `snpEff_genes.txt` | Arquivo tabular com o número de variantes e impacto estimado por gene |
| `wgs` | Arquivo textual com métricas de mapeamento, incluindo profundidade por região |
| `prefixo_coveragePlot.png` | Arquivo PNG com visualização gráfica da cobertura do genoma |
| `prefixo_coveragePlot.svg` | Arquivo SVG com visualização gráfica da cobertura do genoma |
| `prefixo_coveragePlot.html` | Arquivo HTML com visualização gráfica da cobertura do genoma |
| `prefixo_snpPlot.png` | Arquivo PNG com visualização gráfica dos SNPs |
| `prefixo_snpPlot.svg` | Arquivo SVG com visualização gráfica dos SNPs |
| `prefixo_snpPlot.html` | Arquivo HTML com visualização gráfica dos SNPs |
| `prefixo_depthPlot.png` | Arquivo PNG com visualização gráfica da profundidade |
| `prefixo_depthPlot.svg` | Arquivo SVG com visualização gráfica da profundidade |
| `prefixo_depthPlot.html` | Arquivo HTML com visualização gráfica da profundidade |

## Saídas Compiladas

O ViralFlow também gera saídas compiladas que agregam resultados de todas as amostras:

| Arquivo | Descrição |
|---------|-----------|
| `assembly_statistics_summary.tsv` | Resumo das estatísticas de montagem para todas as amostras |
| `alignment_metrics_summary.tsv` | Resumo das métricas de alinhamento para todas as amostras |
| `all_consensus.fasta` | Sequências consenso concatenadas para todas as amostras |
| `intrahost_summary.tsv` | Informações compiladas de variantes intrahospedeiro para todas as amostras |
