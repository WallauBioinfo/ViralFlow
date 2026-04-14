# Archivos de Salida

La herramienta ViralFlow genera 2 tipos de salida: **específica** y **compilada**.

## Salidas Específicas

Para cada muestra en el análisis, se creará un directorio de resultados con el patrón `prefijo_results` donde el prefijo es un código creado a partir del nombre del archivo fastq de la muestra. Cada directorio tiene los siguientes resultados:

| Archivo | Descripción |
|---------|-------------|
| `prefijo.all.fa.pango.out.csv` | Archivo tabular con los resultados de la herramienta pangolin |
| `prefijo.ann.vcf` | Archivo en formato TSV (Tab separated value) con anotación de las variantes hechas por la herramienta snpEff |
| `prefijo.depth5.all.fa.nextclade.csv` | Archivo tabular con los resultados de la herramienta nextclade |
| `prefijo.depth5.amb.fa` | Genoma consenso con nucleótidos ambiguos en posiciones de múltiples alelos |
| `prefijo.depth5.fa` | Genoma consenso con nucleótidos mayoritarios. Consenso mayoritario, normalmente depositado en GISAID |
| `prefijo.depth5.fa.algn` | Genoma consenso alineado con genoma de referencia (siguiendo el mismo tamaño, mafft --keeplenght) |
| `prefijo.depth5.fa.algn.minor.fa` | Genoma consenso con nucleótidos minoritarios |
| `prefijo.depth5.fa.bc.intrahost.short.tsv` | Archivo tabular resumiendo las posiciones genómicas donde se soportan variantes intra-huésped |
| `prefijo.depth5.fa.bc.intrahost.tsv` | Archivo tabular con toda la información sobre posiciones de variantes intra-huésped |
| `prefijo.mapped.R1.fq.gz` | Archivo FASTQ R1 con reads mapeadas |
| `prefijo.mapped.R2.fq.gz` | Archivo FASTQ R2 con reads mapeadas |
| `prefijo.metrics.genome.tsv` | Archivo tabular con métricas de profundidad de mapeo y cobertura |
| `prefijo.fastp.html` | Archivo HTML resumiendo los resultados de la herramienta fastp |
| `prefijo_snpEff_summary.html` | Archivo HTML resumiendo los resultados de la herramienta snpEff |
| `prefijo.sorted.bam` | Archivo con el mapeo ordenado de las reads de la muestra contra el genoma de referencia |
| `prefijo.tsv` | Archivo similar al VCF generado por iVar |
| `prefijo.unmapped.R1.bam.fq` | Archivo FASTQ R1 con reads no mapeadas |
| `prefijo.unmapped.R2.bam.fq` | Archivo FASTQ R2 con reads no mapeadas |
| `prefijo.vcf` | Archivo VCF |
| `metrics.alignment_summary_metrics` | Archivo de texto con un resumen de varias métricas del mapeo |
| `nextclade.errors.csv` | Archivo informando errores de nextclade |
| `nextclade_gene_*.translation.fasta` | Archivo con las proteínas de cada gen |
| `snpEff_genes.txt` | Archivo tabular con el número de variantes e impacto estimado por gen |
| `wgs` | Archivo textual con métricas de mapeo, incluyendo profundidad por región |
| `prefijo_coveragePlot.png` | Archivo PNG con visualización gráfica de la cobertura del genoma |
| `prefijo_coveragePlot.svg` | Archivo SVG con visualización gráfica de la cobertura del genoma |
| `prefijo_coveragePlot.html` | Archivo HTML con visualización gráfica de la cobertura del genoma |
| `prefijo_snpPlot.png` | Archivo PNG con visualización gráfica de los SNPs |
| `prefijo_snpPlot.svg` | Archivo SVG con visualización gráfica de los SNPs |
| `prefijo_snpPlot.html` | Archivo HTML con visualización gráfica de los SNPs |
| `prefijo_depthPlot.png` | Archivo PNG con visualización gráfica de la profundidad |
| `prefijo_depthPlot.svg` | Archivo SVG con visualización gráfica de la profundidad |
| `prefijo_depthPlot.html` | Archivo HTML con visualización gráfica de la profundidad |

## Salidas Compiladas

ViralFlow también genera salidas compiladas que agregan resultados de todas las muestras:

| Archivo | Descripción |
|---------|-------------|
| `assembly_statistics_summary.tsv` | Resumen de las estadísticas de ensamblaje para todas las muestras |
| `alignment_metrics_summary.tsv` | Resumen de las métricas de alineamiento para todas las muestras |
| `all_consensus.fasta` | Secuencias consenso concatenadas para todas las muestras |
| `intrahost_summary.tsv` | Información compilada de variantes intra-huésped para todas las muestras |
