# Parámetros

Esta página describe todos los argumentos de línea de comando y parámetros disponibles en ViralFlow.

## Archivo de Parámetros

ViralFlow requiere un archivo de parámetros que contiene todas las opciones de configuración. Ejemplos de archivos de parámetros se pueden encontrar en el [directorio test_files](https://github.com/WallauBioinfo/ViralFlow/tree/main/test_files).

## Referencia de Parámetros

| Argumento | Valor Predeterminado | Descripción |
|-----------|----------------------|-------------|
| `virus` | sars-cov2 | Tipo de análisis (sars-cov2 o custom) |
| `primersBED` | null | Ruta absoluta al archivo bed con información de primers usados en la amplificación genómica (opcional) |
| `outDir` | launchDir/output/ | Ruta absoluta al directorio donde se almacenarán los resultados |
| `inDir` | launchDir/input/ | Ruta absoluta al directorio con los datos de entrada (directorio con los archivos FASTQ) |
| `runSnpEff` | true | Necesario para ejecutar la herramienta snpEff (true o false) |
| `writeMappedReads` | true | Necesario para generar los archivos FASTQ que contienen las reads de secuenciación que mapearon al genoma de referencia |
| `minLen` | 75 | Tamaño mínimo que las reads deben tener. Reads por debajo de este umbral serán eliminadas por FastP |
| `depth` | 5 | Profundidad de cobertura mínima para llamar bases de consenso. Posiciones con profundidad de cobertura menor no serán llamadas y se agregará un "-" a la respectiva posición genómica de consenso |
| `mapping_quality` | 30 | Umbral de calidad de mapeo usado para variant calling |
| `base_quality` | 30 | Umbral de calidad de base usado para variant calling |
| `minDpIntrahost` | 100 | Profundidad de cobertura mínima por sitio genómico para ser considerado en el análisis intra-huésped |
| `trimLen` | 0 | Número de bases que deben ser recortadas en ambos extremos de las reads |
| `refGenomeCode` | null | Código del genoma a ser usado en el análisis custom |
| `referenceGFF` | null | Archivo GFF del genoma a ser usado en el análisis custom |
| `referenceGenome` | null | Archivo Fasta del genoma a ser usado en el análisis custom |
| `nextflowSimCalls` | null | Número de llamadas simultáneas que nextflow puede realizar |
| `fastp_threads` | 1 | Número de threads a ser usados en la etapa de filtrado de reads de fastp |
| `bwa_threads` | 1 | Número de threads a ser usados en la etapa de mapeo de bwa |
| `mafft_threads` | 1 | Número de threads a ser usados en la etapa de alineamiento de mafft |
| `dedup` | false | Este argumento habilita el modo dedup de fastp. Para activarlo, cambie el valor a true en el archivo de parámetros de prueba |
| `ndedup` | 3 | Cuando el modo dedup está activo, puede usar niveles de precisión (1 - 6). Puede cambiar este valor, pero recomendamos el estándar. Cuanto mayor, más RAM y tiempo se consumen. Para activar, cambie el valor de 1 a 6 en el archivo de parámetros de prueba |
