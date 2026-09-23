# Parámetros

Esta página describe todos los argumentos de línea de comando y parámetros disponibles en ViralFlow.

## Archivo de Parámetros

ViralFlow requiere un archivo de parámetros que contiene todas las opciones de configuración. Ejemplos de archivos de parámetros se pueden encontrar en el [directorio test_files](https://github.com/WallauBioinfo/ViralFlow/tree/main/test_files).

## Referencia de Parámetros

| Argumento | Valor Predeterminado | Descripción |
|-----------|----------------------|-------------|
| `mode` | ILLUMINA | Tecnología de secuenciación de los datos de entrada (ILLUMINA o NANOPORE) |
| `virus` | sars-cov2 | Tipo de análisis (sars-cov2 o custom) |
| `primersBED` | null | Ruta absoluta al archivo bed con información de primers usados en la amplificación genómica (opcional). Proporcionarlo habilita el recorte de primers tanto en modo ILLUMINA como NANOPORE; sin él no se realiza ningún recorte |
| `outDir` | launchDir/output/ | Directorio de salida de ViralFlow donde se almacenarán los resultados y los metadatos de ejecución |
| `samplesheet` | null | Archivo CSV con las columnas `sample_id`, `fastq_1` y `fastq_2`. Repita un sample ID para proporcionar chunks o lanes ordenados. No puede combinarse con `inDir` |
| `inDir` | null | Obsoleto, será eliminado en v3; use `samplesheet`. Ruta absoluta al directorio con los datos de entrada (directorio con los archivos FASTQ). Si no se define, ViralFlow usa `launchDir/input/` |
| `runSnpEff` | true | Necesario para ejecutar la herramienta snpEff (true o false) |
| `writeMappedReads` | true | Necesario para generar los archivos FASTQ que contienen las reads de secuenciación que mapearon al genoma de referencia |
| `minLen` | 75 | Tamaño mínimo que las reads deben tener. Reads por debajo de este umbral serán eliminadas por FastP |
| `depth` | 25 | Profundidad de cobertura mínima para llamar bases de consenso. Posiciones con profundidad de cobertura menor no serán llamadas y se agregará un "-" a la respectiva posición genómica de consenso |
| `mapping_quality` | 30 | Umbral de calidad de mapeo usado para variant calling |
| `base_quality` | 30 | Umbral de calidad de base usado para variant calling |
| `minDpIntrahost` | 100 | Profundidad de cobertura mínima por sitio genómico para ser considerado en el análisis intra-huésped |
| `trimLen` | 0 | Número de bases recortadas en ambos extremos de las reads; 0 lo desactiva. ILLUMINA lo aplica en fastp antes del alineamiento; NANOPORE enmascara las bases en el BAM alineado con bamUtil |
| `refGenomeCode` | null | Código del genoma a ser usado en el análisis custom |
| `referenceGFF` | null | Archivo GFF del genoma a ser usado en el análisis custom |
| `referenceGenome` | null | Archivo Fasta del genoma a ser usado en el análisis custom |
| `nextflowSimCalls` | 6 | Número de llamadas simultáneas que nextflow puede realizar |
| `fastp_threads` | 1 | Número de threads a ser usados en la etapa de filtrado de reads de fastp |
| `bwa_threads` | 1 | Número de threads a ser usados en la etapa de mapeo de bwa |
| `mafft_threads` | 1 | Número de threads a ser usados en la etapa de alineamiento de mafft |
| `dedup` | false | Este argumento habilita el modo dedup de fastp. Para activarlo, cambie el valor a true en el archivo de parámetros de prueba |
| `ndedup` | 3 | Cuando el modo dedup está activo, puede usar niveles de precisión (1 - 6). Puede cambiar este valor, pero recomendamos el estándar. Cuanto mayor, más RAM y tiempo se consumen. Para activar, cambie el valor de 1 a 6 en el archivo de parámetros de prueba |

## Parámetros de NANOPORE

Se aplican únicamente cuando `mode` es `NANOPORE` y se ignoran en caso
contrario. `virus`, `refGenomeCode`, `referenceGFF`, `runSnpEff` y las opciones
de fastp/bwa/mafft anteriores pertenecen al modo ILLUMINA; `base_quality` es
también exclusivo de ILLUMINA, mientras que `mapping_quality` se usa en ambos.

| Argumento | Valor Predeterminado | Descripción |
|-----------|----------------------|-------------|
| `base_container` | projectDir/containers/baseContainer.sif | Contenedor que provee Porechop_ABI, Minimap2, Samtools, BCFtools y bamUtil. Una ruta local `.sif` con Singularity/Apptainer; el perfil `docker` lo reemplaza por una referencia de imagen |
| `clair3_container` | docker://hkubal/clair3@sha256:1430f7b5… | Imagen de Clair3, fijada por digest para que una reconstrucción no pueda cambiar el llamador de variantes. Corresponde a la versión v1.2.0 |
| `clair3_model` | r941_prom_sup_g5014 | Modelo de basecalling que usa Clair3, pasado como `--model_path`. Debe corresponder a un directorio presente en `/opt/models` dentro de la imagen de Clair3, y debe coincidir con el basecaller y la química que generaron las reads |
| `clair3_qual` | 10 | Calidad mínima para que Clair3 reporte una variante, pasada como `--qual` |
| `clair3_chunk_size` | 10000 | Tamaño en bases de los fragmentos en que Clair3 divide la referencia para el llamado paralelo, pasado como `--chunk_size`. Afecta el tiempo de ejecución y la memoria, no los resultados |
| `mapping_quality` | 30 | Calidad de mapeo mínima para que una read se use en el llamado de variantes, pasada a Clair3 como `--min_mq` |
| `af_threshold` | 0.51 | Umbral de frecuencia alélica aplicado a la salida de Clair3: BCFtools conserva una variante cuando `FORMAT/AF >= af_threshold`. No se aplica ninguna condición adicional de profundidad ni de `FILTER=PASS`. El valor por defecto, superior a 0.5, conserva el alelo mayoritario en cada sitio |
| `np_min_depth` | 20 | Umbral de enmascaramiento del consenso. La cobertura proviene de `samtools depth -J -aa`, y toda posición cuya profundidad sea **menor o igual** a este valor se escribe como `N`. Con el valor por defecto, una posición necesita al menos 21 reads para ser llamada |
| `porechop_cpus` | 4 | CPUs para la etapa de remoción de adaptadores con Porechop_ABI |
| `porechop_memory` | 4.GB | Memoria para la etapa de Porechop_ABI |
| `minimap_cpus` | 4 | CPUs para la etapa de alineamiento con Minimap2 |
| `minimap_memory` | 4.GB | Memoria para la etapa de Minimap2 |
| `clair3_cpus` | 4 | CPUs para Clair3, pasados como `--threads` |
| `clair3_memory` | 4.GB | Memoria para Clair3 |

`af_threshold` y `np_min_depth` se aplican de forma independiente, por lo que una
variante de baja profundidad puede superar el filtro de frecuencia alélica
mientras esa misma posición queda enmascarada en el consenso. Cada directorio de
muestra contiene un archivo `<muestra>.nanopore_summary.tsv` que reporta los umbrales
configurados, los conteos de variantes y el total de bases enmascaradas, de modo
que esto puede inspeccionarse en cada ejecución. Tenga en cuenta que
`masked_bases` y `consensus_n_bases` cuentan cosas distintas; consulte
`NANOPORE.md` para saber por qué difieren con una referencia que contiene códigos
de ambigüedad.
