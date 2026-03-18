# Guía Rápida

El ViralFlow dispone de dos modos de uso: **sars-cov2** y **custom**. Independientemente del modo de uso, el usuario debe realizar los caminos absolutos (todo el camino hasta el directorio o archivo que será indicado para el pipeline) para cada archivo de entrada, en el caso de no hacerse, el pipeline será interrumpido durante la ejecución.

```{important}
Siempre use **rutas absolutas** al especificar directorios y archivos, ej.: `/home/usuario/datos/` en lugar de `./datos/`
```

## Customizando el Banco de snpEff

La herramienta snpEff presente en el ViralFlow ya viene configurada solamente con el genoma NC_045512.2 del virus SARS-CoV-2. En el caso de querer el análisis del snpEff para otros virus, deberá actualizarse el  banco del snpEff con la siguiente línea de comando (ejemplo Dengue):

```bash
viralflow -add_entry_to_snpeff --org_name Dengue --genome_code NC_001474.2
```

## Modo SARS-CoV-2

En este modelo, el análisis es realizado considerando el genoma de referencia NC_045512.2, y posee como análisis adicionales la inclusión de linajes con la herramienta Pangolin, y la inclusión de clados y mutaciones con la herramienta Nextclade.

Para ello, el usuario necesita construir un archivo con los parámetros del análisis. Un ejemplo puede ser [verificado aquí](https://github.com/WallauBioinfo/ViralFlow/blob/main/test_files/sars-cov-2.params).

```bash
viralflow -run --params_file test_files/sars-cov-2.params
```

## Modo Custom

En este modelo, el análisis es realizado considerando los archivos para el virus que el usuario desea analizar. De esta forma, el usuario es responsable por pasar cada uno de los archivos necesarios para el análisis. En el caso de que el usuario desee realizar el análisis del snpEff, debe pasar el código del refseq del genoma viral para el análisis.

Un ejemplo de archivo de parámetros puede ser [verificado aquí](https://github.com/WallauBioinfo/ViralFlow/blob/main/test_files/denv.params).

```bash
viralflow -run --params_file test_files/denv.params
```

## Actualización Pangolin

Periódicamente, la herramienta Pangolin actualiza el banco de linajes, así como la filogenia de la clasificación de Usher, las constelaciones de mutaciones de scorpio y el modelo utilizado del pangoLearn. Para realizar la actualización de la herramienta y/o de sus bases de datos, basta ejecutar el ViralFlow con alguno de los comandos:

```bash
# Actualiza la herramienta y las bases de datos
viralflow -update_pangolin

# Solo actualiza las bases de datos
viralflow -update_pangolin_data
```
