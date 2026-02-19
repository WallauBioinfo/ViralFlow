# Documentación de ViralFlow

```{admonition} 🌐 Language / Idioma / Lengua
:class: tip

This documentation is also available in **English**: [Click here](https://viralflow.readthedocs.io/en/latest/)

Esta documentação também está disponível em **Português**: [Clique aqui](https://viralflow.readthedocs.io/pt-br/latest/)
```

**ViralFlow** es un workflow elaborado para profesionales de la salud, con el objetivo de ejecutar todas las etapas de un análisis genómico viral por referencia.

El código fue escrito en el lenguaje de workflow Nextflow, y puede ser aplicado para diferentes virus, donde el usuario, después de instalar la herramienta, necesita rodar solamente una línea de código.

El grupo que desarrolló el ViralFlow no ofrece soporte, por el momento, para errores relacionados con la aplicación del código en otras plataformas. Hasta el presente, el código ha sido utilizado solamente en datos generados por secuenciadores Illumina, empleando estrategia de reads por pares (paired-end) o no emparejados (single-end).

## Visión General del Workflow

![Workflow de ViralFlow](_static/images/flowchart_es.png)

## Análisis Intra-huésped

El ViralFlow presenta un algoritmo propio para detectar regiones de iSNV (intrahost Single Nucleotide Variant).

Para que un determinado locus multi alélico tenga un alelo con menor frecuencia, considerado iSNV, tres condiciones deben cumplirse, estas condiciones están esquematizadas en la parte superior de la figura que se muestra a continuación.

Siguiendo esta lógica, la parte inferior de la figura a continuación representa 3 sitios multi alélicos y ejemplifica la lógica para considerar cuáles de ellos serían considerados iSNVs.

![Lógica de Detección de iSNV](_static/images/viralflow_snv_es.png)

## Contenido de la Documentación

```{toctree}
:maxdepth: 2

installation
dependencies
quickstart
parameters
outputs
```

## Publicaciones

- **Versión 0.1**: Dezordi, F. Z., et al. (2022). ViralFlow: A Versatile Automated Workflow for SARS-CoV-2 Genome Assembly, Lineage Assignment, Mutations and Intrahost Variant Detection. *Viruses*, 14(2), 217. [https://www.mdpi.com/1999-4915/14/2/217](https://www.mdpi.com/1999-4915/14/2/217)

- **Versión 1.0**: da Silva, A. F., et al. (2024). ViralFlow v1.0: characterization and annotation of viral genomes. *NAR Genomics and Bioinformatics*, 6(2), lqae056. [https://academic.oup.com/nargab/article/6/2/lqae056/7682253](https://academic.oup.com/nargab/article/6/2/lqae056/7682253)

## Enlaces Rápidos

- [Repositorio GitHub](https://github.com/WallauBioinfo/ViralFlow)
- [Versiones de Herramientas](https://github.com/WallauBioinfo/ViralFlow/tree/main/versions)
- [Issues/Bugs](https://github.com/WallauBioinfo/ViralFlow/issues)
- **Licencia:** MIT
