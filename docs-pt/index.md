# Documentação do ViralFlow

```{admonition} 🌐 Language / Idioma / Lengua
:class: tip

This documentation is also available in **English**: [Click here](https://viralflow.readthedocs.io/en/latest/)

Esta documentación también está disponible en **Español**: [Haga clic aquí](https://viralflow.readthedocs.io/es/latest/)
```

**ViralFlow** é um workflow feito para profissionais da saúde, cujo objetivo é realizar todas as etapas de uma análise genômica viral por referência.

O código foi escrito na linguagem de workflow Nextflow, e pode ser aplicado para diferentes vírus, onde o usuário, após instalar a ferramenta, precisa rodar apenas 1 linha de código.

A equipe de desenvolvedores do ViralFlow não oferece suporte, até o presente momento, para erros relacionados à aplicação do código para outras plataformas. Atualmente, o código foi testado apenas para dados gerados por sequenciadores Illumina, utilizando estratégia de reads pareados (paired-end) ou não pareados (single-end).

## Visão Geral do Workflow

![Workflow do ViralFlow](_static/images/flowchart_pt-br.png)

## Análise Intrahospedeiro

O ViralFlow apresenta um algoritmo próprio para detectar regiões de iSNV (intrahost Single Nucleotide Variant).

Para que um determinado loci multi alélico tenha um alelo em menor frequência considerado como iSNV, três condições precisam ser satisfeitas, estas condições estão esquematizadas na parte superior da figura abaixo.

Nesta lógica, a parte inferior da figura abaixo representa 3 sítios multi alélicos e exemplifica a lógica para considerar quais deles seriam considerados como iSNVs.

![Lógica de Detecção de iSNV](_static/images/viralflow_snv_pt-br.png)

## Conteúdo da Documentação

```{toctree}
:maxdepth: 2

installation
dependencies
quickstart
parameters
outputs
```

## Publicações

- **Versão 0.1**: Dezordi, F. Z., et al. (2022). ViralFlow: A Versatile Automated Workflow for SARS-CoV-2 Genome Assembly, Lineage Assignment, Mutations and Intrahost Variant Detection. *Viruses*, 14(2), 217. [https://www.mdpi.com/1999-4915/14/2/217](https://www.mdpi.com/1999-4915/14/2/217)

- **Versão 1.0**: da Silva, A. F., et al. (2024). ViralFlow v1.0: characterization and annotation of viral genomes. *NAR Genomics and Bioinformatics*, 6(2), lqae056. [https://academic.oup.com/nargab/article/6/2/lqae056/7682253](https://academic.oup.com/nargab/article/6/2/lqae056/7682253)

## Links Rápidos

- [Repositório GitHub](https://github.com/WallauBioinfo/ViralFlow)
- [Versões das Ferramentas](https://github.com/WallauBioinfo/ViralFlow/tree/main/versions)
- [Issues/Bugs](https://github.com/WallauBioinfo/ViralFlow/issues)
- **Licença:** Não especificada
