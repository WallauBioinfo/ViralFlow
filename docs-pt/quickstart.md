# Guia Rápido

O ViralFlow fornece 2 modos de uso: **sars-cov2** e **custom**. Independente do modo, o usuário deve passar os caminhos absolutos (todo o caminho até o diretório ou arquivo a ser indicado para o pipeline) para cada arquivo de entrada, caso contrário o pipeline será interrompido durante a execução.

```{important}
Sempre use **caminhos absolutos** ao especificar diretórios e arquivos, ex.: `/home/usuario/dados/` ao invés de `./dados/`
```

## Customizando o Banco do snpEff

Por padrão, a ferramenta snpEff presente no ViralFlow vem configurada apenas com o genoma NC_045512.2 do vírus SARS-CoV-2. Caso você queira incluir a análise do snpEff para outros vírus, você deve atualizar o banco do snpEff com a seguinte linha (exemplo com Dengue):

```bash
viralflow -add_entry_to_snpeff --org_name Dengue --genome_code NC_001474.2
```

## Modo SARS-CoV-2

Neste modelo, a análise é realizada com base no genoma de referência NC_045512.2, e possui, como análises adicionais, a assinatura de linhagens com a ferramenta Pangolin, e a assinatura de clados e mutações com a ferramenta Nextclade.

Para isto, o usuário precisa construir um arquivo com os parâmetros da análise. Um exemplo pode ser [conferido aqui](https://github.com/WallauBioinfo/ViralFlow/blob/main/test_files/sars-cov-2.params).

```bash
viralflow -run --params_file test_files/sars-cov-2.params
```

## Modo Custom

Neste modelo, a análise é realizada com base nos arquivos para o vírus que o usuário deseja analisar. Neste modo, o usuário é responsável por passar cada um dos arquivos necessários para a análise. Caso o usuário deseje realizar a análise do snpEff, deve passar o código do refseq do genoma viral para a análise.

Um exemplo de arquivo de parâmetros pode ser [conferido aqui](https://github.com/WallauBioinfo/ViralFlow/blob/main/test_files/denv.params).

```bash
viralflow -run --params_file test_files/denv.params
```

## Atualização Pangolin

Periodicamente a ferramenta pangolin atualiza o banco de linhagens, bem como a filogenia de classificação do usher, as constelações de mutações do scorpio, e o modelo treinado do pangoLearn. Para realizar a atualização da ferramenta e/ou de suas bases de dados, basta rodar o ViralFlow com algum dos comandos:

```bash
# Atualiza a ferramenta e as bases de dados
viralflow -update_pangolin

# Atualiza apenas as bases de dados
viralflow -update_pangolin_data
```
