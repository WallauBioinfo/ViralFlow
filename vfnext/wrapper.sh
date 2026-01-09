#!/bin/bash
# Wrapper script para ativar o ambiente viralflow e executar o Nextflow
# Salve como: /home/pinho/ViralFlow/vfnext/wrapper.sh
# Dê permissão de execução: chmod +x /home/pinho/ViralFlow/vfnext/wrapper.sh

set -e  # Sair se houver erro

# Inicializar micromamba (ajuste o caminho conforme sua instalação)
eval "$(micromamba shell hook --shell bash)"

# Ativar o ambiente viralflow
micromamba activate viralflow

# Executar o nextflow com todos os argumentos passados
exec nextflow "$@"
