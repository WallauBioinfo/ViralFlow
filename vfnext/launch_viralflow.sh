#!/bin/bash
# Script launcher completo para ViralFlow via Macworp
# Salve como: /home/pinho/ViralFlow/vfnext/launch_viralflow.sh
# Dê permissão: chmod +x /home/pinho/ViralFlow/vfnext/launch_viralflow.sh

set -e

# Cores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== ViralFlow Launcher ===${NC}"

# 1. Verificar se micromamba está disponível
if ! command -v micromamba &> /dev/null; then
    echo -e "${RED}ERRO: micromamba não encontrado no PATH${NC}"
    echo "Tentando inicializar micromamba..."
    
    # Tentar inicializar micromamba (ajuste o caminho se necessário)
    if [ -f "${HOME}/.bashrc" ]; then
        source "${HOME}/.bashrc"
    fi
    
    # Tentar hook do micromamba
    eval "$(micromamba shell hook --shell bash 2>/dev/null || true)"
fi

# 2. Verificar novamente
if ! command -v micromamba &> /dev/null; then
    echo -e "${RED}ERRO: Não foi possível inicializar micromamba${NC}"
    exit 1
fi

echo -e "${YELLOW}Micromamba encontrado: $(which micromamba)${NC}"

# 3. Ativar ambiente viralflow
echo -e "${YELLOW}Ativando ambiente viralflow...${NC}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate viralflow

if [ $? -ne 0 ]; then
    echo -e "${RED}ERRO: Não foi possível ativar o ambiente viralflow${NC}"
    exit 1
fi

echo -e "${GREEN}Ambiente ativado com sucesso!${NC}"

# 4. Verificar se nextflow está disponível
if ! command -v nextflow &> /dev/null; then
    echo -e "${RED}ERRO: nextflow não encontrado no ambiente viralflow${NC}"
    exit 1
fi

echo -e "${YELLOW}Nextflow encontrado: $(which nextflow)${NC}"
echo -e "${YELLOW}Versão: $(nextflow -version 2>&1 | head -n 1)${NC}"

# 5. Executar nextflow
echo -e "${GREEN}Executando ViralFlow...${NC}"
echo -e "${YELLOW}Comando: nextflow $@${NC}"
echo ""

# Executar o nextflow passando todos os argumentos
exec nextflow "$@"
