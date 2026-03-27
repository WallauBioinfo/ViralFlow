"""Configurações da ViralFlow API — carregadas via variáveis de ambiente."""

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    # Chave secreta compartilhada com o MacWorp
    API_KEY: str = "change-me-in-production"

    # Caminho para o executável do Nextflow
    NEXTFLOW_BIN: str = "nextflow"

    # Caminho para o script principal do ViralFlow
    VIRALFLOW_MAIN: str = "/opt/viralflow/vfnext/main.nf"

    # Diretório base para outputs dos jobs
    JOBS_BASE_DIR: str = "/data/viralflow_jobs"

    # Ambiente conda/micromamba a ativar antes de rodar Nextflow
    CONDA_ENV: str = "viralflow"

    NEXTFLOW_PROFILE: str = "micromamba"

    # Timeout máximo de execução por job (segundos) — 0 = sem limite
    JOB_TIMEOUT: int = 0

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


settings = Settings()
