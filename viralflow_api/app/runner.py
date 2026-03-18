"""
Runner assíncrono: constrói o comando Nextflow, executa, coleta resultados
e dispara o callback para o MacWorp ao término.

Os resultados são escritos DIRETAMENTE em job.out_dir (pasta do projeto
no NFS compartilhado), sem cópia — o MacWorp lê do mesmo caminho.
"""

import os
import httpx
import asyncio
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Optional

from app.models import JobStatus, JobResult, OutputFile, PipelineJob
from app.job_store import job_store
from app.config import settings


def _build_nextflow_cmd(job: PipelineJob, work_dir: Path) -> list[str]:
    """Monta o comando nextflow run com todos os parâmetros do job."""

    # out_dir é a pasta do projeto MacWorp no NFS — os resultados
    # vão diretamente para lá sem precisar de cópia posterior
    out_dir = job.out_dir if job.out_dir else str(work_dir / "results")

    cmd = [
        settings.NEXTFLOW_BIN, "run", settings.VIRALFLOW_MAIN,
        "-w",        str(work_dir / "work"),
        "--inDir",            job.in_dir,
        "--outDir",           out_dir,
        "--virus",            job.virus,
        "--depth",            str(job.depth),
        "--minLen",           str(job.min_len),
        "--minDpIntrahost",   str(job.min_dp_intrahost),
    ]

    if job.run_snpeff:
        cmd.append("--runSnpEff")
    if job.write_mapped_reads:
        cmd.append("--writeMappedReads")

    if job.extra_params:
        for key, value in job.extra_params.items():
            cmd.extend([f"--{key}", str(value)])

    return cmd, out_dir


def _collect_output_files(results_dir: Path) -> list[OutputFile]:
    files = []
    if not results_dir.exists():
        return files
    for f in results_dir.rglob("*"):
        if f.is_file():
            files.append(OutputFile(
                name=f.name,
                path=str(f),
                size_bytes=f.stat().st_size,
                description=_describe_file(f.suffix),
            ))
    return files


def _describe_file(suffix: str) -> Optional[str]:
    descriptions = {
        ".fasta": "Sequência consenso",
        ".fa":    "Sequência consenso",
        ".vcf":   "Variantes chamadas (VCF)",
        ".bam":   "Alinhamento (BAM)",
        ".html":  "Relatório HTML",
        ".csv":   "Tabela de métricas",
        ".tsv":   "Tabela de métricas",
        ".json":  "Métricas JSON",
        ".log":   "Log de execução",
        ".png":   "Gráfico de cobertura",
        ".pdf":   "Relatório PDF",
    }
    return descriptions.get(suffix.lower())


async def _send_callback(job_id: str, job: PipelineJob, result: JobResult):
    if not job.callback_url:
        return
    headers = {"Content-Type": "application/json"}
    if job.callback_api_key:
        headers["X-Api-Key"] = job.callback_api_key
    try:
        async with httpx.AsyncClient(timeout=30) as client:
            response = await client.post(
                job.callback_url,
                json=result.model_dump(mode="json"),
                headers=headers,
            )
            response.raise_for_status()
            print(f"[callback] MacWorp notificado para job {job_id}: HTTP {response.status_code}")
    except Exception as exc:
        print(f"[callback] Erro ao notificar MacWorp para job {job_id}: {exc}")


async def run_viralflow_pipeline(job_id: str, job: PipelineJob):
    """
    Executa o pipeline ViralFlow em background.
    Grava resultados diretamente em job.out_dir (pasta do projeto no NFS).
    """
    started_at = datetime.utcnow()

    # Diretório de trabalho temporário do Nextflow (work/ e .nextflow/)
    work_dir = Path(settings.JOBS_BASE_DIR) / job_id
    work_dir.mkdir(parents=True, exist_ok=True)

    # out_dir = pasta do projeto MacWorp no NFS
    cmd, out_dir = _build_nextflow_cmd(job, work_dir)
    out_dir_path = Path(out_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)

    job_store.update_status(
        job_id,
        status=JobStatus.RUNNING,
        started_at=started_at,
        progress_message=f"Iniciando ViralFlow → resultados em {out_dir}",
    )

    log_path = work_dir / "nextflow.log"

    env = os.environ.copy()
    #if settings.CONDA_ENV:
    #    cmd = ["micromamba", "run", "-n", settings.CONDA_ENV] + cmd

    exit_code = -1
    error_message = None

    try:
        with open(log_path, "w") as log_file:
            proc = subprocess.Popen(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                cwd=str(work_dir),
                env=env,
            )
            job_store.register_process(job_id, proc)

            timeout = settings.JOB_TIMEOUT if settings.JOB_TIMEOUT > 0 else None
            loop = asyncio.get_event_loop()
            exit_code = await loop.run_in_executor(
                None, lambda: proc.wait(timeout=timeout)
            )

    except subprocess.TimeoutExpired:
        proc.kill()
        exit_code = -1
        error_message = "Timeout: pipeline excedeu o tempo máximo."
    except FileNotFoundError as exc:
        exit_code = -1
        error_message = f"Executável não encontrado: {exc}"
    except Exception as exc:
        exit_code = -1
        error_message = f"Erro inesperado: {exc}"

    finished_at = datetime.utcnow()
    duration = (finished_at - started_at).total_seconds()

    # Log das últimas 100 linhas
    nextflow_log_tail = None
    if log_path.exists():
        try:
            lines = log_path.read_text(errors="replace").splitlines()
            nextflow_log_tail = "\n".join(lines[-100:])
        except Exception:
            pass

    # Coleta arquivos gerados na pasta do projeto MacWorp
    output_files = _collect_output_files(out_dir_path)

    status = JobStatus.COMPLETED if exit_code == 0 else JobStatus.FAILED

    result = JobResult(
        job_id=job_id,
        macworp_project_id=job.macworp_project_id,
        macworp_run_id=job.macworp_run_id,
        status=status,
        started_at=started_at,
        finished_at=finished_at,
        duration_seconds=duration,
        exit_code=exit_code,
        output_dir=out_dir,       # caminho NFS direto
        output_files=output_files,
        nextflow_log=nextflow_log_tail,
        error_message=error_message,
    )

    job_store.update_status(
        job_id,
        status=status,
        finished_at=finished_at,
        result=result,
        progress_message="Concluído." if status == JobStatus.COMPLETED
                         else f"Falhou: {error_message}",
    )

    await _send_callback(job_id, job, result)