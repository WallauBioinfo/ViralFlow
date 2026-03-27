"""
ViralFlow API — servidor FastAPI dedicado para execução do pipeline ViralFlow.
Recebe jobs do MacWorp, executa o pipeline Nextflow e retorna os resultados.
"""

import uuid
import asyncio
import subprocess
from pathlib import Path
from datetime import datetime
from typing import Optional

from fastapi import FastAPI, HTTPException, BackgroundTasks, Header
from fastapi.responses import JSONResponse

from app.models import (
    PipelineJob, PipelineJobResponse, JobStatus,
    JobResult, JobStatusResponse
)
from app.job_store import job_store
from app.config import settings
from app.runner import run_viralflow_pipeline

app = FastAPI(
    title="ViralFlow API",
    description="API REST para execução remota do pipeline ViralFlow pelo MacWorp.",
    version="1.0.0",
)


# ──────────────────────────────────────────
# Autenticação simples por API Key
# ──────────────────────────────────────────

def verify_api_key(x_api_key: str = Header(...)):
    if x_api_key != settings.API_KEY:
        raise HTTPException(status_code=401, detail="API Key inválida.")
    return x_api_key


# ──────────────────────────────────────────
# Rotas
# ──────────────────────────────────────────

@app.get("/health")
def health_check():
    """Verifica se o servidor ViralFlow está online."""
    return {"status": "ok", "timestamp": datetime.utcnow().isoformat()}


@app.post("/jobs", response_model=PipelineJobResponse, status_code=202)
async def submit_job(
    job: PipelineJob,
    background_tasks: BackgroundTasks,
    x_api_key: str = Header(...),
):
    """
    Submete um job de pipeline ao ViralFlow.
    O MacWorp envia os parâmetros; a execução ocorre em background.
    Retorna um job_id para acompanhamento.
    """
    verify_api_key(x_api_key)

    job_id = str(uuid.uuid4())
    job_store.create(job_id, job)

    background_tasks.add_task(run_viralflow_pipeline, job_id, job)

    return PipelineJobResponse(
        job_id=job_id,
        status=JobStatus.QUEUED,
        message="Job submetido com sucesso. Use /jobs/{job_id} para acompanhar.",
        submitted_at=datetime.utcnow(),
    )


@app.get("/jobs/{job_id}", response_model=JobStatusResponse)
def get_job_status(job_id: str, x_api_key: str = Header(...)):
    """
    Retorna o status atual de um job.
    MacWorp faz polling neste endpoint para saber quando o job terminou.
    """
    verify_api_key(x_api_key)

    record = job_store.get(job_id)
    if not record:
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' não encontrado.")

    return record


@app.get("/jobs/{job_id}/results", response_model=JobResult)
def get_job_results(job_id: str, x_api_key: str = Header(...)):
    """
    Retorna os resultados completos de um job concluído.
    Chamado pelo MacWorp após detectar status COMPLETED.
    """
    verify_api_key(x_api_key)

    record = job_store.get(job_id)
    if not record:
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' não encontrado.")

    if record.status not in (JobStatus.COMPLETED, JobStatus.FAILED):
        raise HTTPException(
            status_code=409,
            detail=f"Job ainda em execução (status: {record.status}).",
        )

    return record.result


@app.delete("/jobs/{job_id}", status_code=204)
def cancel_job(job_id: str, x_api_key: str = Header(...)):
    """Cancela um job em execução (envia SIGTERM ao processo Nextflow)."""
    verify_api_key(x_api_key)

    record = job_store.get(job_id)
    if not record:
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' não encontrado.")

    if record.status in (JobStatus.COMPLETED, JobStatus.FAILED):
        raise HTTPException(status_code=409, detail="Job já finalizado.")

    job_store.cancel(job_id)
    return JSONResponse(status_code=204, content={})


@app.get("/jobs")
def list_jobs(x_api_key: str = Header(...)):
    """Lista todos os jobs (útil para debug/monitoramento)."""
    verify_api_key(x_api_key)
    return {"jobs": job_store.list_all()}