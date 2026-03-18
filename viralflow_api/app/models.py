"""Modelos Pydantic para a ViralFlow API."""

from enum import Enum
from typing import Optional, Dict, Any, List
from datetime import datetime
from pydantic import BaseModel, Field


class JobStatus(str, Enum):
    QUEUED    = "queued"
    RUNNING   = "running"
    COMPLETED = "completed"
    FAILED    = "failed"
    CANCELLED = "cancelled"


# ── Request: MacWorp → ViralFlow ──────────────────────────────────────

class PipelineJob(BaseModel):
    """Payload enviado pelo MacWorp proxy (main.nf) para a API ViralFlow."""

    macworp_project_id: str
    macworp_run_id:     str

    # Caminhos NFS
    in_dir:  str = Field(..., description="Pasta com FASTQs (caminho NFS)")
    out_dir: str = Field(..., description="Pasta de saída no NFS (projeto MacWorp)")

    # Parâmetros principais
    virus:              str  = "sars-cov2"
    depth:              int  = 25
    min_len:            int  = 75
    run_snpeff:         bool = True
    write_mapped_reads: bool = True
    min_dp_intrahost:   int  = 100

    # Parâmetros extras — todos os campos adicionais do workflow definition
    extra_params: Optional[Dict[str, Any]] = Field(
        default=None,
        description=(
            "Parâmetros adicionais: dedup, trimLen, ndedup, nextflowSimCalls, "
            "fastp_threads, bwa_threads, mafft_threads, "
            "refGenomeCode, referenceGFF, referenceGenome"
        )
    )

    # Callback opcional
    callback_url:     Optional[str] = None
    callback_api_key: Optional[str] = None


# ── Responses ─────────────────────────────────────────────────────────

class PipelineJobResponse(BaseModel):
    job_id:       str
    status:       JobStatus
    message:      str
    submitted_at: datetime


class OutputFile(BaseModel):
    name:        str
    path:        str
    size_bytes:  Optional[int] = None
    description: Optional[str] = None


class JobResult(BaseModel):
    job_id:             str
    macworp_project_id: str
    macworp_run_id:     str
    status:             JobStatus
    started_at:         Optional[datetime] = None
    finished_at:        Optional[datetime] = None
    duration_seconds:   Optional[float]    = None
    exit_code:          Optional[int]      = None
    output_dir:         Optional[str]      = None
    output_files:       List[OutputFile]   = []
    nextflow_log:       Optional[str]      = None
    error_message:      Optional[str]      = None
    metrics:            Optional[Dict[str, Any]] = None


class JobStatusResponse(BaseModel):
    job_id:             str
    macworp_project_id: str
    macworp_run_id:     str
    status:             JobStatus
    submitted_at:       datetime
    started_at:         Optional[datetime] = None
    finished_at:        Optional[datetime] = None
    progress_message:   Optional[str]      = None
    result:             Optional[JobResult] = None
