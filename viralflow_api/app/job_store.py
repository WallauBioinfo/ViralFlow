"""
Store de jobs em memória (thread-safe).
Para produção, substitua por Redis ou banco de dados.
"""

import threading
from typing import Dict, Optional, List
from datetime import datetime

from app.models import JobStatus, JobStatusResponse, PipelineJob, JobResult


class JobStore:
    def __init__(self):
        self._store: Dict[str, JobStatusResponse] = {}
        self._processes: Dict[str, object] = {}  # job_id → subprocess.Popen
        self._lock = threading.Lock()

    def create(self, job_id: str, job: PipelineJob) -> JobStatusResponse:
        record = JobStatusResponse(
            job_id=job_id,
            macworp_project_id=job.macworp_project_id,
            macworp_run_id=job.macworp_run_id,
            status=JobStatus.QUEUED,
            submitted_at=datetime.utcnow(),
        )
        with self._lock:
            self._store[job_id] = record
        return record

    def get(self, job_id: str) -> Optional[JobStatusResponse]:
        with self._lock:
            return self._store.get(job_id)

    def update_status(
        self,
        job_id: str,
        status: JobStatus,
        progress_message: Optional[str] = None,
        started_at: Optional[datetime] = None,
        finished_at: Optional[datetime] = None,
        result: Optional[JobResult] = None,
    ):
        with self._lock:
            record = self._store.get(job_id)
            if record:
                record.status = status
                if progress_message:
                    record.progress_message = progress_message
                if started_at:
                    record.started_at = started_at
                if finished_at:
                    record.finished_at = finished_at
                if result:
                    record.result = result

    def register_process(self, job_id: str, process):
        with self._lock:
            self._processes[job_id] = process

    def cancel(self, job_id: str):
        with self._lock:
            proc = self._processes.get(job_id)
            if proc:
                try:
                    proc.terminate()
                except Exception:
                    pass
            record = self._store.get(job_id)
            if record:
                record.status = JobStatus.CANCELLED
                record.finished_at = datetime.utcnow()

    def list_all(self) -> List[dict]:
        with self._lock:
            return [
                {
                    "job_id": r.job_id,
                    "macworp_project_id": r.macworp_project_id,
                    "status": r.status,
                    "submitted_at": r.submitted_at.isoformat(),
                }
                for r in self._store.values()
            ]


# Singleton compartilhado pela aplicação
job_store = JobStore()
