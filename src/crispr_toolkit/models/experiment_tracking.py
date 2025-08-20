"""
MLflow experiment tracking for CRISPR aging intervention models.

This module provides experiment tracking, model versioning, and performance
monitoring capabilities using MLflow.
"""

import logging
import mlflow
import mlflow.sklearn
import mlflow.pytorch
import pandas as pd
import numpy as np
from typing import Dict, Any, Optional, List
from pathlib import Path
import json
from datetime import datetime
import hashlib

logger = logging.getLogger(__name__)


class ExperimentTracker:
    """MLflow-based experiment tracking for aging intervention models."""
    
    def __init__(self, experiment_name: str = "crispr_aging_interventions"):
        """Initialize experiment tracker."""
        self.experiment_name = experiment_name
        self.setup_mlflow()
        
    def setup_mlflow(self):
        """Set up MLflow tracking."""
        # Set tracking URI (can be local or remote)
        mlflow_dir = Path("mlruns")
        mlflow_dir.mkdir(exist_ok=True)
        mlflow.set_tracking_uri(f"file://{mlflow_dir.absolute()}")
        
        # Set or create experiment
        try:
            experiment = mlflow.get_experiment_by_name(self.experiment_name)
            if experiment is None:
                experiment_id = mlflow.create_experiment(self.experiment_name)
                logger.info(f"Created new experiment: {self.experiment_name}")
            else:
                experiment_id = experiment.experiment_id
                logger.info(f"Using existing experiment: {self.experiment_name}")
        except Exception as e:
            logger.warning(f"MLflow setup issue: {e}")
            experiment_id = mlflow.create_experiment(self.experiment_name)
            
        mlflow.set_experiment(self.experiment_name)
        
    def start_run(self, run_name: Optional[str] = None, 
                  tags: Optional[Dict[str, str]] = None) -> str:
        """Start a new MLflow run."""
        if run_name is None:
            run_name = f"aging_model_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            
        run = mlflow.start_run(run_name=run_name)
        
        # Add default tags
        default_tags = {
            "model_type": "aging_intervention",
            "framework": "scikit-learn",
            "version": "1.0.0"
        }
        
        if tags:
            default_tags.update(tags)
            
        for key, value in default_tags.items():
            mlflow.set_tag(key, value)
            
        logger.info(f"Started MLflow run: {run.info.run_id}")
        return run.info.run_id
        
    def log_parameters(self, params: Dict[str, Any]):
        """Log model parameters."""
        for key, value in params.items():
            if isinstance(value, (dict, list)):
                # Convert complex types to JSON string
                mlflow.log_param(key, json.dumps(value))
            else:
                mlflow.log_param(key, value)
                
    def log_metrics(self, metrics: Dict[str, float], step: Optional[int] = None):
        """Log model metrics."""
        for key, value in metrics.items():
            if isinstance(value, (int, float)) and not np.isnan(value):
                mlflow.log_metric(key, value, step=step)
            else:
                logger.warning(f"Skipping invalid metric {key}: {value}")
                
    def log_model(self, model, model_name: str = "aging_model", 
                  signature=None, input_example=None):
        """Log trained model."""
        try:
            mlflow.sklearn.log_model(
                sk_model=model,
                artifact_path=model_name,
                signature=signature,
                input_example=input_example
            )
            logger.info(f"Logged model: {model_name}")
        except Exception as e:
            logger.error(f"Failed to log model: {e}")
            
    def log_dataset_info(self, dataset_info: Dict[str, Any]):
        """Log dataset information."""
        # Create dataset hash for versioning
        dataset_str = json.dumps(dataset_info, sort_keys=True)
        dataset_hash = hashlib.md5(dataset_str.encode()).hexdigest()[:8]
        
        mlflow.set_tag("dataset_version", dataset_hash)
        mlflow.log_params({
            "dataset_samples": dataset_info.get("n_samples", 0),
            "dataset_features": dataset_info.get("n_features", 0),
            "dataset_source": dataset_info.get("source", "unknown")
        })
        
    def log_artifacts(self, artifacts: Dict[str, Any]):
        """Log artifacts like plots, data files, etc."""
        temp_dir = Path("temp_artifacts")
        temp_dir.mkdir(exist_ok=True)
        
        try:
            for name, artifact in artifacts.items():
                if isinstance(artifact, pd.DataFrame):
                    # Save DataFrame as CSV
                    file_path = temp_dir / f"{name}.csv"
                    artifact.to_csv(file_path, index=False)
                    mlflow.log_artifact(str(file_path))
                elif isinstance(artifact, dict):
                    # Save dict as JSON
                    file_path = temp_dir / f"{name}.json"
                    with open(file_path, 'w') as f:
                        json.dump(artifact, f, indent=2, default=str)
                    mlflow.log_artifact(str(file_path))
                elif isinstance(artifact, (str, Path)):
                    # Log file directly
                    mlflow.log_artifact(str(artifact))
                    
        except Exception as e:
            logger.error(f"Failed to log artifacts: {e}")
        finally:
            # Cleanup temp files
            if temp_dir.exists():
                import shutil
                shutil.rmtree(temp_dir)
                
    def end_run(self):
        """End current MLflow run."""
        mlflow.end_run()
        logger.info("Ended MLflow run")
        
    def compare_runs(self, run_ids: List[str]) -> pd.DataFrame:
        """Compare multiple runs."""
        runs_data = []
        
        for run_id in run_ids:
            run = mlflow.get_run(run_id)
            run_data = {
                "run_id": run_id,
                "start_time": run.info.start_time,
                "status": run.info.status,
                **run.data.params,
                **run.data.metrics
            }
            runs_data.append(run_data)
            
        return pd.DataFrame(runs_data)
        
    def get_best_run(self, metric_name: str, 
                     ascending: bool = False) -> mlflow.entities.Run:
        """Get best run based on a metric."""
        experiment = mlflow.get_experiment_by_name(self.experiment_name)
        runs = mlflow.search_runs(
            experiment_ids=[experiment.experiment_id],
            order_by=[f"metrics.{metric_name} {'ASC' if ascending else 'DESC'}"]
        )
        
        if len(runs) == 0:
            raise ValueError("No runs found in experiment")
            
        best_run_id = runs.iloc[0]['run_id']
        return mlflow.get_run(best_run_id)
        
    def load_model(self, run_id: str, model_name: str = "aging_model"):
        """Load model from MLflow run."""
        try:
            model_uri = f"runs:/{run_id}/{model_name}"
            model = mlflow.sklearn.load_model(model_uri)
            logger.info(f"Loaded model from run {run_id}")
            return model
        except Exception as e:
            logger.error(f"Failed to load model: {e}")
            return None


class ModelVersioning:
    """Semantic versioning for aging intervention models."""
    
    def __init__(self, model_registry_name: str = "aging_intervention_model"):
        """Initialize model versioning."""
        self.model_registry_name = model_registry_name
        
    def register_model(self, run_id: str, model_name: str = "aging_model",
                      description: str = None) -> str:
        """Register model in MLflow Model Registry."""
        try:
            model_uri = f"runs:/{run_id}/{model_name}"
            
            # Register model
            model_version = mlflow.register_model(
                model_uri=model_uri,
                name=self.model_registry_name,
                description=description
            )
            
            logger.info(f"Registered model version {model_version.version}")
            return model_version.version
            
        except Exception as e:
            logger.error(f"Failed to register model: {e}")
            return None
            
    def transition_model_stage(self, version: str, stage: str):
        """Transition model to different stage (Staging, Production, Archived)."""
        try:
            client = mlflow.tracking.MlflowClient()
            client.transition_model_version_stage(
                name=self.model_registry_name,
                version=version,
                stage=stage
            )
            logger.info(f"Transitioned model v{version} to {stage}")
        except Exception as e:
            logger.error(f"Failed to transition model: {e}")
            
    def get_production_model(self):
        """Get current production model."""
        try:
            client = mlflow.tracking.MlflowClient()
            production_models = client.get_latest_versions(
                self.model_registry_name, 
                stages=["Production"]
            )
            
            if production_models:
                latest_prod = production_models[0]
                model_uri = f"models:/{self.model_registry_name}/{latest_prod.version}"
                return mlflow.sklearn.load_model(model_uri)
            else:
                logger.warning("No production model found")
                return None
                
        except Exception as e:
            logger.error(f"Failed to get production model: {e}")
            return None


class PerformanceDriftDetector:
    """Detect performance drift in aging intervention models."""
    
    def __init__(self, baseline_metrics: Dict[str, float], 
                 drift_threshold: float = 0.05):
        """Initialize drift detector."""
        self.baseline_metrics = baseline_metrics
        self.drift_threshold = drift_threshold
        self.drift_history = []
        
    def detect_drift(self, current_metrics: Dict[str, float]) -> Dict[str, Any]:
        """Detect performance drift."""
        drift_detected = False
        drift_details = {}
        
        for metric_name, baseline_value in self.baseline_metrics.items():
            if metric_name in current_metrics:
                current_value = current_metrics[metric_name]
                drift_magnitude = abs(current_value - baseline_value) / baseline_value
                
                drift_details[metric_name] = {
                    "baseline": baseline_value,
                    "current": current_value,
                    "drift_magnitude": drift_magnitude,
                    "drift_detected": drift_magnitude > self.drift_threshold
                }
                
                if drift_magnitude > self.drift_threshold:
                    drift_detected = True
                    
        # Record drift check
        drift_check = {
            "timestamp": datetime.now().isoformat(),
            "drift_detected": drift_detected,
            "details": drift_details
        }
        
        self.drift_history.append(drift_check)
        
        if drift_detected:
            logger.warning(f"Performance drift detected: {drift_details}")
        else:
            logger.info("No significant performance drift detected")
            
        return drift_check
        
    def get_drift_summary(self) -> Dict[str, Any]:
        """Get summary of drift detection history."""
        if not self.drift_history:
            return {"total_checks": 0, "drift_incidents": 0}
            
        total_checks = len(self.drift_history)
        drift_incidents = sum(1 for check in self.drift_history 
                             if check["drift_detected"])
        
        return {
            "total_checks": total_checks,
            "drift_incidents": drift_incidents,
            "drift_rate": drift_incidents / total_checks,
            "last_check": self.drift_history[-1]["timestamp"],
            "recent_drift": self.drift_history[-1]["drift_detected"]
        }


def setup_experiment_tracking(experiment_name: str = "crispr_aging_interventions") -> ExperimentTracker:
    """Convenience function to set up experiment tracking."""
    return ExperimentTracker(experiment_name)
