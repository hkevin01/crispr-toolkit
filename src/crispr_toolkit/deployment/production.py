"""
Production Deployment Infrastructure for CRISPR Toolkit

This module provides cloud deployment capabilities, REST API endpoints,
and production monitoring for the aging intervention research platform.
"""

import json
import logging
import os
import pickle
from contextlib import asynccontextmanager
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

# Production dependencies
try:
    import uvicorn
    from fastapi import BackgroundTasks, Depends, FastAPI, HTTPException
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
    from pydantic import BaseModel, Field
    FASTAPI_AVAILABLE = True
except ImportError:
    FASTAPI_AVAILABLE = False

import numpy as np

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class DeploymentConfig:
    """Configuration for production deployment"""
    host: str = "0.0.0.0"
    port: int = 8000
    workers: int = 4
    log_level: str = "info"
    reload: bool = False
    max_request_size: int = 50 * 1024 * 1024  # 50MB
    cors_origins: List[str] = None
    auth_enabled: bool = True
    model_cache_size: int = 100
    rate_limit_per_minute: int = 60


class ProductionModelManager:
    """Manages ML models in production environment"""

    def __init__(self, model_dir: str = "models/production"):
        self.model_dir = Path(model_dir)
        self.model_dir.mkdir(parents=True, exist_ok=True)
        self.loaded_models = {}
        self.model_metadata = {}

    def load_model(self, model_name: str, version: str = "latest") -> Any:
        """Load a model for production use"""
        model_key = f"{model_name}_{version}"

        if model_key in self.loaded_models:
            return self.loaded_models[model_key]

        model_path = self.model_dir / model_name / version / "model.pkl"

        if not model_path.exists():
            raise FileNotFoundError(f"Model {model_name} v{version} not found")

        try:
            with open(model_path, 'rb') as f:
                model = pickle.load(f)

            self.loaded_models[model_key] = model

            # Load metadata
            metadata_path = model_path.parent / "metadata.json"
            if metadata_path.exists():
                with open(metadata_path, 'r') as f:
                    self.model_metadata[model_key] = json.load(f)

            logger.info(f"Loaded model {model_name} v{version}")
            return model

        except Exception as e:
            logger.error(f"Failed to load model {model_name}: {e}")
            raise

    def predict(self, model_name: str, data: Dict[str, Any],
                version: str = "latest") -> Dict[str, Any]:
        """Make predictions using loaded model"""
        model = self.load_model(model_name, version)

        try:
            # Convert input data to appropriate format
            if isinstance(data.get('features'), list):
                features = np.array(data['features']).reshape(1, -1)
            elif isinstance(data.get('features'), dict):
                # Convert dict to feature vector based on model metadata
                features = self._dict_to_features(data['features'], model_name, version)
            else:
                raise ValueError("Invalid feature format")

            # Make prediction
            prediction = model.predict(features)
            probability = None

            # Get probability if available
            if hasattr(model, 'predict_proba'):
                probability = model.predict_proba(features)

            return {
                'prediction': prediction.tolist() if hasattr(prediction, 'tolist') else prediction,
                'probability': probability.tolist() if probability is not None else None,
                'model_name': model_name,
                'model_version': version,
                'timestamp': datetime.now().isoformat()
            }

        except Exception as e:
            logger.error(f"Prediction failed for {model_name}: {e}")
            raise

    def _dict_to_features(self, feature_dict: Dict[str, float],
                         model_name: str, version: str) -> np.ndarray:
        """Convert feature dictionary to numpy array"""
        model_key = f"{model_name}_{version}"
        metadata = self.model_metadata.get(model_key, {})
        feature_names = metadata.get('feature_names', [])

        if not feature_names:
            raise ValueError(f"No feature names found for {model_name}")

        features = np.array([feature_dict.get(name, 0.0) for name in feature_names])
        return features.reshape(1, -1)


class ProductionMonitor:
    """Monitors production system health and performance"""

    def __init__(self, log_dir: str = "logs/production"):
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.metrics = {
            'requests_count': 0,
            'errors_count': 0,
            'average_response_time': 0.0,
            'models_loaded': 0,
            'cache_hit_rate': 0.0
        }

    def log_request(self, endpoint: str, response_time: float,
                   status_code: int, user_id: Optional[str] = None):
        """Log API request metrics"""
        self.metrics['requests_count'] += 1

        if status_code >= 400:
            self.metrics['errors_count'] += 1

        # Update average response time
        current_avg = self.metrics['average_response_time']
        count = self.metrics['requests_count']
        self.metrics['average_response_time'] = (
            (current_avg * (count - 1) + response_time) / count
        )

        # Log to file
        log_entry = {
            'timestamp': datetime.now().isoformat(),
            'endpoint': endpoint,
            'response_time': response_time,
            'status_code': status_code,
            'user_id': user_id
        }

        log_file = self.log_dir / f"requests_{datetime.now().strftime('%Y%m%d')}.json"
        with open(log_file, 'a') as f:
            f.write(json.dumps(log_entry, default=str) + '\n')

    def get_health_status(self) -> Dict[str, Any]:
        """Get current system health status"""
        error_rate = (self.metrics['errors_count'] / self.metrics['requests_count']
                     if self.metrics['requests_count'] > 0 else 0)

        health_status = "healthy"
        if error_rate > 0.1:  # >10% error rate
            health_status = "unhealthy"
        elif error_rate > 0.05:  # >5% error rate
            health_status = "degraded"

        return {
            'status': health_status,
            'metrics': self.metrics,
            'error_rate': error_rate,
            'timestamp': datetime.now().isoformat()
        }


if FASTAPI_AVAILABLE:
    # Pydantic models for API
    class PredictionRequest(BaseModel):
        model_name: str = Field(..., description="Name of the model to use")
        features: Dict[str, float] = Field(..., description="Feature values")
        version: str = Field("latest", description="Model version")

    class ValidationRequest(BaseModel):
        dataset_accession: str = Field(..., description="GEO dataset accession")
        predictions: Dict[str, float] = Field(..., description="Gene predictions")

    class HealthResponse(BaseModel):
        status: str
        metrics: Dict[str, Any]
        timestamp: str

    # Global instances
    model_manager = ProductionModelManager()
    monitor = ProductionMonitor()
    security = HTTPBearer()

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        # Startup
        logger.info("Starting CRISPR Toolkit Production API")
        yield
        # Shutdown
        logger.info("Shutting down CRISPR Toolkit Production API")

    # Create FastAPI app
    app = FastAPI(
        title="CRISPR Toolkit for Aging Research",
        description="Production API for aging intervention predictions",
        version="2.0.0",
        lifespan=lifespan
    )

    # CORS middleware
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # Configure properly for production
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    async def verify_token(credentials: HTTPAuthorizationCredentials = Depends(security)):
        """Verify authentication token"""
        # In production, implement proper JWT verification
        token = credentials.credentials
        if token != os.getenv("API_TOKEN", "demo_token_12345"):
            raise HTTPException(status_code=401, detail="Invalid authentication token")
        return token

    @app.get("/health", response_model=HealthResponse)
    async def health_check():
        """Health check endpoint"""
        start_time = datetime.now()
        health_status = monitor.get_health_status()
        response_time = (datetime.now() - start_time).total_seconds()

        monitor.log_request("/health", response_time, 200)
        return health_status

    @app.post("/predict")
    async def predict_intervention(
        request: PredictionRequest,
        background_tasks: BackgroundTasks,
        token: str = Depends(verify_token)
    ):
        """Predict aging intervention effects"""
        start_time = datetime.now()

        try:
            result = model_manager.predict(
                model_name=request.model_name,
                data={'features': request.features},
                version=request.version
            )

            response_time = (datetime.now() - start_time).total_seconds()
            background_tasks.add_task(
                monitor.log_request, "/predict", response_time, 200
            )

            return result

        except Exception as e:
            response_time = (datetime.now() - start_time).total_seconds()
            background_tasks.add_task(
                monitor.log_request, "/predict", response_time, 500
            )
            raise HTTPException(status_code=500, detail=str(e))

    @app.post("/validate")
    async def validate_predictions(
        request: ValidationRequest,
        background_tasks: BackgroundTasks,
        token: str = Depends(verify_token)
    ):
        """Validate predictions against real datasets"""
        start_time = datetime.now()

        try:
            # Import validation framework
            from ..validation.geo_datasets import AgingValidationFramework

            validator = AgingValidationFramework()
            result = validator.validate_predictions(
                request.dataset_accession,
                request.predictions
            )

            response_time = (datetime.now() - start_time).total_seconds()
            background_tasks.add_task(
                monitor.log_request, "/validate", response_time, 200
            )

            return result

        except Exception as e:
            response_time = (datetime.now() - start_time).total_seconds()
            background_tasks.add_task(
                monitor.log_request, "/validate", response_time, 500
            )
            raise HTTPException(status_code=500, detail=str(e))

    @app.get("/models")
    async def list_models(token: str = Depends(verify_token)):
        """List available models"""
        try:
            models = []
            for model_dir in model_manager.model_dir.iterdir():
                if model_dir.is_dir():
                    versions = [v.name for v in model_dir.iterdir() if v.is_dir()]
                    models.append({
                        'name': model_dir.name,
                        'versions': versions
                    })

            return {'models': models}

        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))

    @app.get("/datasets")
    async def list_datasets(token: str = Depends(verify_token)):
        """List available validation datasets"""
        try:
            from ..validation.geo_datasets import GEODatasetLoader

            loader = GEODatasetLoader()
            datasets = []

            for accession in loader.list_available_datasets():
                dataset_info = loader.get_dataset_info(accession)
                datasets.append({
                    'accession': accession,
                    'title': dataset_info.title,
                    'intervention_type': dataset_info.intervention_type,
                    'species': dataset_info.species,
                    'tissue': dataset_info.tissue
                })

            return {'datasets': datasets}

        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))


class CloudDeployment:
    """Handles cloud deployment configuration and management"""

    def __init__(self, provider: str = "aws"):
        self.provider = provider
        self.config = DeploymentConfig()

    def generate_docker_config(self) -> str:
        """Generate Dockerfile for containerized deployment"""
        dockerfile = """
FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \\
    build-essential \\
    && rm -rf /var/lib/apt/lists/*

# Copy requirements
COPY requirements_production.txt .
RUN pip install --no-cache-dir -r requirements_production.txt

# Copy application
COPY . .

# Create necessary directories
RUN mkdir -p models/production logs/production data/geo_cache

# Expose port
EXPOSE 8000

# Run application
CMD ["python", "-m", "uvicorn", "src.crispr_toolkit.deployment.production:app", "--host", "0.0.0.0", "--port", "8000"]
"""
        return dockerfile

    def generate_docker_compose(self) -> str:
        """Generate docker-compose.yml for multi-service deployment"""
        compose = """
version: '3.8'

services:
  crispr-api:
    build: .
    ports:
      - "8000:8000"
    environment:
      - API_TOKEN=${API_TOKEN}
      - LOG_LEVEL=info
    volumes:
      - ./models:/app/models
      - ./logs:/app/logs
      - ./data:/app/data
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3

  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf
      - ./ssl:/etc/nginx/ssl
    depends_on:
      - crispr-api
    restart: unless-stopped

  monitoring:
    image: grafana/grafana:latest
    ports:
      - "3000:3000"
    environment:
      - GF_SECURITY_ADMIN_PASSWORD=${GRAFANA_PASSWORD}
    volumes:
      - grafana-data:/var/lib/grafana
    restart: unless-stopped

volumes:
  grafana-data:
"""
        return compose

    def generate_kubernetes_config(self) -> Dict[str, str]:
        """Generate Kubernetes deployment configuration"""

        deployment = """
apiVersion: apps/v1
kind: Deployment
metadata:
  name: crispr-toolkit-api
  labels:
    app: crispr-toolkit
spec:
  replicas: 3
  selector:
    matchLabels:
      app: crispr-toolkit
  template:
    metadata:
      labels:
        app: crispr-toolkit
    spec:
      containers:
      - name: api
        image: crispr-toolkit:latest
        ports:
        - containerPort: 8000
        env:
        - name: API_TOKEN
          valueFrom:
            secretKeyRef:
              name: crispr-secrets
              key: api-token
        resources:
          requests:
            memory: "1Gi"
            cpu: "500m"
          limits:
            memory: "2Gi"
            cpu: "1000m"
        livenessProbe:
          httpGet:
            path: /health
            port: 8000
          initialDelaySeconds: 30
          periodSeconds: 10
        readinessProbe:
          httpGet:
            path: /health
            port: 8000
          initialDelaySeconds: 5
          periodSeconds: 5
"""

        service = """
apiVersion: v1
kind: Service
metadata:
  name: crispr-toolkit-service
spec:
  selector:
    app: crispr-toolkit
  ports:
  - protocol: TCP
    port: 80
    targetPort: 8000
  type: LoadBalancer
"""

        return {
            'deployment.yaml': deployment,
            'service.yaml': service
        }


def create_production_requirements():
    """Create production requirements file"""
    requirements = """
# Production dependencies for CRISPR Toolkit
fastapi==0.104.1
uvicorn[standard]==0.24.0
pydantic==2.5.0
pandas==2.1.3
numpy==1.24.3
scikit-learn==1.3.2
mlflow==2.8.1
optuna==3.4.0
requests==2.31.0
python-multipart==0.0.6
python-jose[cryptography]==3.3.0
passlib[bcrypt]==1.7.4

# Monitoring and logging
prometheus-client==0.19.0
structlog==23.2.0

# Database (if needed)
sqlalchemy==2.0.23
psycopg2-binary==2.9.9

# Caching
redis==5.0.1

# Development and testing
pytest==7.4.3
pytest-asyncio==0.21.1
httpx==0.25.2
"""
    return requirements


if __name__ == "__main__":
    if not FASTAPI_AVAILABLE:
        print("FastAPI not available. Install with: pip install fastapi uvicorn")
        exit(1)

    # Run the production server
    config = DeploymentConfig()
    uvicorn.run(
        "src.crispr_toolkit.deployment.production:app",
        host=config.host,
        port=config.port,
        workers=config.workers,
        log_level=config.log_level,
        reload=config.reload
    )
