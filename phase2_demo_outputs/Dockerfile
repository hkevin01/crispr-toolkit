
FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
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
