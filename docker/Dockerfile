# Use Python 3.10 slim image as base
FROM python:3.10-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create and set working directory
WORKDIR /app

# Create non-root user
RUN groupadd -r crispr && useradd -r -g crispr crispr

# Copy requirements first for better caching
COPY requirements.txt requirements-dev.txt ./

# Install Python dependencies
RUN pip install --upgrade pip && \
    pip install -r requirements.txt && \
    pip install -r requirements-dev.txt

# Copy project files
COPY . .

# Install the package
RUN pip install -e .

# Create virtual environment for isolated execution
RUN python -m venv /app/venv

# Change ownership to non-root user
RUN chown -R crispr:crispr /app

# Switch to non-root user
USER crispr

# Expose port
EXPOSE 8000

# Set default command
CMD ["uvicorn", "crispr_toolkit.api.main:app", "--host", "0.0.0.0", "--port", "8000"]
