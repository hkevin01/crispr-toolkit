#!/bin/bash

# CRISPR Toolkit Development Setup Script
# This script sets up the development environment

set -e

echo "🧬 Setting up CRISPR Toolkit development environment..."

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is required but not installed."
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
REQUIRED_VERSION="3.8"

if [ "$(printf '%s\n' "$REQUIRED_VERSION" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED_VERSION" ]; then
    echo "❌ Python $REQUIRED_VERSION or higher is required. Found: $PYTHON_VERSION"
    exit 1
fi

echo "✅ Python $PYTHON_VERSION found"

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "📦 Creating virtual environment..."
    python3 -m venv venv
fi

echo "🔧 Activating virtual environment..."
source venv/bin/activate

echo "📥 Installing dependencies..."
pip install --upgrade pip
pip install -r requirements.txt
pip install -r requirements-dev.txt

echo "🔨 Installing package in development mode..."
pip install -e .

echo "🧪 Setting up pre-commit hooks..."
pre-commit install

echo "🐳 Setting up Docker environment..."
if command -v docker &> /dev/null; then
    if [ -f "docker-compose.yml" ]; then
        docker-compose build
        echo "✅ Docker environment ready"
    else
        echo "⚠️  docker-compose.yml not found, skipping Docker setup"
    fi
else
    echo "⚠️  Docker not found, skipping Docker setup"
fi

echo "🎯 Running initial tests..."
pytest tests/ -v || echo "⚠️  Some tests failed - this is expected for a new setup"

echo ""
echo "🎉 Setup complete! Next steps:"
echo "  1. Activate the environment: source venv/bin/activate"
echo "  2. Run tests: pytest"
echo "  3. Start development: code ."
echo "  4. Check the CLI: crispr-aging --help"
echo ""
echo "📚 Documentation: docs/"
echo "🐛 Issues: https://github.com/username/crispr-toolkit/issues"
