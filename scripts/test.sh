#!/bin/bash

# CRISPR Toolkit Test Script
# Runs comprehensive tests for the project

set -e

echo "🧪 Running CRISPR Toolkit test suite..."

# Activate virtual environment if it exists
if [ -d "venv" ]; then
    source venv/bin/activate
fi

echo "🔍 Running code quality checks..."

echo "  📏 Checking code formatting with black..."
black --check src/ tests/ || {
    echo "❌ Code formatting issues found. Run 'black src/ tests/' to fix."
    exit 1
}

echo "  📐 Checking import sorting with isort..."
isort --check-only src/ tests/ || {
    echo "❌ Import sorting issues found. Run 'isort src/ tests/' to fix."
    exit 1
}

echo "  🔧 Running flake8 linting..."
flake8 src/ tests/

echo "  🏷️  Running type checking with mypy..."
mypy src/

echo "🔒 Running security checks..."
bandit -r src/

echo "📊 Running unit tests with coverage..."
pytest tests/ --cov=src/ --cov-report=term-missing --cov-report=html

echo "🧬 Running integration tests..."
# Add integration test commands here

echo ""
echo "✅ All tests passed!"
echo "📈 Coverage report generated in htmlcov/"
