#!/bin/bash

# CRISPR Toolkit Aging Research Setup Script
# This script sets up the complete development environment

set -e

echo "🧬 Setting up CRISPR Toolkit for aging research..."

# Check Python version
python_version=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
required_version="3.8"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" != "$required_version" ]; then
    echo "❌ Python $required_version or higher is required. Found $python_version"
    exit 1
fi

echo "✅ Python version: $python_version"

# Create virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "📦 Creating virtual environment..."
    python3 -m venv venv
fi

# Activate virtual environment
echo "🔄 Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "⬆️ Upgrading pip..."
pip install --upgrade pip

# Install dependencies
echo "📚 Installing dependencies..."
pip install -r requirements.txt

# Install package in development mode
echo "🔧 Installing CRISPR Toolkit in development mode..."
pip install -e .

# Create necessary directories
echo "📁 Creating output directories..."
mkdir -p outputs/prioritization
mkdir -p outputs/predictions
mkdir -p data/expression
mkdir -p data/literature
mkdir -p models/aging
mkdir -p logs

# Create .env file with examples
if [ ! -f ".env" ]; then
    echo "⚙️ Creating .env configuration file..."
    cat > .env << 'EOF'
# CRISPR Toolkit Configuration

# API Keys (replace with your actual keys)
PUBMED_API_KEY=your_pubmed_api_key_here
STRING_API_KEY=your_string_api_key_here

# Database Configuration
DATABASE_URL=postgresql://user:password@localhost/crispr_toolkit
REDIS_URL=redis://localhost:6379

# Model Paths
AGING_MODEL_PATH=models/aging/rejuvenation_predictor.joblib
PRIORITIZATION_MODEL_PATH=models/aging/target_prioritizer.joblib

# Logging
LOG_LEVEL=INFO
LOG_FILE=logs/crispr_toolkit.log

# Development Settings
DEBUG=false
ENVIRONMENT=development
EOF
    echo "📝 Created .env file with example configuration"
    echo "⚠️  Please update the API keys and database URLs in .env"
fi

# Test basic imports
echo "🔬 Running functionality tests..."

# Test target prioritization
echo "  Testing target prioritization..."
python -c "
from crispr_toolkit.analysis.aging.prioritization import AgingTargetPrioritizer
prioritizer = AgingTargetPrioritizer()
genes = ['TP53', 'CDKN2A', 'SIRT1']
results = prioritizer.rank_targets(genes, tissue='liver', n_targets=3)
print(f'✅ Prioritized {len(results)} targets')
" 2>/dev/null || echo "⚠️  Target prioritization test failed (dependencies may be missing)"

# Test rejuvenation prediction
echo "  Testing rejuvenation prediction..."
python -c "
from crispr_toolkit.analysis.aging.rejuvenation import predict_rejuvenation
config = {'intervention': 'OSK', 'targets': ['TP53'], 'delivery': 'AAV'}
results = predict_rejuvenation(config, tissue='liver')
print(f'✅ Predicted epigenetic clock delta: {results[\"predictions\"][\"epigenetic_clock_delta\"]:.2f}')
" 2>/dev/null || echo "⚠️  Rejuvenation prediction test failed (dependencies may be missing)"

echo ""
echo "🎉 CRISPR Toolkit setup completed!"
echo ""
echo "📋 Next steps:"
echo "1. 🔑 Update API keys in .env file"
echo "2. 🐍 Activate environment: source venv/bin/activate"
echo "3. 🚀 Test CLI: python -m crispr_toolkit.cli.main --help"
echo ""
echo "🧪 Quick examples:"
echo "  Prioritize targets:  python -m crispr_toolkit.cli.main aging prioritize --tissue liver --output outputs/test.json"
echo "  Predict outcomes:    python -m crispr_toolkit.cli.main aging predict --intervention OSK --tissue liver --output outputs/prediction.json"
echo ""
echo "📖 See README.md for detailed documentation"
