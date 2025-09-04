#!/bin/bash

# CRISPR Aging Toolkit - Main Runner Script
# Quick access to key functionality and testing

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print header
echo -e "${BLUE}🧬 CRISPR Aging Toolkit${NC}"
echo -e "${BLUE}========================${NC}"
echo ""

# Function to display usage
show_usage() {
    echo -e "${YELLOW}Usage: ./run.sh [command]${NC}"
    echo ""
    echo "Available commands:"
    echo "  setup           - Set up the virtual environment and install dependencies"
    echo "  test            - Run comprehensive test suite"
    echo "  test-simple     - Run simple functionality tests"
    echo "  clinical        - Run clinical translator demo"
    echo "  aging-analysis  - Run aging biomarker analysis"
    echo "  mageck          - Run MAGeCK integration demo"
    echo "  docs            - Generate documentation"
    echo "  clean           - Clean up temporary files and caches"
    echo "  help            - Show this help message"
    echo ""
    echo "Examples:"
    echo "  ./run.sh setup     # First-time setup"
    echo "  ./run.sh test      # Run all tests"
    echo "  ./run.sh clinical  # Demo clinical aging assessment"
}

# Activate virtual environment if it exists
activate_venv() {
    if [ -d "venv" ]; then
        echo -e "${YELLOW}🔧 Activating virtual environment...${NC}"
        source venv/bin/activate
    else
        echo -e "${YELLOW}⚠️  Virtual environment not found. Run './run.sh setup' first.${NC}"
        exit 1
    fi
}

# Setup function
setup() {
    echo -e "${YELLOW}🚀 Setting up CRISPR Aging Toolkit...${NC}"

    # Create virtual environment if it doesn't exist
    if [ ! -d "venv" ]; then
        echo -e "${YELLOW}📦 Creating virtual environment...${NC}"
        python3 -m venv venv
    fi

    # Activate and install dependencies
    source venv/bin/activate
    echo -e "${YELLOW}📚 Installing dependencies...${NC}"
    pip install --upgrade pip
    pip install -r requirements.txt

    # Install package in development mode
    echo -e "${YELLOW}🔧 Installing package in development mode...${NC}"
    pip install -e .

    echo -e "${GREEN}✅ Setup complete!${NC}"
}

# Test functions
run_tests() {
    echo -e "${YELLOW}🧪 Running comprehensive test suite...${NC}"
    activate_venv
    python tests/scripts/final_integration_test.py
    echo -e "${GREEN}✅ Tests complete!${NC}"
}

run_simple_tests() {
    echo -e "${YELLOW}🧪 Running simple functionality tests...${NC}"
    activate_venv
    python tests/scripts/quick_test.py
    echo -e "${GREEN}✅ Simple tests complete!${NC}"
}

# Demo functions
run_clinical_demo() {
    echo -e "${YELLOW}🏥 Running clinical translator demo...${NC}"
    activate_venv
    python tests/clinical/test_clinical_translator.py
    echo -e "${GREEN}✅ Clinical demo complete!${NC}"
}

run_aging_analysis() {
    echo -e "${YELLOW}🧬 Running aging biomarker analysis...${NC}"
    activate_venv
    python tests/demos/test_senescence_demo.py
    echo -e "${GREEN}✅ Aging analysis complete!${NC}"
}

run_mageck_demo() {
    echo -e "${YELLOW}📊 Running MAGeCK integration demo...${NC}"
    activate_venv
    python examples/mageck_screen_analysis.py
    echo -e "${GREEN}✅ MAGeCK demo complete!${NC}"
}

# Utility functions
generate_docs() {
    echo -e "${YELLOW}📖 Generating documentation...${NC}"
    activate_venv
    cd docs && make html
    echo -e "${GREEN}✅ Documentation generated in docs/_build/html/index.html${NC}"
}

clean_project() {
    echo -e "${YELLOW}🧹 Cleaning up temporary files...${NC}"
    find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
    find . -name "*.pyc" -delete 2>/dev/null || true
    find . -name "*.pyo" -delete 2>/dev/null || true
    rm -rf .pytest_cache/ .mypy_cache/ *.egg-info/ 2>/dev/null || true
    echo -e "${GREEN}✅ Cleanup complete!${NC}"
}

# Main command handling
case "${1:-help}" in
    "setup")
        setup
        ;;
    "test")
        run_tests
        ;;
    "test-simple")
        run_simple_tests
        ;;
    "clinical")
        run_clinical_demo
        ;;
    "aging-analysis")
        run_aging_analysis
        ;;
    "mageck")
        run_mageck_demo
        ;;
    "docs")
        generate_docs
        ;;
    "clean")
        clean_project
        ;;
    "help"|*)
        show_usage
        ;;
esac
