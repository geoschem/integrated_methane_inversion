#!/bin/bash
# Helper script to run tests
# Sets up conda environment and PYTHONPATH automatically
# NOTE: I (Vedant) can't actually get this script to run because of some `conda init` error

set -e

# Activate imi_env
conda activate imi_env

# Set up Python path
export PYTHONPATH="$(pwd):$PYTHONPATH"

# Parse arguments
if [ "$1" = "coverage" ]; then
    # Run with coverage report
    pytest src --cov=src --cov-report=html --cov-report=term -v
    echo "✓ Coverage report generated in htmlcov/index.html"
elif [ "$1" = "quick" ]; then
    # Quick run without verbose output
    pytest src -q
elif [ "$1" = "watch" ]; then
    # Watch mode (requires pytest-watch: pip install pytest-watch)
    ptw src
elif [ "$1" = "help" ] || [ "$1" = "-h" ]; then
    cat << EOF
Usage: ./run_tests.sh [option]

Options:
    (no option)  Run all tests with verbose output
    coverage     Run tests and generate coverage report
    quick        Run tests with minimal output
    watch        Watch mode (re-run on file changes)
    help, -h     Show this help message

Examples:
    ./run_tests.sh              # Run all tests
    ./run_tests.sh coverage     # Run with coverage
    ./run_tests.sh quick        # Quick test run

Environment:
    Automatically activates 'imi_env' conda environment
    Sets PYTHONPATH to workspace root

EOF
else
    # Default: run all tests with verbose output
    pytest src -v "$@"
fi
