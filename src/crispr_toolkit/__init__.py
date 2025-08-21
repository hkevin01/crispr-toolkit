"""
CRISPR Toolkit: AI/ML-powered CRISPR analysis platform for aging research.

This package provides comprehensive tools for CRISPR guide design, analysis,
and optimization with specialized focus on aging and rejuvenation applications.
"""

__version__ = "0.1.0"
__author__ = "CRISPR Toolkit Team"
__email__ = "contact@crisprToolkit.org"

from .analysis.aging.rejuvenation_prediction import predict_rejuvenation
from .analysis.aging.target_prioritization import prioritize_targets
from .analysis.screens import (
    MAGeCKAnalyzer,
    PathwayEnrichment,
    ScreenQualityControl,
    SenescenceScreenAnalyzer,
)
from .cli.aging_cli import main as aging_cli_main

__all__ = [
    "prioritize_targets",
    "predict_rejuvenation",
    "aging_cli_main",
    "MAGeCKAnalyzer",
    "ScreenQualityControl",
    "PathwayEnrichment",
    "SenescenceScreenAnalyzer"
]
