"""
CRISPR Screen Analysis Module

This module provides comprehensive tools for analyzing CRISPR screens with
specialized focus on aging and senescence pathways. Integrates MAGeCK and
other industry-standard tools for robust screen analysis.

Key Features:
- MAGeCK integration for screen analysis
- Quality control metrics
- Pathway enrichment analysis
- Senescence marker analysis
- ReHMGB1 pathway investigation
"""

from .mageck_wrapper import MAGeCKAnalyzer
from .pathway_analysis import PathwayEnrichment
from .screen_qc import ScreenQualityControl
from .senescence_analysis import SenescenceScreenAnalyzer

__all__ = [
    "MAGeCKAnalyzer",
    "ScreenQualityControl",
    "PathwayEnrichment",
    "SenescenceScreenAnalyzer"
]
