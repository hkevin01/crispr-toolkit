"""
Clinical Intelligence Module for CRISPR Toolkit Phase 3
=======================================================

Advanced clinical trial support and patient stratification for
aging intervention research.
"""

from .adaptive_trials import AdaptiveTrialEngine, TrialOptimization
from .biomarkers import BiomarkerDiscoveryEngine, BiomarkerPanel
from .stratification import PatientStratificationEngine, StratificationResult

__version__ = "3.0.0"

__all__ = [
    "PatientStratificationEngine",
    "StratificationResult",
    "BiomarkerDiscoveryEngine",
    "BiomarkerPanel",
    "AdaptiveTrialEngine",
    "TrialOptimization"
]
