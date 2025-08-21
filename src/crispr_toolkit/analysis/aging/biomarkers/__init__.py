"""
Aging Biomarker Integration Module

This module integrates pyaging and biolearn libraries to provide comprehensive
aging biomarker analysis capabilities, with specific focus on ReHMGB1/RAGE
signaling pathways and senescence research.

Key Features:
- 100+ aging clocks via pyaging integration
- Standardized aging biomarkers via biolearn
- ReHMGB1-specific aging biomarker analysis
- Clinical significance scoring
- Multi-modal aging biomarker analysis
- Senescence pathway-specific biomarkers

Author: CRISPR Toolkit Development Team
"""

from .biolearn_integration import BiolearnAnalyzer
from .biomarker_analysis import AgingBiomarkerAnalyzer
from .clinical_scoring import ClinicalAgingScorer
from .pyaging_integration import PyAgingAnalyzer
from .rehmgb1_biomarkers import ReHMGB1BiomarkerAnalyzer
from .visualization import AgingBiomarkerVisualizer

__all__ = [
    'PyAgingAnalyzer',
    'BiolearnAnalyzer',
    'AgingBiomarkerAnalyzer',
    'ReHMGB1BiomarkerAnalyzer',
    'ClinicalAgingScorer',
    'AgingBiomarkerVisualizer'
]

__version__ = '1.0.0'
