"""
Advanced Aging Biomarker Analysis - Phase 2B

This module provides cutting-edge single-cell multi-omics and AI-driven strategies
for precision aging biology, implementing the latest 2025 research developments.

Features:
- Single-cell aging clocks (scAge)
- Spatial aging mapping
- AI-powered biomarker discovery
- Non-linear aging models
- Senescence heterogeneity analysis
- Multi-omics integration
- Clinical translation tools

Components:
- SingleCellAgingAnalyzer: Cell-type specific aging analysis
- SpatialAgingMapper: Tissue architecture aging assessment
- AIBiomarkerDiscovery: ML-driven biomarker identification
- SenescenceClassifier: Advanced senescence detection
- MultiOmicsIntegrator: Cross-platform data fusion
- ClinicalTranslator: Research to clinic translation

Author: CRISPR Toolkit Team
Version: 2B.0.0
Date: August 21, 2025
"""

from .single_cell_aging import SingleCellAgingAnalyzer
from .spatial_aging import SpatialAgingMapper
from .ai_biomarker_discovery import AIBiomarkerDiscovery
from .senescence_classifier import SenescenceClassifier
from .multi_omics_integrator import MultiOmicsIntegrator
from .clinical_translator import ClinicalTranslator

__all__ = [
    'SingleCellAgingAnalyzer',
    'SpatialAgingMapper', 
    'AIBiomarkerDiscovery',
    'SenescenceClassifier',
    'MultiOmicsIntegrator',
    'ClinicalTranslator'
]

__version__ = "2B.0.0"
