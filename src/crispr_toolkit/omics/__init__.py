"""
Multi-Omics Integration Framework for CRISPR Toolkit Phase 3
============================================================

Advanced multi-omics analysis capabilities for aging intervention research
including proteomics, metabolomics, epigenomics, multi-omics fusion,
pathway enrichment analysis, temporal modeling, cross-omics correlation,
and comprehensive quality control.
"""

from .cross_omics_correlation import (
    CrossOmicsCorrelationAnalyzer,
    IntegrationQualityMetrics,
    MultiOmicsNetworkResult,
)
from .epigenomics import EpigeneticAgeResult, EpigenomicsAnalyzer
from .fusion import MultiOmicsFusion, OmicsDataLayer, TemporalOmicsModel
from .metabolomics import MetabolicPathwayModel, MetabolomicsAnalyzer
from .pathway_enrichment import AgingPathwayScore, PathwayEnrichmentAnalyzer
from .proteomics import ProteinExpressionModel, ProteomicsAnalyzer
from .quality_control import (
    BatchEffectResult,
    MultiOmicsQualityControl,
    NormalizationResult,
    OutlierDetectionResult,
    QualityMetrics,
)
from .temporal_modeling import (
    InterventionTrajectory,
    TemporalModelResult,
    TemporalMultiOmicsModel,
    TemporalOmicsData,
)

__version__ = "3.0.0"

__all__ = [
    # Proteomics
    "ProteomicsAnalyzer",
    "ProteinExpressionModel",
    # Metabolomics
    "MetabolomicsAnalyzer",
    "MetabolicPathwayModel",
    # Epigenomics
    "EpigenomicsAnalyzer",
    "EpigeneticAgeResult",
    # Multi-omics Fusion
    "MultiOmicsFusion",
    "OmicsDataLayer",
    "TemporalOmicsModel",
    # Pathway Analysis
    "PathwayEnrichmentAnalyzer",
    "AgingPathwayScore",
    # Temporal Modeling
    "TemporalMultiOmicsModel",
    "InterventionTrajectory",
    "TemporalModelResult",
    "TemporalOmicsData",
    # Cross-Omics Correlation
    "CrossOmicsCorrelationAnalyzer",
    "MultiOmicsNetworkResult",
    "IntegrationQualityMetrics",
    # Quality Control
    "MultiOmicsQualityControl",
    "QualityMetrics",
    "BatchEffectResult",
    "OutlierDetectionResult",
    "NormalizationResult"
]


__version__ = "3.0.0"

__all__ = [
    # Proteomics
    "ProteomicsAnalyzer",
    "ProteinExpressionModel",
    # Metabolomics
    "MetabolomicsAnalyzer",
    "MetabolicPathwayModel",
    # Epigenomics
    "EpigenomicsAnalyzer",
    "EpigeneticAgeResult",
    # Multi-omics Fusion
    "MultiOmicsFusion",
    "OmicsDataLayer",
    "TemporalOmicsModel",
    # Pathway Analysis
    "PathwayEnrichmentAnalyzer",
    "AgingPathwayScore"
]


__version__ = "3.0.0"

__all__ = [
    # Proteomics
    "ProteomicsAnalyzer",
    "ProteinExpressionModel",
    # Metabolomics
    "MetabolomicsAnalyzer",
    "MetabolicPathwayModel",
    # Epigenomics
    "EpigenomicsAnalyzer",
    "EpigeneticAgeResult",
    # Multi-omics Fusion
    "MultiOmicsFusion",
    "OmicsDataLayer",
    "TemporalOmicsModel",
    # Pathway Analysis
    "PathwayEnrichmentAnalyzer",
    "AgingPathwayScore"
]
