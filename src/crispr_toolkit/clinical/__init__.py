"""
Clinical Intelligence Module for CRISPR Toolkit Phase 3
=======================================================

Advanced clinical trial support and intelligent automation for
aging intervention research with comprehensive EHR integration,
adaptive trial design, and real-time safety monitoring.
"""

from .adaptive_trials import (
    AdaptationDecision,
    AdaptiveTrialEngine,
    TrialOptimizationResult,
)
from .adverse_events import AdverseEventDetector, SafetyAssessment, SafetySignal
from .biomarkers import (
    BiomarkerPanel,
    BiomarkerValidationResult,
    RealTimeBiomarkerEngine,
)
from .ehr_integration import EHRIntegrationEngine, FHIRConnector, PatientRecord
from .stratification import PatientStratificationEngine, StratificationResult

__version__ = "3.0.0"

__all__ = [
    # Patient Stratification
    "PatientStratificationEngine",
    "StratificationResult",

    # Biomarker Discovery
    "RealTimeBiomarkerEngine",
    "BiomarkerPanel",
    "BiomarkerValidationResult",

    # Adaptive Trial Design
    "AdaptiveTrialEngine",
    "TrialOptimizationResult",
    "AdaptationDecision",

    # Safety Monitoring
    "AdverseEventDetector",
    "SafetySignal",
    "SafetyAssessment",

    # EHR Integration
    "EHRIntegrationEngine",
    "FHIRConnector",
    "PatientRecord"
]
    "SafetySignal",
    "SafetyAssessment",

    # EHR Integration
    "EHRIntegrationEngine",
    "FHIRConnector",
    "PatientRecord"
]
