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
- SenescenceClassifier: Advanced senescence detection (planned)
- MultiOmicsIntegrator: Cross-platform data fusion (planned)
- ClinicalTranslator: Research to clinic translation (planned)

Author: CRISPR Toolkit Team
Version: 2B.0.1
Date: August 21, 2025
"""

# Import available components with graceful error handling
try:
    from .single_cell_aging import SingleCellAgingAnalyzer
    SINGLE_CELL_AVAILABLE = True
except ImportError:
    SingleCellAgingAnalyzer = None
    SINGLE_CELL_AVAILABLE = False

try:
    from .spatial_aging import SpatialAgingMapper
    SPATIAL_AVAILABLE = True
except ImportError:
    SpatialAgingMapper = None
    SPATIAL_AVAILABLE = False

try:
    from .ai_biomarker_discovery import AIBiomarkerDiscovery
    AI_DISCOVERY_AVAILABLE = True
except ImportError:
    AIBiomarkerDiscovery = None
    AI_DISCOVERY_AVAILABLE = False

try:
    from .senescence_classifier import SenescenceClassifier
    SENESCENCE_AVAILABLE = True
except ImportError:
    SenescenceClassifier = None
    SENESCENCE_AVAILABLE = False

# Planned components (not yet implemented)
MultiOmicsIntegrator = None
ClinicalTranslator = None

# Export available components
__all__ = []

if SINGLE_CELL_AVAILABLE:
    __all__.append('SingleCellAgingAnalyzer')

if SPATIAL_AVAILABLE:
    __all__.append('SpatialAgingMapper')

if AI_DISCOVERY_AVAILABLE:
    __all__.append('AIBiomarkerDiscovery')

if SENESCENCE_AVAILABLE:
    __all__.append('SenescenceClassifier')

__all__.extend([
    'get_available_components', 'get_component_info', '__version__'
])

__version__ = "2B.0.1"


def get_available_components():
    """Get list of currently available advanced aging analysis components."""
    available = []

    if SINGLE_CELL_AVAILABLE:
        available.append('SingleCellAgingAnalyzer')

    if SPATIAL_AVAILABLE:
        available.append('SpatialAgingMapper')

    if AI_DISCOVERY_AVAILABLE:
        available.append('AIBiomarkerDiscovery')

    if SENESCENCE_AVAILABLE:
        available.append('SenescenceClassifier')

    return available


def get_component_info():
    """Get detailed information about Phase 2B components."""

    info = {
        'version': __version__,
        'phase': '2B - Advanced Aging Analysis',
        'implementation_status': {
            'implemented': len(get_available_components()),
            'planned': 3,
            'total': 6
        },
        'components': {
            'SingleCellAgingAnalyzer': {
                'status': 'implemented' if SINGLE_CELL_AVAILABLE else 'import_error',
                'description': 'scRNA-seq aging clocks and cellular heterogeneity',
                'features': [
                    'Cell-type specific biological age estimation',
                    'Aging trajectory analysis',
                    'Single-cell senescence scoring',
                    'Cellular aging heterogeneity quantification'
                ]
            },
            'SpatialAgingMapper': {
                'status': 'implemented' if SPATIAL_AVAILABLE else 'import_error',
                'description': 'Spatial transcriptomics aging pattern analysis',
                'features': [
                    'Spatial aging hotspot detection',
                    'Senescence-associated spatial patterns',
                    'Tissue architecture aging analysis',
                    'Spatial aging gradient mapping'
                ]
            },
            'AIBiomarkerDiscovery': {
                'status': 'implemented' if AI_DISCOVERY_AVAILABLE else 'import_error',
                'description': 'AI-driven biomarker discovery with explainable AI',
                'features': [
                    'Non-linear aging models',
                    'SHAP-based explanations',
                    'Sudden aging detection',
                    'Multi-modal biomarker integration'
                ]
            },
            'SenescenceClassifier': {
                'status': 'planned',
                'description': 'Advanced senescence detection and classification',
                'features': [
                    'Multi-modal senescence scoring',
                    'Senolytic target identification',
                    'SASP analysis',
                    'Senescent cell type classification'
                ]
            },
            'MultiOmicsIntegrator': {
                'status': 'planned',
                'description': 'Multi-omics aging data integration',
                'features': [
                    'Cross-omics aging signatures',
                    'Pathway-level aging analysis',
                    'Network-based aging modules',
                    'Regulatory aging circuits'
                ]
            },
            'ClinicalTranslator': {
                'status': 'planned',
                'description': 'Clinical aging assessment and translation',
                'features': [
                    'Clinical aging scores',
                    'Risk stratification models',
                    'Intervention target identification',
                    'Longitudinal aging tracking'
                ]
            }
        }
    }

    return info


def demo_advanced_aging():
    """Run demonstrations of available Phase 2B components."""

    print("🧬 Advanced Aging Analysis - Phase 2B Demo")
    print("=" * 50)

    info = get_component_info()
    available = get_available_components()

    print("📊 Implementation Status:")
    print(f"   • Implemented: {info['implementation_status']['implemented']}")
    print(f"   • Planned: {info['implementation_status']['planned']}")
    print(f"   • Total: {info['implementation_status']['total']}")

    if not available:
        print("\n❌ No components available for demo")
        print("💡 Components may require additional dependencies:")
        print("   • scanpy: pip install scanpy")
        print("   • squidpy: pip install squidpy")
        print("   • scikit-learn: pip install scikit-learn")
        print("   • xgboost: pip install xgboost")
        print("   • shap: pip install shap")
        return

    print(f"\n✅ Available for demo: {len(available)}")
    for component in available:
        print(f"   • {component}")

    # Run available demos
    if SINGLE_CELL_AVAILABLE:
        print("\n" + "="*50)
        print("🧬 Single-Cell Aging Analysis Demo")
        print("="*50)
        try:
            from .single_cell_aging import demo_single_cell_aging
            demo_single_cell_aging()
        except Exception as e:
            print(f"❌ Single-cell demo failed: {e}")

    if SPATIAL_AVAILABLE:
        print("\n" + "="*50)
        print("🗺️  Spatial Aging Analysis Demo")
        print("="*50)
        try:
            from .spatial_aging import demo_spatial_aging
            demo_spatial_aging()
        except Exception as e:
            print(f"❌ Spatial demo failed: {e}")

    if AI_DISCOVERY_AVAILABLE:
        print("\n" + "="*50)
        print("🤖 AI Biomarker Discovery Demo")
        print("="*50)
        try:
            from .ai_biomarker_discovery import demo_ai_biomarker_discovery
            demo_ai_biomarker_discovery()
        except Exception as e:
            print(f"❌ AI discovery demo failed: {e}")

    if SENESCENCE_AVAILABLE:
        print("\n" + "="*50)
        print("🧬 Senescence Classification Demo")
        print("="*50)
        try:
            from .senescence_classifier import demo_senescence_classifier
            demo_senescence_classifier()
        except Exception as e:
            print(f"❌ Senescence demo failed: {e}")

    print("\n" + "="*50)
    print("🎉 Phase 2B Demo Complete!")
    print("="*50)


if __name__ == "__main__":
    # Display component information and run demos
    info = get_component_info()
    print(f"Phase {info['version']} - {info['phase']}")
    print(f"Components: {info['implementation_status']['implemented']}/{info['implementation_status']['total']} implemented")

    demo_advanced_aging()
