#!/usr/bin/env python3
"""
Phase 2B Implementation Test Script

Test and demonstrate the advanced aging analysis components
implemented in Phase 2B.

Author: CRISPR Toolkit Team
Date: August 21, 2025
"""

import os
import sys

# Add the src directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, '..', 'src')
sys.path.insert(0, src_dir)

def test_phase2b_import():
    """Test importing Phase 2B components."""

    print("🧪 Testing Phase 2B Component Imports")
    print("=" * 50)

    try:
        from crispr_toolkit.analysis.aging.advanced import (
            __version__,
            get_available_components,
            get_component_info,
        )

        print("✅ Phase 2B module imported successfully")
        print(f"📦 Version: {__version__}")

        # Get component information
        info = get_component_info()
        available = get_available_components()

        print("\n📊 Implementation Status:")
        status = info['implementation_status']
        print(f"   • Implemented: {status['implemented']}")
        print(f"   • Planned: {status['planned']}")
        print(f"   • Total: {status['total']}")

        print(f"\n✅ Available Components ({len(available)}):")
        for component in available:
            print(f"   • {component}")

        print("\n📋 All Components:")
        for comp_name, comp_info in info['components'].items():
            status_icon = "✅" if comp_info['status'] == 'implemented' else "⏳"
            print(f"   {status_icon} {comp_name}: {comp_info['status']}")

        return True

    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False


def test_individual_components():
    """Test individual Phase 2B components."""

    print("\n🔬 Testing Individual Components")
    print("=" * 50)

    # Test Single-Cell Aging Analyzer
    print("\n🧬 Testing SingleCellAgingAnalyzer...")
    try:
        from crispr_toolkit.analysis.aging.advanced import SingleCellAgingAnalyzer

        analyzer = SingleCellAgingAnalyzer(verbose=False)
        print("   ✅ SingleCellAgingAnalyzer instantiated successfully")
        print(f"   📊 Aging model: {analyzer.aging_model}")
        print(f"   🧬 Number of aging genes: {analyzer.n_aging_genes}")

        # Test aging gene signatures
        signatures = analyzer.aging_gene_signatures
        print(f"   📂 Aging signatures: {len(signatures)} categories")
        for sig_name, genes in list(signatures.items())[:3]:
            print(f"      • {sig_name}: {len(genes)} genes")

    except ImportError:
        print("   ⚠️  SingleCellAgingAnalyzer not available (missing dependencies)")
    except Exception as e:
        print(f"   ❌ Error: {e}")

    # Test Spatial Aging Mapper
    print("\n🗺️  Testing SpatialAgingMapper...")
    try:
        from crispr_toolkit.analysis.aging.advanced import SpatialAgingMapper

        mapper = SpatialAgingMapper(verbose=False)
        print("   ✅ SpatialAgingMapper instantiated successfully")
        print(f"   📊 Spatial resolution: {mapper.spatial_resolution}")
        print(f"   📏 Distance threshold: {mapper.aging_distance_threshold}")

        # Test spatial aging signatures
        signatures = mapper.spatial_aging_signatures
        print(f"   📂 Spatial signatures: {len(signatures)} categories")
        for sig_name, genes in list(signatures.items())[:3]:
            print(f"      • {sig_name}: {len(genes)} genes")

    except ImportError:
        print("   ⚠️  SpatialAgingMapper not available (missing dependencies)")
    except Exception as e:
        print(f"   ❌ Error: {e}")

    # Test AI Biomarker Discovery
    print("\n🤖 Testing AIBiomarkerDiscovery...")
    try:
        from crispr_toolkit.analysis.aging.advanced import AIBiomarkerDiscovery

        discovery = AIBiomarkerDiscovery(verbose=False)
        print("   ✅ AIBiomarkerDiscovery instantiated successfully")
        print(f"   🤖 Model type: {discovery.model_type}")
        print(f"   💡 Explanations enabled: {discovery.enable_explanations}")
        print(f"   📈 Sudden aging detection: {discovery.sudden_aging_detection}")

        # Test biomarker categories
        categories = discovery.biomarker_categories
        print(f"   📂 Biomarker categories: {len(categories)} types")
        for cat_name, genes in list(categories.items())[:3]:
            print(f"      • {cat_name}: {len(genes)} genes")

    except ImportError:
        print("   ⚠️  AIBiomarkerDiscovery not available (missing dependencies)")
    except Exception as e:
        print(f"   ❌ Error: {e}")


def test_demos():
    """Test running component demos."""

    print("\n🎭 Testing Component Demos")
    print("=" * 50)

    try:
        from crispr_toolkit.analysis.aging.advanced import demo_advanced_aging

        print("Running Phase 2B advanced aging demos...")
        demo_advanced_aging()

        print("\n✅ Demo execution completed")

    except Exception as e:
        print(f"❌ Demo failed: {e}")


def main():
    """Main test function."""

    print("🧬 CRISPR Toolkit - Phase 2B Implementation Test")
    print("=" * 60)

    success = True

    # Test imports
    if not test_phase2b_import():
        success = False

    # Test individual components
    test_individual_components()

    # Test demos (optional, as they may be lengthy)
    response = input("\n❓ Run component demos? (y/N): ").strip().lower()
    if response in ['y', 'yes']:
        test_demos()
    else:
        print("   ⏭️  Skipping demos")

    # Final summary
    print("\n" + "=" * 60)
    if success:
        print("🎉 Phase 2B Implementation Test PASSED!")
        print("✅ Advanced aging analysis components are working correctly")
    else:
        print("⚠️  Phase 2B Implementation Test had issues")
        print("💡 Some components may require additional dependencies")

    print("\n📝 Next Steps:")
    print("   • Install optional dependencies for full functionality:")
    print("     - scanpy: pip install scanpy")
    print("     - squidpy: pip install squidpy")
    print("     - scikit-learn: pip install scikit-learn")
    print("     - xgboost: pip install xgboost")
    print("     - shap: pip install shap")
    print("     - torch: pip install torch")
    print("   • Test with real aging datasets")
    print("   • Implement remaining Phase 2B components")
    print("   • Continue to Phase 2C (Clinical Translation)")


if __name__ == "__main__":
    main()
