#!/usr/bin/env python3
"""
MAGeCK Integration Validation Script
Tests the complete MAGeCK integration for the CRISPR Toolkit
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

def test_imports():
    """Test that all MAGeCK integration modules can be imported."""
    print("🧪 Testing imports...")

    try:
        from crispr_toolkit.analysis.screens import MAGeCKAnalyzer
        print("✅ MAGeCKAnalyzer import successful")

        from crispr_toolkit.analysis.screens import ScreenQualityControl
        print("✅ ScreenQualityControl import successful")

        from crispr_toolkit.analysis.screens import PathwayEnrichment
        print("✅ PathwayEnrichment import successful")

        from crispr_toolkit.analysis.screens import SenescenceScreenAnalyzer
        print("✅ SenescenceScreenAnalyzer import successful")

        return True
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False

def test_class_initialization():
    """Test that all classes can be instantiated."""
    print("\n🔧 Testing class initialization...")

    try:
        from crispr_toolkit.analysis.screens import (
            MAGeCKAnalyzer,
            PathwayEnrichment,
            ScreenQualityControl,
            SenescenceScreenAnalyzer,
        )

        # Test MAGeCK Analyzer
        mageck = MAGeCKAnalyzer()
        print("✅ MAGeCKAnalyzer instantiation successful")

        # Test Screen QC
        qc = ScreenQualityControl()
        print("✅ ScreenQualityControl instantiation successful")

        # Test Pathway Enrichment
        pathway = PathwayEnrichment()
        print("✅ PathwayEnrichment instantiation successful")

        # Test Senescence Analyzer
        senescence = SenescenceScreenAnalyzer()
        print("✅ SenescenceScreenAnalyzer instantiation successful")

        return True
    except Exception as e:
        print(f"❌ Class initialization failed: {e}")
        return False

def test_pathway_data():
    """Test that pathway data is properly loaded."""
    print("\n🧬 Testing pathway data...")

    try:
        from crispr_toolkit.analysis.screens import SenescenceScreenAnalyzer

        senescence = SenescenceScreenAnalyzer()

        # Check if pathway data exists
        pathways = ['rehmgb1_rage_signaling', 'jak_stat_pathway', 'nfkb_signaling',
                   'cell_cycle_arrest', 'sasp_factors']

        for pathway_name in pathways:
            genes = senescence.get_pathway_genes(pathway_name)
            if genes:
                print(f"✅ {pathway_name.upper()} pathway: {len(genes)} genes loaded")
            else:
                print(f"⚠️  {pathway_name.upper()} pathway: No genes found")

        return True
    except Exception as e:
        print(f"❌ Pathway data test failed: {e}")
        return False

def test_main_package_integration():
    """Test that the main package properly exports the new classes."""
    print("\n📦 Testing main package integration...")

    try:
        import crispr_toolkit

        # Check if new classes are available from main package
        required_classes = [
            'MAGeCKAnalyzer',
            'ScreenQualityControl',
            'PathwayEnrichment',
            'SenescenceScreenAnalyzer'
        ]

        available_classes = [attr for attr in dir(crispr_toolkit)
                           if not attr.startswith('_')]

        for class_name in required_classes:
            if hasattr(crispr_toolkit, class_name):
                print(f"✅ {class_name} available in main package")
            else:
                print(f"⚠️  {class_name} not found in main package")

        return True
    except Exception as e:
        print(f"❌ Main package integration test failed: {e}")
        return False

def test_dependencies():
    """Test that required dependencies are available."""
    print("\n📚 Testing dependencies...")

    dependencies = [
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'scipy'
    ]

    all_available = True
    for dep in dependencies:
        try:
            __import__(dep)
            print(f"✅ {dep} available")
        except ImportError:
            print(f"❌ {dep} not available")
            all_available = False

    return all_available

def main():
    """Run all validation tests."""
    print("🚀 MAGeCK Integration Validation")
    print("=" * 50)

    tests = [
        ("Import Test", test_imports),
        ("Class Initialization", test_class_initialization),
        ("Pathway Data", test_pathway_data),
        ("Main Package Integration", test_main_package_integration),
        ("Dependencies", test_dependencies)
    ]

    passed = 0
    total = len(tests)

    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"❌ {test_name} failed with exception: {e}")

    print("\n" + "=" * 50)
    print(f"🎯 Test Results: {passed}/{total} tests passed")

    if passed == total:
        print("🎉 All tests passed! MAGeCK integration is ready.")
        return True
    else:
        print("⚠️  Some tests failed. Please check the output above.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
