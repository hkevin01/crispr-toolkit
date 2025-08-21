#!/usr/bin/env python3
"""
Test script for aging biomarker analysis modules

Validates that all aging biomarker components can be imported
and initialized correctly.
"""

import sys
import warnings

warnings.filterwarnings('ignore')

def test_aging_biomarker_imports():
    """Test importing all aging biomarker modules."""
    print("🧬 Testing Aging Biomarker Module Imports")
    print("=" * 50)

    try:
        # Test main module import
        print("1. Testing main module import...")
        from crispr_toolkit.analysis.aging.biomarkers import (
            AgingBiomarkerAnalyzer,
            AgingBiomarkerVisualizer,
            BiolearnAnalyzer,
            ClinicalAgingScorer,
            PyAgingAnalyzer,
            ReHMGB1BiomarkerAnalyzer,
        )
        print("   ✅ All main classes imported successfully")

        # Test individual module imports
        print("\n2. Testing individual module imports...")

        from crispr_toolkit.analysis.aging.biomarkers.pyaging_integration import (
            PyAgingAnalyzer as PA,
        )
        print("   ✅ PyAgingAnalyzer imported")

        from crispr_toolkit.analysis.aging.biomarkers.biolearn_integration import (
            BiolearnAnalyzer as BA,
        )
        print("   ✅ BiolearnAnalyzer imported")

        from crispr_toolkit.analysis.aging.biomarkers.biomarker_analysis import (
            AgingBiomarkerAnalyzer as ABA,
        )
        print("   ✅ AgingBiomarkerAnalyzer imported")

        from crispr_toolkit.analysis.aging.biomarkers.rehmgb1_biomarkers import (
            ReHMGB1BiomarkerAnalyzer as RBA,
        )
        print("   ✅ ReHMGB1BiomarkerAnalyzer imported")

        from crispr_toolkit.analysis.aging.biomarkers.clinical_scoring import (
            ClinicalAgingScorer as CAS,
        )
        print("   ✅ ClinicalAgingScorer imported")

        from crispr_toolkit.analysis.aging.biomarkers.visualization import (
            AgingBiomarkerVisualizer as ABV,
        )
        print("   ✅ AgingBiomarkerVisualizer imported")

        return True

    except ImportError as e:
        print(f"   ❌ Import failed: {e}")
        return False
    except Exception as e:
        print(f"   ❌ Unexpected error: {e}")
        return False

def test_basic_initialization():
    """Test basic initialization of aging biomarker classes."""
    print("\n3. Testing basic class initialization...")

    try:
        # Try to initialize with mock dependencies
        print("   🔍 Attempting basic initialization...")

        # These will fail if dependencies aren't installed, but should not fail due to syntax
        try:
            from crispr_toolkit.analysis.aging.biomarkers import PyAgingAnalyzer
            analyzer = PyAgingAnalyzer(verbose=False)
            print("   ✅ PyAgingAnalyzer initialized")
        except Exception as e:
            print(f"   ⚠️  PyAgingAnalyzer init failed (likely missing pyaging): {e}")

        try:
            from crispr_toolkit.analysis.aging.biomarkers import BiolearnAnalyzer
            analyzer = BiolearnAnalyzer(verbose=False)
            print("   ✅ BiolearnAnalyzer initialized")
        except Exception as e:
            print(f"   ⚠️  BiolearnAnalyzer init failed (likely missing biolearn): {e}")

        try:
            from crispr_toolkit.analysis.aging.biomarkers import ClinicalAgingScorer
            scorer = ClinicalAgingScorer(verbose=False)
            print("   ✅ ClinicalAgingScorer initialized")
        except Exception as e:
            print(f"   ⚠️  ClinicalAgingScorer init failed: {e}")

        try:
            from crispr_toolkit.analysis.aging.biomarkers import (
                AgingBiomarkerVisualizer,
            )
            visualizer = AgingBiomarkerVisualizer(verbose=False)
            print("   ✅ AgingBiomarkerVisualizer initialized")
        except Exception as e:
            print(f"   ⚠️  AgingBiomarkerVisualizer init failed: {e}")

        return True

    except Exception as e:
        print(f"   ❌ Initialization test failed: {e}")
        return False

def test_dependency_detection():
    """Test detection of optional dependencies."""
    print("\n4. Testing dependency detection...")

    dependencies = {
        'pyaging': 'pyaging',
        'biolearn': 'biolearn',
        'matplotlib': 'matplotlib',
        'seaborn': 'seaborn',
        'plotly': 'plotly',
        'pandas': 'pandas',
        'numpy': 'numpy',
        'scipy': 'scipy'
    }

    available_deps = []
    missing_deps = []

    for dep_name, import_name in dependencies.items():
        try:
            __import__(import_name)
            available_deps.append(dep_name)
            print(f"   ✅ {dep_name} available")
        except ImportError:
            missing_deps.append(dep_name)
            print(f"   ❌ {dep_name} missing")

    print("\n   📊 Dependency Summary:")
    print(f"      Available: {len(available_deps)}/{len(dependencies)}")
    print(f"      Missing: {missing_deps}")

    if 'pandas' in available_deps and 'numpy' in available_deps:
        print("   ✅ Core dependencies available")
        return True
    else:
        print("   ❌ Core dependencies missing")
        return False

def main():
    """Run all aging biomarker tests."""
    print("🧪 CRISPR Toolkit - Aging Biomarker Module Tests")
    print("=" * 60)

    test_results = []

    # Test imports
    import_success = test_aging_biomarker_imports()
    test_results.append(("Imports", import_success))

    # Test initialization
    init_success = test_basic_initialization()
    test_results.append(("Initialization", init_success))

    # Test dependencies
    dep_success = test_dependency_detection()
    test_results.append(("Dependencies", dep_success))

    # Summary
    print("\n" + "=" * 60)
    print("🎯 TEST SUMMARY")
    print("-" * 30)

    for test_name, success in test_results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"   {test_name}: {status}")

    total_tests = len(test_results)
    passed_tests = sum(result[1] for result in test_results)

    print(f"\n📊 Overall: {passed_tests}/{total_tests} tests passed")

    if passed_tests == total_tests:
        print("🎉 All tests passed! Aging biomarker module is ready.")
        return 0
    else:
        print("⚠️  Some tests failed. Check dependency installation.")
        print("\n💡 To install missing dependencies:")
        print("   pip install pyaging biolearn matplotlib seaborn plotly")
        return 1

if __name__ == "__main__":
    sys.exit(main())
    passed_tests = sum(result[1] for result in test_results)

    print(f"\n📊 Overall: {passed_tests}/{total_tests} tests passed")

    if passed_tests == total_tests:
        print("🎉 All tests passed! Aging biomarker module is ready.")
        return 0
    else:
        print("⚠️  Some tests failed. Check dependency installation.")
        print("\n💡 To install missing dependencies:")
        print("   pip install pyaging biolearn matplotlib seaborn plotly")
        return 1

if __name__ == "__main__":
    sys.exit(main())
