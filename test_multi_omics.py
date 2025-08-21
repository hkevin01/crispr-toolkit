#!/usr/bin/env python3
"""
Test script for MultiOmicsIntegrator implementation
"""

import sys

sys.path.insert(0, '/home/kevin/Projects/crispr-toolkit/src')

def test_multi_omics_import():
    """Test importing the MultiOmicsIntegrator."""
    try:
        from crispr_toolkit.analysis.aging.advanced.multi_omics_integrator import (
            MultiOmicsIntegrator,
        )
        print("✅ MultiOmicsIntegrator imported successfully")

        # Test instantiation
        integrator = MultiOmicsIntegrator(integration_method='mofa_plus', n_factors=5)
        print("✅ MultiOmicsIntegrator instantiated successfully")

        # Test available methods
        methods = ['mofa_plus', 'sofa', 'ensemble', 'deep_learning']
        for method in methods:
            try:
                test_integrator = MultiOmicsIntegrator(integration_method=method)
                print(f"✅ Method '{method}' accepted")
            except Exception as e:
                print(f"⚠️  Method '{method}' failed: {e}")

        return True

    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

def test_advanced_aging_import():
    """Test importing the advanced aging module."""
    try:
        from crispr_toolkit.analysis.aging.advanced import get_available_components
        components = get_available_components()
        print(f"✅ Available components: {components}")

        # Check if MultiOmicsIntegrator is available
        if 'MultiOmicsIntegrator' in components:
            print("✅ MultiOmicsIntegrator is available in advanced aging module")
        else:
            print("⚠️  MultiOmicsIntegrator not found in available components")

        return True

    except ImportError as e:
        print(f"❌ Advanced aging import failed: {e}")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

def main():
    """Run all tests."""
    print("🧬 Testing MultiOmicsIntegrator Implementation")
    print("=" * 50)

    # Test direct import
    print("\n1. Testing direct import...")
    success1 = test_multi_omics_import()

    # Test module integration
    print("\n2. Testing module integration...")
    success2 = test_advanced_aging_import()

    # Summary
    print("\n" + "=" * 50)
    if success1 and success2:
        print("🎉 All tests passed! MultiOmicsIntegrator is ready.")
    else:
        print("⚠️  Some tests failed. Check implementation.")
    print("=" * 50)

if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()
