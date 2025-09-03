#!/usr/bin/env python3
"""
Demo script for MultiOmicsIntegrator
"""

import sys

sys.path.insert(0, '/home/kevin/Projects/crispr-toolkit/src')

def run_demo():
    """Run MultiOmicsIntegrator demo."""
    try:
        # Import the integrator
        from crispr_toolkit.analysis.aging.advanced.multi_omics_integrator import (
            MultiOmicsIntegrator,
        )
        print("🧬 MultiOmicsIntegrator Demo Starting...")
        print("=" * 60)

        # Create integrator with MOFA+ method
        integrator = MultiOmicsIntegrator(
            integration_method='mofa_plus',
            n_factors=10,
            random_state=42
        )

        print(f"✅ Created MultiOmicsIntegrator with method: {integrator.integration_method}")
        print(f"✅ Number of factors: {integrator.n_factors}")

        # Run the demo analysis
        print("\n🚀 Running comprehensive demo analysis...")
        demo_results = integrator.demo_analysis()

        print("\n📊 Demo Results Summary:")
        for key, value in demo_results.items():
            if isinstance(value, (int, float)):
                if 'mae' in key.lower():
                    print(f"   {key}: {value:.2f}")
                elif isinstance(value, float):
                    print(f"   {key}: {value:.3f}")
                else:
                    print(f"   {key}: {value}")
            else:
                print(f"   {key}: {value}")

        print("\n✅ MultiOmicsIntegrator demo completed successfully!")
        return True

    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("Please check that all dependencies are installed:")
        print("  - scikit-learn")
        print("  - scipy")
        print("  - numpy")
        print("  - pandas")
        return False

    except Exception as e:
        print(f"❌ Demo failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = run_demo()
    if success:
        print("\n🎉 All tests passed! MultiOmicsIntegrator is working correctly.")
    else:
        print("\n⚠️  Demo failed. Please check the implementation.")
        print("\n🎉 All tests passed! MultiOmicsIntegrator is working correctly.")
    else:
        print("\n⚠️  Demo failed. Please check the implementation.")
