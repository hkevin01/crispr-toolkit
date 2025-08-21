#!/usr/bin/env python3
"""
Simple test for MultiOmicsIntegrator
"""

import sys

sys.path.insert(0, '/home/kevin/Projects/crispr-toolkit/src')

# Test import
try:
    from crispr_toolkit.analysis.aging.advanced.multi_omics_integrator import (
        MultiOmicsIntegrator,
    )
    print("✅ MultiOmicsIntegrator imported successfully")

    # Test basic instantiation
    integrator = MultiOmicsIntegrator(integration_method='mofa_plus')
    print("✅ MultiOmicsIntegrator instantiated successfully")

except Exception as e:
    print(f"❌ Error: {e}")

# Test advanced aging module
try:
    from crispr_toolkit.analysis.aging.advanced import get_available_components
    components = get_available_components()
    print(f"✅ Available components: {components}")

except Exception as e:
    print(f"❌ Module import error: {e}")

print("🎉 Test completed!")
    print(f"❌ Module import error: {e}")

print("🎉 Test completed!")
