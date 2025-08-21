#!/usr/bin/env python3

import sys

sys.path.insert(0, '/home/kevin/Projects/crispr-toolkit/src')

# Quick test
print("Testing MultiOmicsIntegrator...")

try:
    exec(open('/home/kevin/Projects/crispr-toolkit/src/crispr_toolkit/analysis/aging/advanced/multi_omics_integrator.py').read())
    print("✅ MultiOmicsIntegrator code is valid")

    # Create instance
    integrator = MultiOmicsIntegrator()
    print("✅ MultiOmicsIntegrator created successfully")

    # Check methods
    print(f"✅ Integration method: {integrator.integration_method}")
    print(f"✅ Number of factors: {integrator.n_factors}")

except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("✅ Test completed!")
print("✅ Test completed!")
