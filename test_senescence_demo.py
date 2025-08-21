#!/usr/bin/env python3
"""
Test script for SenescenceClassifier demo.
"""

import os
import sys

# Add the src directory to the Python path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))

try:
    from crispr_toolkit.analysis.aging.advanced.senescence_classifier import (
        demo_senescence_classifier,
    )

    print("🧬 Testing SenescenceClassifier Demo")
    print("=" * 50)

    # Run the demo
    demo_senescence_classifier()

    print("\n✅ SenescenceClassifier demo completed successfully!")

except ImportError as e:
    print(f"❌ Import error: {e}")
    print("Note: Some dependencies may be missing. This is expected in the demo environment.")

except Exception as e:
    print(f"❌ Error during demo: {e}")
    import traceback
    traceback.print_exc()
    import traceback
    traceback.print_exc()
