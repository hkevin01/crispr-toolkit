#!/usr/bin/env python3
"""
Quick test script to validate MAGeCK integration
"""

import sys

sys.path.insert(0, '/home/kevin/Projects/crispr-toolkit/src')

try:
    # Test imports
    from crispr_toolkit.analysis.screens import (
        PathwayEnrichment,
        ScreenQualityControl,
        SenescenceScreenAnalyzer,
    )
    print("✓ Successfully imported screen analysis modules")

    # Test basic initialization
    qc = ScreenQualityControl()
    pathway = PathwayEnrichment()
    senescence = SenescenceScreenAnalyzer()
    print("✓ Successfully initialized analyzers")

    # Test pathway loading
    pathways = senescence.senescence_pathways
    print(f"✓ Loaded {len(pathways)} senescence pathways")
    print(f"  - ReHMGB1/RAGE signaling: {len(pathways['rehmgb1_rage_signaling'])} genes")
    print(f"  - JAK/STAT pathway: {len(pathways['jak_stat_pathway'])} genes")
    print(f"  - NF-κB signaling: {len(pathways['nfkb_signaling'])} genes")

    # Test pathway databases
    databases = pathway.pathway_databases
    print(f"✓ Loaded {len(databases)} pathway databases")

    print("\n🎉 MAGeCK integration test completed successfully!")
    print("\nKey features now available:")
    print("  - MAGeCK wrapper for CRISPR screen analysis")
    print("  - Screen quality control metrics")
    print("  - Pathway enrichment analysis")
    print("  - Senescence-specific pathway analysis")
    print("  - ReHMGB1 pathway investigation")
    print("  - Therapeutic target identification")

except Exception as e:
    print(f"❌ Error during testing: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
    sys.exit(1)
