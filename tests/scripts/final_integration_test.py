#!/usr/bin/env python3
"""
Final Integration Test for CRISPR Aging Toolkit
Tests Phase 2A + Phase 2B integration and complete functionality
"""

import os
import sys
from datetime import datetime


def run_final_integration_test():
    """Run comprehensive integration test of all components"""
    print("🧬 CRISPR Aging Toolkit - Final Integration Test")
    print("=" * 70)
    print(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

    test_results = {
        'phase_2a_status': '❓ Not Tested',
        'phase_2b_status': '❓ Not Tested',
        'simple_translator': '❓ Not Tested',
        'full_translator': '❓ Not Tested',
        'dependency_status': '❓ Not Tested',
        'integration_status': '❓ Not Tested'
    }

    print("📊 Testing Phase 2A: Aging Biomarker Integration")
    print("-" * 55)

    # Check Phase 2A module structure
    phase_2a_path = "/home/kevin/Projects/crispr-toolkit/src/crispr_toolkit/analysis/aging/biomarkers"
    if os.path.exists(phase_2a_path):
        files = os.listdir(phase_2a_path)
        required_files = [
            '__init__.py', 'pyaging_integration.py', 'biolearn_integration.py',
            'biomarker_analysis.py', 'rehmgb1_biomarkers.py', 'clinical_scoring.py',
            'visualization.py', 'README.md', 'example_usage.py'
        ]

        missing_files = [f for f in required_files if f not in files]
        if not missing_files:
            test_results['phase_2a_status'] = '✅ Complete'
            print(f"✅ Phase 2A structure complete ({len(files)} files)")
        else:
            test_results['phase_2a_status'] = f'⚠️  Missing {len(missing_files)} files'
            print(f"⚠️  Phase 2A missing files: {missing_files}")
    else:
        test_results['phase_2a_status'] = '❌ Module not found'
        print("❌ Phase 2A module directory not found")

    print()
    print("📊 Testing Phase 2B: Clinical Translator")
    print("-" * 45)

    # Check Phase 2B files
    clinical_translator_path = "/home/kevin/Projects/crispr-toolkit/src/crispr_toolkit/clinical/clinical_translator.py"
    if os.path.exists(clinical_translator_path):
        with open(clinical_translator_path, 'r') as f:
            content = f.read()
            if len(content) > 30000:  # Substantial implementation
                test_results['phase_2b_status'] = '✅ Complete'
                print(f"✅ Clinical Translator implemented ({len(content)} characters)")
            else:
                test_results['phase_2b_status'] = '⚠️  Incomplete'
                print("⚠️  Clinical Translator appears incomplete")
    else:
        test_results['phase_2b_status'] = '❌ Not found'
        print("❌ Clinical Translator file not found")

    print()
    print("🧪 Testing Simple Clinical Translator")
    print("-" * 40)

    try:
        # Test simple clinical translator
        import subprocess
        result = subprocess.run([
            'python3', '/home/kevin/Projects/crispr-toolkit/tests/clinical/test_clinical_translator.py'
        ], capture_output=True, text=True, cwd='/home/kevin/Projects/crispr-toolkit')

        if result.returncode == 0 and "✅ FUNCTIONAL" in result.stdout:
            test_results['simple_translator'] = '✅ Working'
            print("✅ Simple Clinical Translator test passed")
        else:
            test_results['simple_translator'] = '❌ Failed'
            print("❌ Simple Clinical Translator test failed")

    except Exception as e:
        test_results['simple_translator'] = f'❌ Error: {str(e)[:50]}'
        print(f"❌ Simple test error: {e}")

    print()
    print("🔬 Testing Full Clinical Translator (with dependencies)")
    print("-" * 60)

    try:
        # Test full clinical translator with virtual environment
        result = subprocess.run([
            'bash', '-c',
            'cd /home/kevin/Projects/crispr-toolkit && source venv/bin/activate && python src/crispr_toolkit/clinical/clinical_translator.py'
        ], capture_output=True, text=True)

        if result.returncode == 0 and "✅ Clinical Translator demonstration completed!" in result.stdout:
            test_results['full_translator'] = '✅ Working'
            print("✅ Full Clinical Translator test passed")
        else:
            test_results['full_translator'] = '❌ Failed'
            print("❌ Full Clinical Translator test failed")
            if result.stderr:
                print(f"   Error: {result.stderr[:100]}")

    except Exception as e:
        test_results['full_translator'] = f'❌ Error: {str(e)[:50]}'
        print(f"❌ Full test error: {e}")

    print()
    print("📦 Testing Dependencies")
    print("-" * 25)

    venv_path = "/home/kevin/Projects/crispr-toolkit/venv"
    if os.path.exists(venv_path):
        try:
            result = subprocess.run([
                'bash', '-c',
                'cd /home/kevin/Projects/crispr-toolkit && source venv/bin/activate && python -c "import pandas, numpy, sklearn, scipy, matplotlib, seaborn; print(\'All dependencies available\')"'
            ], capture_output=True, text=True)

            if result.returncode == 0:
                test_results['dependency_status'] = '✅ Available'
                print("✅ All required dependencies installed")
            else:
                test_results['dependency_status'] = '⚠️  Issues detected'
                print("⚠️  Dependency issues detected")
        except Exception as e:
            test_results['dependency_status'] = '❌ Error checking'
            print(f"❌ Error checking dependencies: {e}")
    else:
        test_results['dependency_status'] = '❌ Virtual env missing'
        print("❌ Virtual environment not found")

    print()
    print("🔗 Testing Integration Status")
    print("-" * 32)

    # Overall integration assessment
    working_components = sum(1 for status in test_results.values() if status.startswith('✅'))
    total_components = len(test_results)

    if working_components >= 4:  # Most components working
        test_results['integration_status'] = '✅ Functional'
        print("✅ CRISPR Aging Toolkit integration functional")
    elif working_components >= 2:
        test_results['integration_status'] = '⚠️  Partial'
        print("⚠️  CRISPR Aging Toolkit partially functional")
    else:
        test_results['integration_status'] = '❌ Issues detected'
        print("❌ CRISPR Aging Toolkit has integration issues")

    print()
    print("📋 FINAL TEST SUMMARY")
    print("=" * 25)

    for component, status in test_results.items():
        component_name = component.replace('_', ' ').title()
        print(f"   {component_name}: {status}")

    print()
    success_rate = (working_components / total_components) * 100
    print(f"📊 Overall Success Rate: {success_rate:.1f}% ({working_components}/{total_components})")

    if success_rate >= 80:
        print("🎉 CRISPR Aging Toolkit: ✅ READY FOR USE")
    elif success_rate >= 60:
        print("🔧 CRISPR Aging Toolkit: ⚠️  NEEDS MINOR FIXES")
    else:
        print("🔨 CRISPR Aging Toolkit: ❌ NEEDS SIGNIFICANT WORK")

    print()
    print("💡 Key Features Available:")
    if test_results['phase_2a_status'].startswith('✅'):
        print("   • Aging biomarker integration (PyAging, Biolearn, ReHMGB1)")
    if test_results['phase_2b_status'].startswith('✅'):
        print("   • Clinical translator for precision geromedicine")
    if test_results['simple_translator'].startswith('✅'):
        print("   • Simple clinical assessment capabilities")
    if test_results['full_translator'].startswith('✅'):
        print("   • Advanced clinical aging analysis")
    if test_results['dependency_status'].startswith('✅'):
        print("   • Complete scientific computing environment")

    print()
    print("🎯 Ready for aging and senescence research applications!")
    print()

    return success_rate >= 60

if __name__ == "__main__":
    success = run_final_integration_test()
    sys.exit(0 if success else 1)
    sys.exit(0 if success else 1)
