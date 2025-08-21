"""
Aging Biomarker Analysis Example

Demonstrates how to use the comprehensive aging biomarker analysis pipeline
with focus on ReHMGB1/RAGE signaling pathways for senescence research.

This example shows integration of:
- pyaging (100+ aging clocks)
- biolearn (standardized biomarkers)
- ReHMGB1 pathway analysis
- Clinical scoring and interpretation
- Comprehensive visualization
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')

def create_example_data():
    """Create example dataset for aging biomarker analysis."""

    print("🧬 Creating example aging biomarker dataset...")

    # Sample metadata
    np.random.seed(42)
    n_samples = 50

    sample_metadata = pd.DataFrame({
        'sample_id': [f"Patient_{i:03d}" for i in range(n_samples)],
        'chronological_age': np.random.uniform(30, 80, n_samples),
        'sex': np.random.choice(['M', 'F'], n_samples),
        'treatment_group': np.random.choice(['Control', 'Treatment'], n_samples),
        'tissue_type': np.random.choice(['Blood', 'Tissue'], n_samples)
    })

    # Simulate gene expression data (focusing on aging/senescence genes)
    aging_genes = [
        'HMGB1', 'AGER', 'STAT1', 'STAT3', 'NFKB1', 'RELA',
        'IL6', 'TNF', 'IL1B', 'CDKN1A', 'CDKN2A', 'TP53',
        'RB1', 'ATM', 'SIRT1', 'FOXO1', 'FOXO3', 'TERT',
        'SOD1', 'SOD2', 'CAT', 'GPX1', 'NRF2', 'KEAP1'
    ]

    # Create expression matrix (log2 normalized counts)
    expression_data = pd.DataFrame(
        np.random.normal(5, 2, (n_samples, len(aging_genes))),
        index=sample_metadata['sample_id'],
        columns=aging_genes
    )

    # Add age-correlated expression patterns
    for i, sample in sample_metadata.iterrows():
        age_factor = (sample['chronological_age'] - 30) / 50  # Normalize to 0-1

        # Age-up genes (increase with age)
        age_up_genes = ['HMGB1', 'IL6', 'TNF', 'CDKN1A', 'CDKN2A', 'TP53']
        for gene in age_up_genes:
            if gene in expression_data.columns:
                expression_data.loc[sample['sample_id'], gene] += age_factor * 2

        # Age-down genes (decrease with age)
        age_down_genes = ['SIRT1', 'TERT', 'SOD1', 'SOD2']
        for gene in age_down_genes:
            if gene in expression_data.columns:
                expression_data.loc[sample['sample_id'], gene] -= age_factor * 1.5

    # Ensure non-negative values
    expression_data = expression_data.clip(lower=0)

    print(f"   ✅ Created dataset with {n_samples} samples and {len(aging_genes)} genes")
    print(f"   📊 Age range: {sample_metadata['chronological_age'].min():.1f} - "
          f"{sample_metadata['chronological_age'].max():.1f} years")

    return sample_metadata, expression_data


def run_aging_biomarker_analysis():
    """Run comprehensive aging biomarker analysis example."""

    print("🔬 CRISPR Toolkit - Aging Biomarker Analysis Example")
    print("=" * 60)

    # Create example data
    sample_metadata, expression_data = create_example_data()

    # Try to import aging biomarker modules
    print("\n📦 Loading aging biomarker analysis modules...")
    try:
        from crispr_toolkit.analysis.aging.biomarkers import (
            AgingBiomarkerAnalyzer,
            AgingBiomarkerVisualizer,
            BiolearnAnalyzer,
            ClinicalAgingScorer,
            PyAgingAnalyzer,
            ReHMGB1BiomarkerAnalyzer,
        )
        print("   ✅ All aging biomarker modules loaded successfully")
    except ImportError as e:
        print(f"   ❌ Module import failed: {e}")
        print("   💡 This is expected in demo mode - modules would be available after installation")
        print("   📦 Install dependencies:")
        print("      pip install pyaging biolearn matplotlib seaborn plotly")
        return

    # Initialize analyzers
    print("\n🛠️  Initializing aging biomarker analyzers...")

    try:
        # Main aging biomarker analyzer (combines pyaging + biolearn)
        aging_analyzer = AgingBiomarkerAnalyzer(verbose=True)
        print("   ✅ Main aging analyzer initialized")

        # ReHMGB1-specific analyzer
        rehmgb1_analyzer = ReHMGB1BiomarkerAnalyzer(verbose=True)
        print("   ✅ ReHMGB1 pathway analyzer initialized")

        # Clinical scoring
        clinical_scorer = ClinicalAgingScorer(verbose=True)
        print("   ✅ Clinical aging scorer initialized")

        # Visualization
        visualizer = AgingBiomarkerVisualizer(style='scientific', verbose=True)
        print("   ✅ Biomarker visualizer initialized")

    except Exception as e:
        print(f"   ❌ Analyzer initialization failed: {e}")
        print("   💡 This would work with proper pyaging/biolearn installation")
        return

    # Analyze aging biomarkers
    print("\n⚗️  Running aging biomarker analysis...")

    try:
        # 1. Comprehensive aging analysis
        print("   🔍 Running comprehensive aging clock analysis...")
        aging_results = aging_analyzer.analyze_aging_biomarkers(
            expression_data=expression_data,
            metadata=sample_metadata,
            chronological_age_col='chronological_age',
            selected_clocks=['Horvath2013', 'PhenoAge', 'GrimAge', 'DunedinPACE']
        )
        print(f"      ✅ Analyzed {len(aging_results.get('predictions', []))} predictions")

        # 2. ReHMGB1 pathway analysis
        print("   🧬 Running ReHMGB1 pathway analysis...")
        rehmgb1_results = rehmgb1_analyzer.analyze_rehmgb1_biomarkers(
            expression_data=expression_data,
            metadata=sample_metadata
        )
        print(f"      ✅ Analyzed {len(rehmgb1_results.get('pathway_scores', {}))} pathways")

        # 3. Clinical scoring
        print("   🏥 Computing clinical aging scores...")
        clinical_results = {}
        for idx, sample in sample_metadata.iterrows():
            sample_id = sample['sample_id']

            # Get aging predictions for this sample
            sample_predictions = [
                pred for pred in aging_results.get('predictions', [])
                if pred.get('sample_id') == sample_id
            ]

            if sample_predictions:
                # Create aging profile
                aging_profile = clinical_scorer.create_aging_profile(
                    chronological_age=sample['chronological_age'],
                    aging_predictions=sample_predictions,
                    pathway_scores=rehmgb1_results.get('pathway_scores', {}),
                    expression_data=expression_data.loc[sample_id].to_dict()
                )

                # Calculate risk scores
                risk_scores = clinical_scorer.calculate_risk_scores(aging_profile)

                clinical_results[sample_id] = {
                    'aging_profile': aging_profile,
                    'risk_scores': risk_scores
                }

        print(f"      ✅ Computed clinical scores for {len(clinical_results)} samples")

    except Exception as e:
        print(f"   ❌ Analysis failed: {e}")
        print("   💡 This would work with proper data and installed dependencies")
        return

    # Create visualizations
    print("\n📊 Creating aging biomarker visualizations...")

    try:
        output_dir = Path('./aging_analysis_output')
        output_dir.mkdir(exist_ok=True)

        # Prepare visualization data
        vis_data = {
            'aging_predictions': pd.DataFrame(aging_results.get('predictions', [])),
            'chronological_age': sample_metadata.set_index('sample_id')['chronological_age'],
            'pathway_scores': rehmgb1_results.get('pathway_scores', {}),
            'aging_profile': list(clinical_results.values())[0]['aging_profile'] if clinical_results else None,
            'risk_scores': list(clinical_results.values())[0]['risk_scores'] if clinical_results else None
        }

        # Create comprehensive visualization report
        plot_files = visualizer.create_aging_report_visualizations(
            aging_data=vis_data,
            output_dir=str(output_dir),
            prefix='example_aging_analysis'
        )

        print(f"   ✅ Created {len(plot_files)} visualization files:")
        for plot_type, file_path in plot_files.items():
            print(f"      📈 {plot_type}: {file_path}")

    except Exception as e:
        print(f"   ⚠️  Visualization creation failed: {e}")
        print("   💡 Install visualization dependencies: matplotlib, seaborn, plotly")

    # Generate summary report
    print("\n📋 Generating aging biomarker analysis summary...")

    try:
        summary_stats = {}

        if aging_results.get('predictions'):
            predictions_df = pd.DataFrame(aging_results['predictions'])

            # Calculate summary statistics
            if 'predicted_age' in predictions_df.columns:
                summary_stats['avg_biological_age'] = predictions_df['predicted_age'].mean()
                summary_stats['biological_age_std'] = predictions_df['predicted_age'].std()

            # Age acceleration (if chronological age available)
            if 'chronological_age' in predictions_df.columns:
                age_accel = predictions_df['predicted_age'] - predictions_df['chronological_age']
                summary_stats['avg_age_acceleration'] = age_accel.mean()
                summary_stats['age_acceleration_std'] = age_accel.std()

        # ReHMGB1 pathway summary
        if rehmgb1_results.get('pathway_scores'):
            pathway_scores = rehmgb1_results['pathway_scores']
            summary_stats['avg_hmgb1_score'] = pathway_scores.get('hmgb1_core', 0)
            summary_stats['avg_rage_score'] = pathway_scores.get('rage_signaling', 0)
            summary_stats['avg_inflammation_score'] = pathway_scores.get('sasp_factors', 0)

        # Clinical risk summary
        if clinical_results:
            mortality_risks = []
            senescence_burdens = []

            for result in clinical_results.values():
                if result.get('aging_profile'):
                    profile = result['aging_profile']
                    if hasattr(profile, 'mortality_risk'):
                        mortality_risks.append(profile.mortality_risk)
                    if hasattr(profile, 'senescence_burden'):
                        senescence_burdens.append(profile.senescence_burden)

            if mortality_risks:
                summary_stats['avg_mortality_risk'] = np.mean(mortality_risks)
            if senescence_burdens:
                summary_stats['avg_senescence_burden'] = np.mean(senescence_burdens)

        # Print summary
        print("\n🎯 AGING BIOMARKER ANALYSIS SUMMARY")
        print("-" * 50)

        for stat_name, stat_value in summary_stats.items():
            formatted_name = stat_name.replace('_', ' ').title()
            if isinstance(stat_value, float):
                print(f"   {formatted_name}: {stat_value:.2f}")
            else:
                print(f"   {formatted_name}: {stat_value}")

    except Exception as e:
        print(f"   ⚠️  Summary generation failed: {e}")

    print("\n🎉 Aging biomarker analysis example completed!")

    print("\n💡 Key Features Demonstrated:")
    print("   ✅ Multi-clock aging analysis (pyaging integration)")
    print("   ✅ Standardized biomarkers (biolearn integration)")
    print("   ✅ ReHMGB1/RAGE pathway-specific analysis")
    print("   ✅ Clinical risk scoring and interpretation")
    print("   ✅ Comprehensive visualization suite")
    print("   ✅ Senescence burden quantification")

    print("\n🔬 Research Applications:")
    print("   • Aging biomarker discovery and validation")
    print("   • ReHMGB1 senescence pathway analysis")
    print("   • Clinical aging assessment tools")
    print("   • Longitudinal aging studies")
    print("   • Therapeutic intervention monitoring")
    print("   • Personalized aging medicine")

    print("\n📚 Next Steps:")
    print("   1. Install required dependencies (pyaging, biolearn)")
    print("   2. Prepare your own expression/methylation data")
    print("   3. Customize aging clocks for your research question")
    print("   4. Integrate with CRISPR screen analysis results")
    print("   5. Develop custom biomarker panels")


if __name__ == "__main__":
    run_aging_biomarker_analysis()
if __name__ == "__main__":
    run_aging_biomarker_analysis()
