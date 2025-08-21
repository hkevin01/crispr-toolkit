"""
MAGeCK Integration Example

This example demonstrates how to use the integrated MAGeCK functionality
for analyzing CRISPR screens, with specific focus on ReHMGB1 pathway analysis
and senescence-related gene discovery.

Example Usage:
    python examples/mageck_screen_analysis.py
"""

import logging

import numpy as np
import pandas as pd

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Import our CRISPR toolkit modules
from crispr_toolkit.analysis.screens import (
    MAGeCKAnalyzer,
    PathwayEnrichment,
    ScreenQualityControl,
    SenescenceScreenAnalyzer,
)


def create_example_screen_data():
    """Create example CRISPR screen data for demonstration."""
    np.random.seed(42)

    # Define gene sets
    senescence_genes = [
        'HMGB1', 'AGER', 'JAK1', 'JAK2', 'STAT1', 'STAT3',
        'NFKB1', 'RELA', 'CDKN1A', 'CDKN2A', 'TP53', 'RB1',
        'IL6', 'IL1B', 'TNF', 'CXCL8', 'CCL2', 'ATM', 'ATR'
    ]

    control_genes = [
        f'CTRL_GENE_{i}' for i in range(100)
    ]

    all_genes = senescence_genes + control_genes

    # Create mock screen results
    n_genes = len(all_genes)

    # Simulate log2 fold changes with some senescence genes having stronger effects
    log2fc = np.random.normal(0, 0.5, n_genes)

    # Make some senescence genes more significant
    for i, gene in enumerate(all_genes):
        if gene in senescence_genes[:10]:  # Top 10 senescence genes
            log2fc[i] = np.random.normal(-1.5, 0.3)  # Stronger negative effect
        elif gene in senescence_genes[10:]:
            log2fc[i] = np.random.normal(-0.8, 0.4)  # Moderate effect

    # Simulate p-values (more significant for senescence genes)
    p_values = np.random.beta(2, 20, n_genes)
    for i, gene in enumerate(all_genes):
        if gene in senescence_genes:
            p_values[i] = np.random.beta(1, 50)  # More significant

    # Calculate FDR (simplified)
    fdr = p_values * n_genes / (np.arange(n_genes) + 1)
    fdr = np.minimum(fdr, 1.0)

    # Create DataFrame
    screen_data = pd.DataFrame({
        'gene': all_genes,
        'log2fc': log2fc,
        'pvalue': p_values,
        'fdr': fdr,
        'rank': np.arange(n_genes) + 1
    })

    # Sort by p-value for ranking
    screen_data = screen_data.sort_values('pvalue').reset_index(drop=True)
    screen_data['rank'] = np.arange(len(screen_data)) + 1

    return screen_data


def create_example_count_matrix():
    """Create example count matrix for QC analysis."""
    np.random.seed(42)

    genes = [f'GENE_{i}' for i in range(1000)]
    samples = ['Control_1', 'Control_2', 'Treatment_1', 'Treatment_2', 'Treatment_3']

    # Create count matrix with some realistic patterns
    count_matrix = pd.DataFrame(
        index=genes,
        columns=samples,
        data=np.random.negative_binomial(50, 0.1, (1000, 5))
    )

    # Add some zero counts
    zero_mask = np.random.random((1000, 5)) < 0.1
    count_matrix[zero_mask] = 0

    return count_matrix


def demonstrate_mageck_analysis():
    """Demonstrate MAGeCK analysis workflow."""
    logger.info("=== MAGeCK Analysis Demonstration ===")

    # Create example data
    screen_data = create_example_screen_data()
    logger.info(f"Created example screen data with {len(screen_data)} genes")

    # Initialize MAGeCK analyzer (without actual MAGeCK installation)
    try:
        mageck = MAGeCKAnalyzer()
        logger.info("MAGeCK analyzer initialized successfully")
    except RuntimeError as e:
        logger.warning(f"MAGeCK not installed: {e}")
        logger.info("Continuing with mock analysis...")

    # Demonstrate senescence pathway analysis using our mock data
    # Convert screen_data to MAGeCK-like format
    mageck_results_mock = {
        'gene_summary': screen_data.rename(columns={
            'gene': 'id',
            'log2fc': 'pos|lfc',
            'pvalue': 'pos|p-value',
            'fdr': 'pos|fdr',
            'rank': 'pos|rank'
        }),
        'quality_metrics': {
            'total_genes': len(screen_data),
            'significant_genes_fdr01': len(screen_data[screen_data['fdr'] < 0.1]),
            'mean_gene_lfc': screen_data['log2fc'].mean()
        }
    }

    logger.info("Analyzing senescence pathways...")

    # Identify significant senescence hits
    significant_hits = screen_data[screen_data['fdr'] < 0.1]
    senescence_hits = significant_hits[
        significant_hits['gene'].str.contains('HMGB1|AGER|JAK|STAT|NFKB|CDKN|TP53')
    ]

    logger.info(f"Found {len(significant_hits)} significant hits overall")
    logger.info(f"Found {len(senescence_hits)} significant senescence-related hits")

    if len(senescence_hits) > 0:
        logger.info("Top senescence hits:")
        for _, row in senescence_hits.head().iterrows():
            logger.info(f"  {row['gene']}: LFC={row['log2fc']:.2f}, FDR={row['fdr']:.3f}")


def demonstrate_qc_analysis():
    """Demonstrate screen quality control analysis."""
    logger.info("\n=== Screen Quality Control Demonstration ===")

    # Create example count matrix
    count_matrix = create_example_count_matrix()
    logger.info(f"Created count matrix: {count_matrix.shape[0]} genes x {count_matrix.shape[1]} samples")

    # Initialize QC analyzer
    qc_analyzer = ScreenQualityControl()

    # Run QC analysis
    qc_results = qc_analyzer.analyze_count_matrix(count_matrix)

    logger.info("QC Analysis Results:")
    logger.info(f"  Sample read count CV: {qc_results['read_distribution']['sample_read_counts']['cv']:.3f}")
    logger.info(f"  Mean detection rate: {qc_results['library_representation']['mean_detection_rate']:.3f}")
    logger.info(f"  Zero count percentage: {qc_results['zero_count_analysis']['zero_count_percentage']:.1%}")

    logger.info("Recommendations:")
    for i, rec in enumerate(qc_results['recommendations'], 1):
        logger.info(f"  {i}. {rec}")


def demonstrate_pathway_enrichment():
    """Demonstrate pathway enrichment analysis."""
    logger.info("\n=== Pathway Enrichment Demonstration ===")

    # Create example screen data
    screen_data = create_example_screen_data()

    # Initialize pathway analyzer
    pathway_analyzer = PathwayEnrichment()

    # Run enrichment analysis
    enrichment_results = pathway_analyzer.run_enrichment_analysis(
        screen_results=screen_data,
        fdr_threshold=0.1
    )

    # Create summary
    summary_df = pathway_analyzer.create_enrichment_summary(enrichment_results)

    logger.info("Pathway enrichment analysis completed")
    logger.info(f"Found {len(summary_df)} enriched pathways")

    if len(summary_df) > 0:
        logger.info("Top enriched pathways:")
        for _, row in summary_df.head().iterrows():
            logger.info(f"  {row['database']}.{row['pathway']}: "
                       f"hits={row['num_hits']}, p={row['p_value']:.3e}")


def demonstrate_senescence_analysis():
    """Demonstrate senescence-specific analysis."""
    logger.info("\n=== Senescence-Specific Analysis Demonstration ===")

    # Create example screen data
    screen_data = create_example_screen_data()

    # Initialize senescence analyzer
    senescence_analyzer = SenescenceScreenAnalyzer()

    # Run senescence analysis
    senescence_results = senescence_analyzer.analyze_senescence_screen(
        screen_results=screen_data,
        fdr_threshold=0.1,
        lfc_threshold=0.5
    )

    logger.info("Senescence analysis completed")

    # ReHMGB1 pathway analysis
    rehmgb1_results = senescence_results['rehmgb1_pathway_hits']
    logger.info("ReHMGB1 pathway analysis:")
    logger.info(f"  Total ReHMGB1 genes: {rehmgb1_results['total_rehmgb1_genes']}")
    logger.info(f"  Significant hits: {rehmgb1_results['significant_hits']}")

    if rehmgb1_results['significant_hits'] > 0:
        logger.info("  Key regulators:")
        for hit in rehmgb1_results['key_regulators']:
            logger.info(f"    {hit['gene']}: LFC={hit['log2fc']:.2f}, FDR={hit['fdr']:.3f}")

    # Therapeutic targets
    therapeutic_targets = senescence_results['therapeutic_targets']
    logger.info("\nTherapeutic target identification:")
    logger.info(f"  Senescence suppressors: {len(therapeutic_targets['senescence_suppressors'])}")

    if therapeutic_targets['senescence_suppressors']:
        logger.info("  Top targets:")
        for target in therapeutic_targets['senescence_suppressors'][:3]:
            logger.info(f"    {target['gene']}: LFC={target['log2fc']:.2f}, FDR={target['fdr']:.3f}")


def main():
    """Main demonstration function."""
    logger.info("CRISPR Toolkit - MAGeCK Integration Demonstration")
    logger.info("=" * 60)

    try:
        # Run all demonstrations
        demonstrate_mageck_analysis()
        demonstrate_qc_analysis()
        demonstrate_pathway_enrichment()
        demonstrate_senescence_analysis()

        logger.info("\n" + "=" * 60)
        logger.info("Demonstration completed successfully!")
        logger.info("\nKey features demonstrated:")
        logger.info("✓ MAGeCK integration for screen analysis")
        logger.info("✓ Comprehensive quality control metrics")
        logger.info("✓ Pathway enrichment analysis")
        logger.info("✓ Senescence-specific pathway analysis")
        logger.info("✓ ReHMGB1 pathway investigation")
        logger.info("✓ Therapeutic target identification")

    except Exception as e:
        logger.error(f"Demonstration failed: {e}")
        raise


if __name__ == "__main__":
    main()
    main()
