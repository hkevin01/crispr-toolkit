"""
Screen Quality Control Module

Provides comprehensive quality control metrics and validation for CRISPR screens.
Includes metrics for library representation, read distribution, and data quality.
"""

import logging
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)


class ScreenQualityControl:
    """
    Quality control analysis for CRISPR screens with focus on
    aging and senescence research applications.
    """

    def __init__(self):
        """Initialize screen QC analyzer."""
        self.qc_metrics = {}
        self.plots = {}

    def analyze_count_matrix(self,
                           count_matrix: pd.DataFrame,
                           library_info: Optional[pd.DataFrame] = None
                           ) -> Dict:
        """
        Perform comprehensive QC analysis on count matrix.

        Args:
            count_matrix: Count matrix with genes as rows, samples as columns
            library_info: Optional library information with gene annotations

        Returns:
            Dictionary containing QC metrics and recommendations
        """
        logger.info("Performing count matrix quality control analysis")

        qc_results = {
            'read_distribution': self._analyze_read_distribution(count_matrix),
            'library_representation': self._analyze_library_representation(
                count_matrix
            ),
            'sample_correlation': self._calculate_sample_correlation(
                count_matrix
            ),
            'zero_count_analysis': self._analyze_zero_counts(count_matrix),
            'replicate_consistency': self._analyze_replicate_consistency(
                count_matrix
            ),
            'recommendations': []
        }

        # Generate recommendations based on QC metrics
        qc_results['recommendations'] = self._generate_recommendations(
            qc_results
        )

        self.qc_metrics = qc_results
        return qc_results

    def _analyze_read_distribution(self,
                                 count_matrix: pd.DataFrame) -> Dict:
        """Analyze read count distribution across samples."""
        sample_totals = count_matrix.sum(axis=0)
        gene_totals = count_matrix.sum(axis=1)

        return {
            'sample_read_counts': {
                'mean': float(sample_totals.mean()),
                'std': float(sample_totals.std()),
                'min': float(sample_totals.min()),
                'max': float(sample_totals.max()),
                'cv': float(sample_totals.std() / sample_totals.mean())
            },
            'gene_read_counts': {
                'mean': float(gene_totals.mean()),
                'std': float(gene_totals.std()),
                'min': float(gene_totals.min()),
                'max': float(gene_totals.max()),
                'cv': float(gene_totals.std() / gene_totals.mean())
            }
        }

    def _analyze_library_representation(self,
                                      count_matrix: pd.DataFrame) -> Dict:
        """Analyze library representation and sgRNA coverage."""
        # Calculate percentage of genes with sufficient coverage
        genes_with_reads = (count_matrix > 0).sum(axis=1)
        total_samples = count_matrix.shape[1]

        coverage_stats = {
            'total_genes': count_matrix.shape[0],
            'total_samples': total_samples,
            'genes_detected_all_samples': int(
                (genes_with_reads == total_samples).sum()
            ),
            'genes_detected_50pct_samples': int(
                (genes_with_reads >= total_samples * 0.5).sum()
            ),
            'mean_detection_rate': float(genes_with_reads.mean() / total_samples)
        }

        return coverage_stats

    def _calculate_sample_correlation(self,
                                    count_matrix: pd.DataFrame) -> Dict:
        """Calculate correlation between samples."""
        # Log transform for correlation calculation
        log_counts = np.log2(count_matrix + 1)
        correlation_matrix = log_counts.corr()

        return {
            'correlation_matrix': correlation_matrix.to_dict(),
            'mean_correlation': float(
                correlation_matrix.values[
                    np.triu_indices_from(correlation_matrix.values, k=1)
                ].mean()
            ),
            'min_correlation': float(
                correlation_matrix.values[
                    np.triu_indices_from(correlation_matrix.values, k=1)
                ].min()
            )
        }

    def _analyze_zero_counts(self, count_matrix: pd.DataFrame) -> Dict:
        """Analyze zero count patterns."""
        zero_counts = (count_matrix == 0).sum()
        total_entries = count_matrix.shape[0] * count_matrix.shape[1]

        return {
            'zero_count_percentage': float(zero_counts.sum() / total_entries),
            'samples_with_zeros': zero_counts.to_dict(),
            'genes_with_all_zeros': int((count_matrix == 0).all(axis=1).sum())
        }

    def _analyze_replicate_consistency(self,
                                     count_matrix: pd.DataFrame) -> Dict:
        """Analyze consistency between replicates."""
        # This is a simplified version - would need sample grouping info
        # for proper replicate analysis
        log_counts = np.log2(count_matrix + 1)
        cv_per_gene = (log_counts.std(axis=1) / log_counts.mean(axis=1))

        return {
            'mean_cv_across_genes': float(cv_per_gene.mean()),
            'high_variability_genes': int((cv_per_gene > 0.5).sum()),
            'low_variability_genes': int((cv_per_gene < 0.1).sum())
        }

    def _generate_recommendations(self, qc_results: Dict) -> List[str]:
        """Generate QC recommendations based on analysis results."""
        recommendations = []

        # Check read count distribution
        read_cv = qc_results['read_distribution']['sample_read_counts']['cv']
        if read_cv > 0.3:
            recommendations.append(
                f"High coefficient of variation in sample read counts "
                f"(CV={read_cv:.3f}). Consider normalization or "
                f"removing outlier samples."
            )

        # Check library representation
        detection_rate = qc_results['library_representation'][
            'mean_detection_rate'
        ]
        if detection_rate < 0.8:
            recommendations.append(
                f"Low library representation (detection rate={detection_rate:.3f}). "
                f"Consider increasing sequencing depth."
            )

        # Check sample correlation
        min_corr = qc_results['sample_correlation']['min_correlation']
        if min_corr < 0.7:
            recommendations.append(
                f"Low sample correlation detected (min={min_corr:.3f}). "
                f"Check for batch effects or sample mix-ups."
            )

        # Check zero counts
        zero_pct = qc_results['zero_count_analysis']['zero_count_percentage']
        if zero_pct > 0.5:
            recommendations.append(
                f"High percentage of zero counts ({zero_pct:.1%}). "
                f"Consider filtering low-count genes."
            )

        if not recommendations:
            recommendations.append("Screen quality metrics look good!")

        return recommendations

    def create_qc_plots(self,
                       count_matrix: pd.DataFrame,
                       output_dir: str = "./qc_plots") -> Dict[str, str]:
        """
        Create quality control plots.

        Args:
            count_matrix: Count matrix for analysis
            output_dir: Directory to save plots

        Returns:
            Dictionary mapping plot names to file paths
        """
        import os
        os.makedirs(output_dir, exist_ok=True)

        plot_files = {}

        # Read distribution plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # Sample read totals
        sample_totals = count_matrix.sum(axis=0)
        axes[0, 0].hist(sample_totals, bins=20, alpha=0.7)
        axes[0, 0].set_title('Sample Read Count Distribution')
        axes[0, 0].set_xlabel('Total Reads')
        axes[0, 0].set_ylabel('Frequency')

        # Gene read totals
        gene_totals = count_matrix.sum(axis=1)
        axes[0, 1].hist(np.log10(gene_totals + 1), bins=30, alpha=0.7)
        axes[0, 1].set_title('Gene Read Count Distribution (log10)')
        axes[0, 1].set_xlabel('Log10(Total Reads + 1)')
        axes[0, 1].set_ylabel('Frequency')

        # Sample correlation heatmap
        log_counts = np.log2(count_matrix + 1)
        correlation_matrix = log_counts.corr()
        sns.heatmap(correlation_matrix, annot=True, cmap='viridis',
                   ax=axes[1, 0])
        axes[1, 0].set_title('Sample Correlation Matrix')

        # Zero count analysis
        zero_counts_per_sample = (count_matrix == 0).sum(axis=0)
        axes[1, 1].bar(range(len(zero_counts_per_sample)),
                      zero_counts_per_sample)
        axes[1, 1].set_title('Zero Counts per Sample')
        axes[1, 1].set_xlabel('Sample Index')
        axes[1, 1].set_ylabel('Number of Zero Counts')

        plt.tight_layout()
        qc_plot_file = os.path.join(output_dir, 'screen_qc_overview.png')
        plt.savefig(qc_plot_file, dpi=300, bbox_inches='tight')
        plt.close()

        plot_files['qc_overview'] = qc_plot_file

        return plot_files

    def export_qc_report(self,
                        output_file: str = "screen_qc_report.txt") -> None:
        """Export QC report to text file."""
        if not self.qc_metrics:
            logger.warning("No QC metrics available. Run analyze_count_matrix first.")
            return

        with open(output_file, 'w') as f:
            f.write("CRISPR Screen Quality Control Report\n")
            f.write("=" * 50 + "\n\n")

            # Write read distribution metrics
            f.write("Read Distribution Analysis:\n")
            f.write("-" * 30 + "\n")
            read_dist = self.qc_metrics['read_distribution']
            f.write(f"Sample read counts - Mean: {read_dist['sample_read_counts']['mean']:.0f}, "
                   f"CV: {read_dist['sample_read_counts']['cv']:.3f}\n")
            f.write(f"Gene read counts - Mean: {read_dist['gene_read_counts']['mean']:.0f}, "
                   f"CV: {read_dist['gene_read_counts']['cv']:.3f}\n\n")

            # Write library representation
            f.write("Library Representation:\n")
            f.write("-" * 25 + "\n")
            lib_rep = self.qc_metrics['library_representation']
            f.write(f"Total genes: {lib_rep['total_genes']}\n")
            f.write(f"Detection rate: {lib_rep['mean_detection_rate']:.3f}\n")
            f.write(f"Genes detected in all samples: {lib_rep['genes_detected_all_samples']}\n\n")

            # Write recommendations
            f.write("Recommendations:\n")
            f.write("-" * 20 + "\n")
            for rec in self.qc_metrics['recommendations']:
                f.write(f"- {rec}\n")

        logger.info(f"QC report exported to {output_file}")
        logger.info(f"QC report exported to {output_file}")
