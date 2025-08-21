"""
MAGeCK Wrapper for CRISPR Screen Analysis

This module provides a Python interface to MAGeCK (Model-based Analysis of
Genome-wide CRISPR/Cas9 Knockout) for robust analysis of CRISPR screens.

Features:
- Count matrix processing
- Essential gene identification
- Drug sensitivity analysis
- Pathway enrichment integration
- Quality control metrics
"""

import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class MAGeCKResults:
    """Container for MAGeCK analysis results."""
    gene_summary: pd.DataFrame
    sgrna_summary: pd.DataFrame
    normalized_counts: pd.DataFrame
    log_fold_changes: pd.DataFrame
    quality_metrics: Dict[str, float]
    pathway_enrichment: Optional[Dict] = None


class MAGeCKAnalyzer:
    """
    Wrapper class for MAGeCK analysis with enhanced functionality for
    aging and senescence research.
    """

    def __init__(self,
                 mageck_path: str = "mageck",
                 temp_dir: Optional[str] = None):
        """
        Initialize MAGeCK analyzer.

        Args:
            mageck_path: Path to MAGeCK executable
            temp_dir: Temporary directory for intermediate files
        """
        self.mageck_path = mageck_path
        self.temp_dir = temp_dir or tempfile.gettempdir()
        self._validate_mageck_installation()

    def _validate_mageck_installation(self) -> None:
        """Validate that MAGeCK is properly installed."""
        try:
            result = subprocess.run(
                [self.mageck_path, "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("MAGeCK version: %s", result.stdout.strip())
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            raise RuntimeError(
                f"MAGeCK not found or not executable: {e}. "
                "Please install MAGeCK or check the path."
            ) from e

    def count_reads(self,
                    fastq_files: List[str],
                    library_file: str,
                    output_prefix: str,
                    sample_labels: Optional[List[str]] = None) -> str:
        """
        Count reads from FASTQ files using MAGeCK count.

        Args:
            fastq_files: List of FASTQ file paths
            library_file: sgRNA library file
            output_prefix: Output file prefix
            sample_labels: Optional sample labels

        Returns:
            Path to count file
        """
        cmd = [
            self.mageck_path, "count",
            "-l", library_file,
            "-n", output_prefix,
            "--fastq"
        ] + fastq_files

        if sample_labels:
            cmd.extend(["--sample-label", ",".join(sample_labels)])

        logger.info("Running MAGeCK count: %s", " ".join(cmd))

        try:
            subprocess.run(cmd, check=True, cwd=self.temp_dir)
            count_file = f"{output_prefix}.count.txt"
            return os.path.join(self.temp_dir, count_file)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MAGeCK count failed: {e}") from e

    def test_essential_genes(self,
                             count_file: str,
                             treatment_samples: List[str],
                             control_samples: List[str],
                             output_prefix: str,
                             **kwargs) -> MAGeCKResults:
        """
        Identify essential genes using MAGeCK test.

        Args:
            count_file: Count matrix file
            treatment_samples: Treatment sample names
            control_samples: Control sample names
            output_prefix: Output prefix
            **kwargs: Additional MAGeCK parameters

        Returns:
            MAGeCKResults object with analysis results
        """
        cmd = [
            self.mageck_path, "test",
            "-k", count_file,
            "-t", ",".join(treatment_samples),
            "-c", ",".join(control_samples),
            "-n", output_prefix
        ]

        # Add optional parameters
        for key, value in kwargs.items():
            if key.startswith("--"):
                cmd.extend([key, str(value)])
            else:
                cmd.extend([f"--{key.replace('_', '-')}", str(value)])

        logger.info("Running MAGeCK test: %s", " ".join(cmd))

        try:
            subprocess.run(cmd, check=True, cwd=self.temp_dir)
            return self._parse_mageck_results(output_prefix)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MAGeCK test failed: {e}") from e

    def drug_sensitivity_analysis(self,
                                  count_file: str,
                                  drug_conditions: Dict[str, List[str]],
                                  control_samples: List[str],
                                  output_prefix: str
                                  ) -> Dict[str, MAGeCKResults]:
        """
        Analyze drug sensitivity across multiple conditions.

        Args:
            count_file: Count matrix file
            drug_conditions: Dict mapping drug names to sample lists
            control_samples: Control sample names
            output_prefix: Output prefix

        Returns:
            Dictionary of drug names to MAGeCKResults
        """
        results = {}

        for drug_name, treatment_samples in drug_conditions.items():
            drug_output = f"{output_prefix}_{drug_name}"
            logger.info("Analyzing drug sensitivity for %s", drug_name)

            results[drug_name] = self.test_essential_genes(
                count_file=count_file,
                treatment_samples=treatment_samples,
                control_samples=control_samples,
                output_prefix=drug_output
            )

        return results

    def senescence_pathway_analysis(self,
                                   results: MAGeCKResults,
                                   senescence_genes: Optional[List[str]] = None
                                   ) -> Dict:
        """
        Analyze senescence-related pathways in screen results.

        Args:
            results: MAGeCK analysis results
            senescence_genes: List of senescence-related genes

        Returns:
            Senescence pathway analysis results
        """
        if senescence_genes is None:
            # Default senescence-related genes including ReHMGB1 pathway
            senescence_genes = [
                "HMGB1", "AGER", "JAK1", "JAK2", "STAT1", "STAT3",
                "NFKB1", "RELA", "CDKN1A", "CDKN2A", "TP53",
                "RB1", "IL6", "IL1B", "TNF", "CXCL8", "CCL2",
                "SASP", "ATM", "ATR", "CHEK1", "CHEK2"
            ]

        gene_summary = results.gene_summary

        # Filter for senescence genes
        senescence_hits = gene_summary[
            gene_summary['id'].isin(senescence_genes)
        ].copy()

        # Calculate enrichment statistics
        significant_mask = senescence_hits['pos|fdr'] < 0.1
        enrichment_stats = {
            'total_senescence_genes': len(senescence_genes),
            'detected_genes': len(senescence_hits),
            'significant_hits': len(senescence_hits[significant_mask]),
            'mean_log2fc': senescence_hits['pos|lfc'].mean(),
            'pathway_enrichment_p': self._calculate_pathway_enrichment(
                gene_summary, senescence_genes
            )
        }

        return {
            'senescence_hits': senescence_hits,
            'enrichment_stats': enrichment_stats,
            'top_senescence_regulators': senescence_hits.nsmallest(
                10, 'pos|fdr'
            )
        }

    def _parse_mageck_results(self, output_prefix: str) -> MAGeCKResults:
        """Parse MAGeCK output files into results object."""
        base_path = Path(self.temp_dir) / output_prefix

        # Load gene summary
        gene_summary = pd.read_csv(f"{base_path}.gene_summary.txt", sep="\t")

        # Load sgRNA summary
        sgrna_summary = pd.read_csv(f"{base_path}.sgrna_summary.txt", sep="\t")

        # Load normalized counts (if available)
        normalized_counts = pd.DataFrame()
        norm_file = f"{base_path}.normalized.txt"
        if os.path.exists(norm_file):
            normalized_counts = pd.read_csv(norm_file, sep="\t")

        # Calculate quality metrics
        quality_metrics = self._calculate_quality_metrics(
            gene_summary, sgrna_summary
        )

        # Calculate log fold changes
        log_fold_changes = self._calculate_log_fold_changes(gene_summary)

        return MAGeCKResults(
            gene_summary=gene_summary,
            sgrna_summary=sgrna_summary,
            normalized_counts=normalized_counts,
            log_fold_changes=log_fold_changes,
            quality_metrics=quality_metrics
        )

    def _calculate_quality_metrics(self,
                                  gene_summary: pd.DataFrame,
                                  sgrna_summary: pd.DataFrame
                                  ) -> Dict[str, float]:
        """Calculate screen quality control metrics."""
        fdr01_mask = gene_summary['pos|fdr'] < 0.1
        fdr05_mask = gene_summary['pos|fdr'] < 0.05

        return {
            'total_genes': len(gene_summary),
            'total_sgrnas': len(sgrna_summary),
            'mean_sgrnas_per_gene': len(sgrna_summary) / len(gene_summary),
            'significant_genes_fdr01': len(gene_summary[fdr01_mask]),
            'significant_genes_fdr05': len(gene_summary[fdr05_mask]),
            'mean_gene_lfc': gene_summary['pos|lfc'].mean(),
            'std_gene_lfc': gene_summary['pos|lfc'].std(),
        }

    def _calculate_log_fold_changes(self,
                                   gene_summary: pd.DataFrame
                                   ) -> pd.DataFrame:
        """Extract and format log fold changes."""
        return gene_summary[['id', 'pos|lfc', 'neg|lfc']].copy()

    def _calculate_pathway_enrichment(self,
                                    gene_summary: pd.DataFrame,
                                    pathway_genes: List[str]) -> float:
        """Calculate pathway enrichment p-value using Fisher's exact test."""
        try:
            from scipy.stats import fisher_exact
        except ImportError:
            logger.warning("scipy not available, returning p-value of 1.0")
            return 1.0

        # Create contingency table
        total_genes = len(gene_summary)
        significant_mask = gene_summary['pos|fdr'] < 0.1
        significant_genes = len(gene_summary[significant_mask])

        pathway_mask = gene_summary['id'].isin(pathway_genes)
        pathway_genes_in_screen = len(gene_summary[pathway_mask])

        significant_pathway_genes = len(
            gene_summary[pathway_mask & significant_mask]
        )

        # Fisher's exact test
        contingency_table = [
            [significant_pathway_genes,
             pathway_genes_in_screen - significant_pathway_genes],
            [significant_genes - significant_pathway_genes,
             total_genes - pathway_genes_in_screen -
             (significant_genes - significant_pathway_genes)]
        ]

        _, p_value = fisher_exact(contingency_table, alternative='greater')
        return p_value

    def create_visualization_data(self, results: MAGeCKResults) -> Dict:
        """
        Prepare data for visualization including volcano plots and heatmaps.

        Args:
            results: MAGeCK analysis results

        Returns:
            Dictionary containing visualization-ready data
        """
        gene_summary = results.gene_summary

        # Volcano plot data
        volcano_data = {
            'genes': gene_summary['id'].tolist(),
            'log2fc': gene_summary['pos|lfc'].tolist(),
            'neg_log10_pval': -np.log10(
                gene_summary['pos|p-value']
            ).tolist(),
            'fdr': gene_summary['pos|fdr'].tolist()
        }

        # Rank plot data
        rank_data = {
            'genes': gene_summary['id'].tolist(),
            'rank': gene_summary['pos|rank'].tolist(),
            'log2fc': gene_summary['pos|lfc'].tolist()
        }

        return {
            'volcano_plot': volcano_data,
            'rank_plot': rank_data,
            'quality_metrics': results.quality_metrics
        }
        }
