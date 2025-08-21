"""
Pathway Enrichment Analysis Module

Provides pathway enrichment analysis for CRISPR screen results
with support for multiple pathway databases and custom gene sets.
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class PathwayEnrichment:
    """
    Pathway enrichment analysis for CRISPR screen results.

    Supports multiple pathway databases and custom gene sets
    with focus on aging and senescence pathways.
    """

    def __init__(self):
        """Initialize pathway enrichment analyzer."""
        self.pathway_databases = self._load_pathway_databases()

    def _load_pathway_databases(self) -> Dict[str, Dict[str, List[str]]]:
        """Load pathway gene sets from various databases."""
        return {
            'kegg_aging': {
                'cellular_senescence': [
                    'TP53', 'CDKN1A', 'CDKN2A', 'RB1', 'ATM', 'ATR',
                    'CHEK1', 'CHEK2', 'H2AFX', 'MDC1'
                ],
                'foxo_signaling': [
                    'FOXO1', 'FOXO3', 'FOXO4', 'FOXO6', 'AKT1', 'AKT2',
                    'PIK3CA', 'PIK3CB', 'PTEN', 'SIRT1'
                ],
                'mtor_signaling': [
                    'MTOR', 'RICTOR', 'RAPTOR', 'MLST8', 'AKT1S1',
                    'RPS6KB1', 'EIF4EBP1', 'ULK1', 'ATG1'
                ]
            },
            'reactome_inflammation': {
                'nfkb_pathway': [
                    'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL',
                    'IKBKA', 'IKBKB', 'IKBKG', 'NFKBIA'
                ],
                'cytokine_signaling': [
                    'IL6', 'IL1A', 'IL1B', 'TNF', 'IFNG', 'IL10',
                    'JAK1', 'JAK2', 'STAT1', 'STAT3'
                ]
            },
            'custom_aging': {
                'longevity_genes': [
                    'SIRT1', 'SIRT3', 'SIRT6', 'KLOTHO', 'FOXO3',
                    'APOE', 'IGF1R', 'GH1', 'TERT', 'TERC'
                ],
                'mitochondrial_dysfunction': [
                    'TFAM', 'NRF1', 'PPARGC1A', 'SIRT3', 'OPA1',
                    'MFN1', 'MFN2', 'DRP1', 'PINK1', 'PARK2'
                ]
            }
        }

    def run_enrichment_analysis(self,
                              screen_results: pd.DataFrame,
                              fdr_threshold: float = 0.1,
                              databases: Optional[List[str]] = None
                              ) -> Dict:
        """
        Run pathway enrichment analysis on screen results.

        Args:
            screen_results: DataFrame with columns ['gene', 'log2fc', 'fdr']
            fdr_threshold: FDR threshold for significant genes
            databases: List of databases to use (default: all)

        Returns:
            Dictionary containing enrichment results
        """
        if databases is None:
            databases = list(self.pathway_databases.keys())

        logger.info("Running pathway enrichment analysis")

        significant_genes = screen_results[
            screen_results['fdr'] < fdr_threshold
        ]['gene'].tolist()

        enrichment_results = {}

        for db_name in databases:
            if db_name not in self.pathway_databases:
                logger.warning(f"Database {db_name} not found, skipping")
                continue

            db_results = {}
            pathways = self.pathway_databases[db_name]

            for pathway_name, pathway_genes in pathways.items():
                enrichment = self._calculate_enrichment(
                    significant_genes=significant_genes,
                    pathway_genes=pathway_genes,
                    total_genes=len(screen_results),
                    screen_results=screen_results
                )
                db_results[pathway_name] = enrichment

            enrichment_results[db_name] = db_results

        return enrichment_results

    def _calculate_enrichment(self,
                            significant_genes: List[str],
                            pathway_genes: List[str],
                            total_genes: int,
                            screen_results: pd.DataFrame) -> Dict:
        """Calculate enrichment statistics for a pathway."""
        # Intersection of significant genes and pathway genes
        pathway_hits = list(set(significant_genes) & set(pathway_genes))

        # Calculate basic statistics
        pathway_size = len(pathway_genes)
        num_significant = len(significant_genes)
        num_hits = len(pathway_hits)

        # Calculate enrichment ratio
        expected_hits = (pathway_size * num_significant) / total_genes
        enrichment_ratio = num_hits / expected_hits if expected_hits > 0 else 0

        # Calculate p-value using hypergeometric test
        p_value = self._hypergeometric_test(
            num_hits, pathway_size, num_significant, total_genes
        )

        # Get effect directions for hits
        hit_effects = screen_results[
            screen_results['gene'].isin(pathway_hits)
        ][['gene', 'log2fc', 'fdr']].to_dict('records')

        return {
            'pathway_size': pathway_size,
            'num_hits': num_hits,
            'hit_genes': pathway_hits,
            'hit_effects': hit_effects,
            'enrichment_ratio': enrichment_ratio,
            'p_value': p_value,
            'expected_hits': expected_hits,
            'fold_enrichment': enrichment_ratio,
            'mean_log2fc': np.mean([h['log2fc'] for h in hit_effects]) if hit_effects else 0
        }

    def _hypergeometric_test(self,
                           num_hits: int,
                           pathway_size: int,
                           num_significant: int,
                           total_genes: int) -> float:
        """Calculate hypergeometric p-value for pathway enrichment."""
        try:
            from scipy.stats import hypergeom

            # Hypergeometric test
            # Population size: total_genes
            # Number of success states in population: pathway_size
            # Number of draws: num_significant
            # Number of observed successes: num_hits

            p_value = hypergeom.sf(
                num_hits - 1, total_genes, pathway_size, num_significant
            )
            return float(p_value)

        except ImportError:
            logger.warning("scipy not available, using Fisher's exact test")
            return self._fisher_exact_test(
                num_hits, pathway_size, num_significant, total_genes
            )

    def _fisher_exact_test(self,
                         num_hits: int,
                         pathway_size: int,
                         num_significant: int,
                         total_genes: int) -> float:
        """Fallback Fisher's exact test for enrichment."""
        try:
            from scipy.stats import fisher_exact

            # Contingency table for Fisher's exact test
            contingency_table = [
                [num_hits, pathway_size - num_hits],
                [num_significant - num_hits,
                 total_genes - pathway_size - (num_significant - num_hits)]
            ]

            _, p_value = fisher_exact(contingency_table, alternative='greater')
            return float(p_value)

        except ImportError:
            logger.warning("scipy not available, returning p-value of 1.0")
            return 1.0

    def create_enrichment_summary(self, enrichment_results: Dict) -> pd.DataFrame:
        """Create a summary table of enrichment results."""
        summary_data = []

        for db_name, db_results in enrichment_results.items():
            for pathway_name, pathway_data in db_results.items():
                if pathway_data['num_hits'] > 0:  # Only include pathways with hits
                    summary_data.append({
                        'database': db_name,
                        'pathway': pathway_name,
                        'pathway_size': pathway_data['pathway_size'],
                        'num_hits': pathway_data['num_hits'],
                        'enrichment_ratio': pathway_data['enrichment_ratio'],
                        'p_value': pathway_data['p_value'],
                        'mean_log2fc': pathway_data['mean_log2fc'],
                        'hit_genes': ';'.join(pathway_data['hit_genes'])
                    })

        summary_df = pd.DataFrame(summary_data)

        if len(summary_df) > 0:
            # Add multiple testing correction
            summary_df['fdr'] = self._multiple_testing_correction(
                summary_df['p_value'].tolist()
            )

            # Sort by significance
            summary_df = summary_df.sort_values('p_value')

        return summary_df

    def _multiple_testing_correction(self, p_values: List[float]) -> List[float]:
        """Apply Benjamini-Hochberg FDR correction."""
        try:
            from scipy.stats import false_discovery_control
            return false_discovery_control(p_values).tolist()
        except ImportError:
            # Simple Bonferroni correction as fallback
            n = len(p_values)
            return [min(p * n, 1.0) for p in p_values]

    def add_custom_pathways(self,
                          database_name: str,
                          pathways: Dict[str, List[str]]) -> None:
        """Add custom pathway gene sets."""
        self.pathway_databases[database_name] = pathways
        logger.info(f"Added custom pathway database: {database_name}")

    def get_pathway_genes(self,
                        database: str,
                        pathway: str) -> Optional[List[str]]:
        """Get genes for a specific pathway."""
        if database in self.pathway_databases:
            return self.pathway_databases[database].get(pathway)
        return None

    def export_enrichment_results(self,
                                enrichment_results: Dict,
                                output_file: str) -> None:
        """Export enrichment results to file."""
        summary_df = self.create_enrichment_summary(enrichment_results)

        if len(summary_df) > 0:
            summary_df.to_csv(output_file, index=False)
            logger.info(f"Enrichment results exported to {output_file}")
        else:
            logger.warning("No enrichment results to export")
            logger.warning("No enrichment results to export")
