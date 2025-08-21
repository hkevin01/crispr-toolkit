"""
Senescence-Specific Screen Analysis Module

Specialized analysis tools for CRISPR screens targeting senescence pathways,
with focus on ReHMGB1, RAGE-mediated signaling, and aging-related processes.
"""

import logging
from typing import Dict, List

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class SenescenceScreenAnalyzer:
    """
    Specialized analyzer for senescence-related CRISPR screens.

    Provides targeted analysis for aging pathways including:
    - ReHMGB1/RAGE signaling
    - JAK/STAT pathway
    - NF-κB signaling
    - Cell cycle arrest markers
    - SASP factors
    """

    def __init__(self):
        """Initialize senescence screen analyzer."""
        self.senescence_pathways = self._load_senescence_pathways()
        self.aging_markers = self._load_aging_markers()

    def _load_senescence_pathways(self) -> Dict[str, List[str]]:
        """Load curated senescence pathway gene sets."""
        return {
            'rehmgb1_rage_signaling': [
                'HMGB1', 'AGER', 'MYD88', 'IRAK1', 'IRAK4',
                'TRAF6', 'TAK1', 'MAPK14', 'MAPK8', 'MAPK9'
            ],
            'jak_stat_pathway': [
                'JAK1', 'JAK2', 'JAK3', 'TYK2', 'STAT1', 'STAT2',
                'STAT3', 'STAT4', 'STAT5A', 'STAT5B', 'STAT6'
            ],
            'nfkb_signaling': [
                'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL',
                'IKBKA', 'IKBKB', 'IKBKG', 'NFKBIA', 'NFKBIB'
            ],
            'cell_cycle_arrest': [
                'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C',
                'CDKN2D', 'RB1', 'RBL1', 'RBL2', 'TP53', 'TP73'
            ],
            'sasp_factors': [
                'IL6', 'IL1A', 'IL1B', 'TNF', 'CXCL8', 'CCL2',
                'CCL20', 'CXCL1', 'CXCL2', 'MMP1', 'MMP3', 'MMP9'
            ],
            'dna_damage_response': [
                'ATM', 'ATR', 'CHEK1', 'CHEK2', 'H2AFX', 'MDC1',
                'TOPBP1', 'CLASPIN', '53BP1', 'BRCA1'
            ],
            'oxidative_stress': [
                'SOD1', 'SOD2', 'CAT', 'GPX1', 'GPX4', 'PRDX1',
                'PRDX3', 'PRDX4', 'NRF2', 'KEAP1'
            ]
        }

    def get_pathway_genes(self, pathway_name: str) -> List[str]:
        """Get genes for a specific pathway.

        Args:
            pathway_name: Name of the pathway (e.g., 'rehmgb1_rage_signaling')

        Returns:
            List of gene symbols for the pathway
        """
        return self.senescence_pathways.get(pathway_name, [])

    def _load_aging_markers(self) -> Dict[str, List[str]]:
        """Load aging-related marker gene sets."""
        return {
            'longevity_genes': [
                'SIRT1', 'SIRT3', 'SIRT6', 'FOXO1', 'FOXO3',
                'KLOTHO', 'IGF1R', 'MTOR', 'AMPK', 'PGC1A'
            ],
            'autophagy_markers': [
                'ATG5', 'ATG7', 'ATG12', 'BECN1', 'MAP1LC3A',
                'MAP1LC3B', 'SQSTM1', 'ULK1', 'PIK3C3'
            ],
            'mitochondrial_function': [
                'TFAM', 'NRF1', 'PPARGC1A', 'SIRT3', 'OPA1',
                'MFN1', 'MFN2', 'DRP1', 'PINK1', 'PARK2'
            ]
        }

    def analyze_senescence_screen(self,
                                screen_results: pd.DataFrame,
                                fdr_threshold: float = 0.1,
                                lfc_threshold: float = 0.5) -> Dict:
        """
        Comprehensive analysis of senescence screen results.

        Args:
            screen_results: Gene-level screen results with columns:
                          ['gene', 'log2fc', 'pvalue', 'fdr']
            fdr_threshold: FDR threshold for significance
            lfc_threshold: Log2 fold change threshold

        Returns:
            Dictionary containing senescence-specific analysis results
        """
        logger.info("Analyzing senescence screen results")

        analysis_results = {
            'pathway_enrichment': self._analyze_pathway_enrichment(
                screen_results, fdr_threshold, lfc_threshold
            ),
            'rehmgb1_pathway_hits': self._analyze_rehmgb1_pathway(
                screen_results, fdr_threshold
            ),
            'senescence_modulators': self._identify_senescence_modulators(
                screen_results, fdr_threshold, lfc_threshold
            ),
            'aging_pathway_analysis': self._analyze_aging_pathways(
                screen_results, fdr_threshold
            ),
            'therapeutic_targets': self._identify_therapeutic_targets(
                screen_results, fdr_threshold, lfc_threshold
            )
        }

        return analysis_results

    def _analyze_pathway_enrichment(self,
                                  screen_results: pd.DataFrame,
                                  fdr_threshold: float,
                                  lfc_threshold: float) -> Dict:
        """Analyze enrichment of senescence pathways."""
        enrichment_results = {}

        significant_hits = screen_results[
            (screen_results['fdr'] < fdr_threshold) &
            (abs(screen_results['log2fc']) > lfc_threshold)
        ]

        for pathway_name, pathway_genes in self.senescence_pathways.items():
            pathway_hits = significant_hits[
                significant_hits['gene'].isin(pathway_genes)
            ]

            # Calculate enrichment statistics
            total_pathway_genes = len(pathway_genes)
            pathway_genes_in_screen = len(
                screen_results[screen_results['gene'].isin(pathway_genes)]
            )
            significant_pathway_hits = len(pathway_hits)

            if pathway_genes_in_screen > 0:
                hit_rate = significant_pathway_hits / pathway_genes_in_screen

                enrichment_results[pathway_name] = {
                    'total_genes': total_pathway_genes,
                    'genes_in_screen': pathway_genes_in_screen,
                    'significant_hits': significant_pathway_hits,
                    'hit_rate': hit_rate,
                    'hit_genes': pathway_hits['gene'].tolist(),
                    'mean_log2fc': pathway_hits['log2fc'].mean() if len(pathway_hits) > 0 else 0,
                    'pathway_significance': self._calculate_pathway_significance(
                        screen_results, pathway_genes, fdr_threshold, lfc_threshold
                    )
                }

        return enrichment_results

    def _analyze_rehmgb1_pathway(self,
                               screen_results: pd.DataFrame,
                               fdr_threshold: float) -> Dict:
        """Focused analysis of ReHMGB1 pathway components."""
        rehmgb1_genes = [
            'HMGB1', 'AGER', 'TLR2', 'TLR4', 'MYD88', 'IRAK1',
            'TRAF6', 'MAPK14', 'NFKB1', 'RELA', 'JAK1', 'JAK2',
            'STAT1', 'STAT3', 'IL6', 'TNF', 'CXCL8'
        ]

        rehmgb1_hits = screen_results[
            (screen_results['gene'].isin(rehmgb1_genes)) &
            (screen_results['fdr'] < fdr_threshold)
        ].copy()

        if len(rehmgb1_hits) > 0:
            rehmgb1_hits = rehmgb1_hits.sort_values('fdr')

        return {
            'total_rehmgb1_genes': len(rehmgb1_genes),
            'significant_hits': len(rehmgb1_hits),
            'hit_details': rehmgb1_hits.to_dict('records'),
            'key_regulators': rehmgb1_hits.head(5).to_dict('records'),
            'pathway_disruption_score': self._calculate_pathway_disruption(
                rehmgb1_hits
            )
        }

    def _identify_senescence_modulators(self,
                                      screen_results: pd.DataFrame,
                                      fdr_threshold: float,
                                      lfc_threshold: float) -> Dict:
        """Identify genes that modulate senescence."""
        # Pro-senescence hits (positive LFC)
        pro_senescence = screen_results[
            (screen_results['fdr'] < fdr_threshold) &
            (screen_results['log2fc'] > lfc_threshold)
        ].copy()

        # Anti-senescence hits (negative LFC)
        anti_senescence = screen_results[
            (screen_results['fdr'] < fdr_threshold) &
            (screen_results['log2fc'] < -lfc_threshold)
        ].copy()

        # Categorize by pathway membership
        all_senescence_genes = set()
        for pathway_genes in self.senescence_pathways.values():
            all_senescence_genes.update(pathway_genes)

        return {
            'pro_senescence_hits': {
                'total': len(pro_senescence),
                'known_senescence_genes': len(
                    pro_senescence[
                        pro_senescence['gene'].isin(all_senescence_genes)
                    ]
                ),
                'novel_hits': pro_senescence[
                    ~pro_senescence['gene'].isin(all_senescence_genes)
                ].head(10).to_dict('records'),
                'top_hits': pro_senescence.nsmallest(10, 'fdr').to_dict('records')
            },
            'anti_senescence_hits': {
                'total': len(anti_senescence),
                'known_senescence_genes': len(
                    anti_senescence[
                        anti_senescence['gene'].isin(all_senescence_genes)
                    ]
                ),
                'novel_hits': anti_senescence[
                    ~anti_senescence['gene'].isin(all_senescence_genes)
                ].head(10).to_dict('records'),
                'top_hits': anti_senescence.nsmallest(10, 'fdr').to_dict('records')
            }
        }

    def _analyze_aging_pathways(self,
                              screen_results: pd.DataFrame,
                              fdr_threshold: float) -> Dict:
        """Analyze aging-related pathway enrichment."""
        aging_results = {}

        for pathway_name, pathway_genes in self.aging_markers.items():
            pathway_hits = screen_results[
                (screen_results['gene'].isin(pathway_genes)) &
                (screen_results['fdr'] < fdr_threshold)
            ]

            aging_results[pathway_name] = {
                'total_genes': len(pathway_genes),
                'significant_hits': len(pathway_hits),
                'hit_genes': pathway_hits['gene'].tolist(),
                'mean_effect': pathway_hits['log2fc'].mean() if len(pathway_hits) > 0 else 0
            }

        return aging_results

    def _identify_therapeutic_targets(self,
                                    screen_results: pd.DataFrame,
                                    fdr_threshold: float,
                                    lfc_threshold: float) -> Dict:
        """Identify potential therapeutic targets for senescence intervention."""
        # Targets that reduce senescence when knocked out (negative LFC)
        senescence_suppressors = screen_results[
            (screen_results['fdr'] < fdr_threshold) &
            (screen_results['log2fc'] < -lfc_threshold)
        ].copy()

        # Prioritize based on druggability and pathway importance
        druggable_targets = self._filter_druggable_targets(senescence_suppressors)

        return {
            'senescence_suppressors': senescence_suppressors.nsmallest(20, 'fdr').to_dict('records'),
            'druggable_targets': druggable_targets,
            'pathway_targets': self._categorize_targets_by_pathway(senescence_suppressors),
            'novel_targets': self._identify_novel_targets(senescence_suppressors)
        }

    def _calculate_pathway_significance(self,
                                      screen_results: pd.DataFrame,
                                      pathway_genes: List[str],
                                      fdr_threshold: float,
                                      lfc_threshold: float) -> float:
        """Calculate pathway-level significance using Fisher's exact test."""
        try:
            from scipy.stats import fisher_exact
        except ImportError:
            logger.warning("scipy not available, returning p-value of 1.0")
            return 1.0

        # Create contingency table
        total_genes = len(screen_results)
        significant_mask = (
            (screen_results['fdr'] < fdr_threshold) &
            (abs(screen_results['log2fc']) > lfc_threshold)
        )
        total_significant = significant_mask.sum()

        pathway_mask = screen_results['gene'].isin(pathway_genes)
        pathway_genes_in_screen = pathway_mask.sum()
        significant_pathway_genes = (significant_mask & pathway_mask).sum()

        # Fisher's exact test
        contingency_table = [
            [significant_pathway_genes,
             pathway_genes_in_screen - significant_pathway_genes],
            [total_significant - significant_pathway_genes,
             total_genes - pathway_genes_in_screen -
             (total_significant - significant_pathway_genes)]
        ]

        _, p_value = fisher_exact(contingency_table, alternative='greater')
        return float(p_value)

    def _calculate_pathway_disruption(self, pathway_hits: pd.DataFrame) -> float:
        """Calculate overall pathway disruption score."""
        if len(pathway_hits) == 0:
            return 0.0

        # Weight by significance and effect size
        disruption_scores = abs(pathway_hits['log2fc']) * -np.log10(pathway_hits['fdr'])
        return float(disruption_scores.mean())

    def _filter_druggable_targets(self, targets: pd.DataFrame) -> List[Dict]:
        """Filter for potentially druggable targets."""
        # Simplified druggability assessment based on known target classes
        druggable_classes = {
            'kinases', 'phosphatases', 'proteases', 'receptors',
            'channels', 'transporters', 'enzymes'
        }

        # This would typically use external druggability databases
        # For now, return top targets based on statistical significance
        return targets.nsmallest(10, 'fdr').to_dict('records')

    def _categorize_targets_by_pathway(self, targets: pd.DataFrame) -> Dict:
        """Categorize targets by senescence pathway membership."""
        pathway_targets = {}

        for pathway_name, pathway_genes in self.senescence_pathways.items():
            pathway_targets[pathway_name] = targets[
                targets['gene'].isin(pathway_genes)
            ].to_dict('records')

        return pathway_targets

    def _identify_novel_targets(self, targets: pd.DataFrame) -> List[Dict]:
        """Identify novel senescence targets not in known pathways."""
        all_known_genes = set()
        for pathway_genes in self.senescence_pathways.values():
            all_known_genes.update(pathway_genes)
        for pathway_genes in self.aging_markers.values():
            all_known_genes.update(pathway_genes)

        novel_targets = targets[
            ~targets['gene'].isin(all_known_genes)
        ]

        return novel_targets.nsmallest(10, 'fdr').to_dict('records')
        return novel_targets.nsmallest(10, 'fdr').to_dict('records')
