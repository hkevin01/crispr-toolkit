"""
Pathway Enrichment Analysis Module for CRISPR Toolkit Phase 3
=============================================================

Advanced pathway enrichment analysis for multi-omics aging intervention
research including Gene Ontology, KEGG, Reactome, and aging-specific
pathway databases.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Set

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)


@dataclass
class PathwayEnrichmentResult:
    """Results from pathway enrichment analysis."""
    pathway_id: str
    pathway_name: str
    p_value: float
    adjusted_p_value: float
    enrichment_score: float
    genes_in_pathway: List[str]
    genes_in_data: List[str]
    overlap_genes: List[str]
    fold_enrichment: float


@dataclass
class AgingPathwayScore:
    """Aging-specific pathway activity scores."""
    pathway_name: str
    activity_score: float
    intervention_response: float
    confidence_score: float
    key_genes: List[str]


class PathwayEnrichmentAnalyzer:
    """
    Comprehensive pathway enrichment analysis for aging intervention research.

    Integrates multiple pathway databases and performs enrichment analysis
    for transcriptomics, proteomics, and metabolomics data.
    """

    def __init__(self):
        """Initialize pathway enrichment analyzer."""
        self.pathway_databases = self._load_pathway_databases()
        self.aging_pathways = self._load_aging_specific_pathways()
        self.intervention_pathways = self._load_intervention_pathways()

        logger.info("Initialized PathwayEnrichmentAnalyzer")

    def _load_pathway_databases(self) -> Dict[str, Dict[str, List[str]]]:
        """Load pathway databases (simplified version)."""

        # Gene Ontology Biological Processes (aging-related)
        go_bp = {
            'GO:0016265': ['TP53', 'CDKN1A', 'CDKN2A', 'RB1', 'ATM'],  # Cell death
            'GO:0006915': ['BCL2', 'BAX', 'CASP3', 'CASP9', 'APAF1'],  # Apoptosis
            'GO:0016446': ['SIRT1', 'SIRT3', 'SIRT6', 'NAMPT', 'PARP1'],  # Somatic mutation
            'GO:0006950': ['HSP90', 'HSP70', 'HSF1', 'HSPA1A', 'DNAJB1'],  # Stress response
            'GO:0016887': ['ATG5', 'ATG7', 'BECN1', 'LC3B', 'SQSTM1'],  # Autophagy
            'GO:0007049': ['CDK1', 'CDK2', 'CCND1', 'CCNE1', 'CDKN1B'],  # Cell cycle
        }

        # KEGG Pathways (aging and longevity related)
        kegg = {
            'hsa04210': ['FOXO1', 'FOXO3', 'SIRT1', 'IGF1', 'IRS1'],  # Apoptosis
            'hsa04068': ['FOXO1', 'FOXO3', 'FOXO4', 'FOXO6', 'PI3K'],  # FoxO signaling
            'hsa04151': ['RHEB', 'TSC1', 'TSC2', 'MTOR', 'RPS6KB1'],  # PI3K-Akt
            'hsa04152': ['AKT1', 'AKT2', 'PIK3CA', 'PTEN', 'GSK3B'],  # AMPK signaling
            'hsa04150': ['PRKAA1', 'PRKAB1', 'PRKAG1', 'PPARA', 'ACADM'],  # mTOR signaling
            'hsa04211': ['TP53', 'MDM2', 'CDKN1A', 'BAX', 'PUMA'],  # Longevity regulating
        }

        # Reactome Pathways
        reactome = {
            'R-HSA-73857': ['TP53', 'ATM', 'CHEK2', 'BRCA1', 'BRCA2'],  # DNA Repair
            'R-HSA-1640170': ['TP53', 'MDM2', 'CDKN1A', 'CDKN2A'],  # Cell Cycle Checkpoints
            'R-HSA-2559583': ['SIRT1', 'SIRT3', 'SIRT6', 'NAD+', 'NAMPT'],  # Cellular Senescence
            'R-HSA-380108': ['CASP3', 'CASP8', 'CASP9', 'BCL2', 'BAX'],  # Apoptosis
            'R-HSA-9006934': ['PINK1', 'PARK2', 'TFEB', 'LAMP1'],  # Mitophagy
            'R-HSA-1655829': ['TOR', 'ULK1', 'ATG13', 'FIP200', 'ATG101'],  # Autophagy regulation
        }

        return {
            'GO_BP': go_bp,
            'KEGG': kegg,
            'Reactome': reactome
        }

    def _load_aging_specific_pathways(self) -> Dict[str, List[str]]:
        """Load aging-specific pathway signatures."""

        aging_pathways = {
            'cellular_senescence': [
                'TP53', 'CDKN1A', 'CDKN2A', 'CDKN2B', 'RB1', 'SASP_IL6',
                'SASP_IL1B', 'SASP_TNF', 'SASP_CXCL1', 'SASP_CCL2'
            ],
            'telomere_maintenance': [
                'TERT', 'TERC', 'TRF1', 'TRF2', 'POT1', 'TINF2',
                'TPP1', 'RAP1', 'RTEL1', 'CTC1'
            ],
            'DNA_damage_response': [
                'ATM', 'ATR', 'CHEK1', 'CHEK2', 'TP53', 'BRCA1',
                'BRCA2', 'RAD51', 'XRCC1', 'PARP1'
            ],
            'mitochondrial_dysfunction': [
                'PGC1A', 'NRF1', 'TFAM', 'COX4I1', 'CYCS',
                'PINK1', 'PARK2', 'MFN1', 'MFN2', 'DRP1'
            ],
            'nutrient_sensing': [
                'MTOR', 'RICTOR', 'RAPTOR', 'PRKAA1', 'PRKAB1',
                'SIRT1', 'SIRT3', 'FOXO1', 'FOXO3', 'IGF1'
            ],
            'proteostasis': [
                'HSP90', 'HSP70', 'HSPA1A', 'HSF1', 'PSMD1',
                'PSMB5', 'UBE2I', 'ATG5', 'ATG7', 'BECN1'
            ],
            'stem_cell_exhaustion': [
                'NANOG', 'OCT4', 'SOX2', 'KLF4', 'CCND1',
                'WNT3A', 'CTNNB1', 'TCF7L2', 'LEF1', 'MYC'
            ],
            'altered_communication': [
                'IL6', 'TNF', 'IL1B', 'CXCL1', 'CCL2',
                'NFKB1', 'RELA', 'STAT3', 'JAK2', 'SOCS3'
            ],
            'epigenetic_alterations': [
                'DNMT1', 'DNMT3A', 'DNMT3B', 'TET1', 'TET2',
                'EZH2', 'SUZ12', 'KDM6A', 'KDM6B', 'SIRT6'
            ]
        }

        return aging_pathways

    def _load_intervention_pathways(self) -> Dict[str, List[str]]:
        """Load intervention-specific pathway signatures."""

        intervention_pathways = {
            'caloric_restriction': [
                'SIRT1', 'FOXO1', 'FOXO3', 'PGC1A', 'ADIPOQ',
                'LEPR', 'PPARA', 'PPARG', 'UCP1', 'FGF21'
            ],
            'exercise_response': [
                'PGC1A', 'NRF1', 'TFAM', 'VEGFA', 'MYOD1',
                'MEF2C', 'PPARA', 'AMPK', 'mTOR', 'IGF1'
            ],
            'rapamycin_target': [
                'MTOR', 'RPS6KB1', 'EIF4EBP1', 'ULK1', 'ATG13',
                'RICTOR', 'RAPTOR', 'DEPTOR', 'MLST8', 'AKT1S1'
            ],
            'metformin_response': [
                'PRKAA1', 'PRKAB1', 'PRKAG1', 'ACC1', 'ACC2',
                'HMGCR', 'SREBF1', 'PPARA', 'FOXO1', 'G6PC'
            ],
            'resveratrol_target': [
                'SIRT1', 'SIRT3', 'SIRT6', 'PGC1A', 'NF1',
                'FOXO1', 'FOXO3', 'PPARA', 'ADIPOQ', 'NAMPT'
            ],
            'nad_precursor': [
                'NAMPT', 'NMNAT1', 'NMNAT2', 'NMNAT3', 'SIRT1',
                'SIRT3', 'SIRT6', 'PARP1', 'CD38', 'NNMT'
            ]
        }

        return intervention_pathways

    def perform_enrichment_analysis(self,
                                  gene_list: List[str],
                                  background_genes: Optional[List[str]] = None,
                                  databases: List[str] = None) -> List[PathwayEnrichmentResult]:
        """
        Perform pathway enrichment analysis.

        Args:
            gene_list: List of genes of interest
            background_genes: Background gene universe
            databases: Databases to use ('GO_BP', 'KEGG', 'Reactome')

        Returns:
            List of enrichment results
        """
        if databases is None:
            databases = ['GO_BP', 'KEGG', 'Reactome']

        logger.info(f"Performing enrichment analysis for {len(gene_list)} genes")

        # Convert to uppercase for consistency
        gene_set = set([gene.upper() for gene in gene_list])

        # Set background if not provided
        if background_genes is None:
            background_genes = self._get_default_background()

        background_set = set([gene.upper() for gene in background_genes])

        enrichment_results = []

        for db_name in databases:
            if db_name in self.pathway_databases:
                db_results = self._analyze_database(
                    gene_set, background_set,
                    self.pathway_databases[db_name], db_name
                )
                enrichment_results.extend(db_results)

        # Sort by p-value
        enrichment_results.sort(key=lambda x: x.p_value)

        # Apply multiple testing correction
        enrichment_results = self._apply_multiple_testing_correction(enrichment_results)

        return enrichment_results

    def _get_default_background(self) -> List[str]:
        """Get default background gene set."""

        all_genes = set()

        # Collect genes from all pathways
        for db in self.pathway_databases.values():
            for pathway_genes in db.values():
                all_genes.update([gene.upper() for gene in pathway_genes])

        # Add aging pathway genes
        for pathway_genes in self.aging_pathways.values():
            all_genes.update([gene.upper() for gene in pathway_genes])

        return list(all_genes)

    def _analyze_database(self,
                         gene_set: Set[str],
                         background_set: Set[str],
                         database: Dict[str, List[str]],
                         db_name: str) -> List[PathwayEnrichmentResult]:
        """Analyze enrichment for a specific database."""

        results = []

        for pathway_id, pathway_genes in database.items():
            pathway_gene_set = set([gene.upper() for gene in pathway_genes])

            # Calculate overlap
            overlap_genes = gene_set.intersection(pathway_gene_set)
            genes_in_background = pathway_gene_set.intersection(background_set)

            if len(overlap_genes) == 0:
                continue

            # Hypergeometric test
            N = len(background_set)  # Total genes in background
            K = len(genes_in_background)  # Pathway genes in background
            n = len(gene_set)  # Query genes
            k = len(overlap_genes)  # Overlap

            if K > 0 and k > 0:
                # Calculate p-value using hypergeometric distribution
                p_value = hypergeom.sf(k - 1, N, K, n)

                # Calculate fold enrichment
                expected = (n * K) / N
                fold_enrichment = k / max(expected, 0.001)

                # Calculate enrichment score
                enrichment_score = -np.log10(max(p_value, 1e-300))

                # Get pathway name
                pathway_name = self._get_pathway_name(pathway_id, db_name)

                result = PathwayEnrichmentResult(
                    pathway_id=f"{db_name}:{pathway_id}",
                    pathway_name=pathway_name,
                    p_value=p_value,
                    adjusted_p_value=p_value,  # Will be corrected later
                    enrichment_score=enrichment_score,
                    genes_in_pathway=list(genes_in_background),
                    genes_in_data=list(gene_set),
                    overlap_genes=list(overlap_genes),
                    fold_enrichment=fold_enrichment
                )

                results.append(result)

        return results

    def _get_pathway_name(self, pathway_id: str, db_name: str) -> str:
        """Get human-readable pathway name."""

        # Mapping of pathway IDs to names (simplified)
        pathway_names = {
            'GO:0016265': 'Death',
            'GO:0006915': 'Apoptotic process',
            'GO:0016446': 'Somatic mutation',
            'GO:0006950': 'Response to stress',
            'GO:0016887': 'ATPase activity',
            'GO:0007049': 'Cell cycle',
            'hsa04210': 'Apoptosis',
            'hsa04068': 'FoxO signaling pathway',
            'hsa04151': 'PI3K-Akt signaling pathway',
            'hsa04152': 'AMPK signaling pathway',
            'hsa04150': 'mTOR signaling pathway',
            'hsa04211': 'Longevity regulating pathway',
            'R-HSA-73857': 'DNA Repair',
            'R-HSA-1640170': 'Cell Cycle Checkpoints',
            'R-HSA-2559583': 'Cellular Senescence',
            'R-HSA-380108': 'Apoptosis',
            'R-HSA-9006934': 'Mitophagy',
            'R-HSA-1655829': 'Regulation of autophagy'
        }

        return pathway_names.get(pathway_id, pathway_id)

    def _apply_multiple_testing_correction(self,
                                         results: List[PathwayEnrichmentResult]) -> List[PathwayEnrichmentResult]:
        """Apply Benjamini-Hochberg correction."""

        if len(results) == 0:
            return results

        # Extract p-values
        p_values = [result.p_value for result in results]

        # Benjamini-Hochberg correction
        n = len(p_values)
        sorted_indices = np.argsort(p_values)

        corrected_p_values = [0] * n

        for i, idx in enumerate(sorted_indices):
            rank = i + 1
            corrected_p = min(1.0, p_values[idx] * n / rank)
            corrected_p_values[idx] = corrected_p

        # Update results with corrected p-values
        for i, result in enumerate(results):
            result.adjusted_p_value = corrected_p_values[i]

        return results

    def analyze_aging_pathways(self,
                             gene_expression: pd.DataFrame,
                             condition_comparison: Optional[Dict[str, List[str]]] = None) -> List[AgingPathwayScore]:
        """
        Analyze aging-specific pathway activity.

        Args:
            gene_expression: Gene expression matrix (genes x samples)
            condition_comparison: Dict mapping conditions to sample lists

        Returns:
            List of aging pathway scores
        """
        logger.info("Analyzing aging-specific pathway activity")

        aging_scores = []

        for pathway_name, pathway_genes in self.aging_pathways.items():
            # Find genes in the expression data
            available_genes = [gene for gene in pathway_genes
                             if gene in gene_expression.index]

            if len(available_genes) == 0:
                continue

            # Calculate pathway activity score
            pathway_expression = gene_expression.loc[available_genes]
            activity_score = pathway_expression.mean().mean()  # Mean across genes and samples

            # Calculate intervention response if conditions provided
            intervention_response = 0.0
            if condition_comparison:
                intervention_response = self._calculate_intervention_response(
                    pathway_expression, condition_comparison
                )

            # Calculate confidence score based on gene coverage
            confidence_score = len(available_genes) / len(pathway_genes)

            aging_score = AgingPathwayScore(
                pathway_name=pathway_name,
                activity_score=activity_score,
                intervention_response=intervention_response,
                confidence_score=confidence_score,
                key_genes=available_genes
            )

            aging_scores.append(aging_score)

        # Sort by activity score
        aging_scores.sort(key=lambda x: abs(x.intervention_response), reverse=True)

        return aging_scores

    def _calculate_intervention_response(self,
                                       pathway_expression: pd.DataFrame,
                                       condition_comparison: Dict[str, List[str]]) -> float:
        """Calculate pathway response to intervention."""

        if 'control' not in condition_comparison or 'treatment' not in condition_comparison:
            return 0.0

        control_samples = condition_comparison['control']
        treatment_samples = condition_comparison['treatment']

        # Find available samples
        available_control = [s for s in control_samples if s in pathway_expression.columns]
        available_treatment = [s for s in treatment_samples if s in pathway_expression.columns]

        if len(available_control) == 0 or len(available_treatment) == 0:
            return 0.0

        # Calculate mean expression for each condition
        control_mean = pathway_expression[available_control].mean(axis=1).mean()
        treatment_mean = pathway_expression[available_treatment].mean(axis=1).mean()

        # Calculate fold change
        if control_mean > 0:
            intervention_response = np.log2(treatment_mean / control_mean)
        else:
            intervention_response = 0.0

        return intervention_response

    def generate_pathway_report(self,
                              enrichment_results: List[PathwayEnrichmentResult],
                              aging_scores: List[AgingPathwayScore]) -> Dict[str, any]:
        """Generate comprehensive pathway analysis report."""

        report = {
            'enrichment_summary': {},
            'aging_pathways': {},
            'intervention_effects': {},
            'key_findings': []
        }

        # Enrichment summary
        significant_pathways = [r for r in enrichment_results if r.adjusted_p_value < 0.05]

        report['enrichment_summary'] = {
            'total_pathways_tested': len(enrichment_results),
            'significant_pathways': len(significant_pathways),
            'top_enriched': [r.pathway_name for r in significant_pathways[:10]]
        }

        # Aging pathway analysis
        highly_active = [s for s in aging_scores if s.activity_score > np.percentile([s.activity_score for s in aging_scores], 75)]

        report['aging_pathways'] = {
            'pathways_analyzed': len(aging_scores),
            'highly_active_pathways': [s.pathway_name for s in highly_active],
            'mean_confidence': np.mean([s.confidence_score for s in aging_scores])
        }

        # Intervention effects
        strong_responders = [s for s in aging_scores if abs(s.intervention_response) > 0.5]

        report['intervention_effects'] = {
            'pathways_with_response': len([s for s in aging_scores if s.intervention_response != 0]),
            'strong_responders': [s.pathway_name for s in strong_responders],
            'mean_response': np.mean([s.intervention_response for s in aging_scores])
        }

        # Key findings
        key_findings = self._generate_key_findings(enrichment_results, aging_scores)
        report['key_findings'] = key_findings

        return report

    def _generate_key_findings(self,
                             enrichment_results: List[PathwayEnrichmentResult],
                             aging_scores: List[AgingPathwayScore]) -> List[str]:
        """Generate key findings from pathway analysis."""

        findings = []

        # Top enriched pathways
        significant = [r for r in enrichment_results if r.adjusted_p_value < 0.05]
        if significant:
            findings.append(
                f"Found {len(significant)} significantly enriched pathways, "
                f"including {significant[0].pathway_name}"
            )

        # Aging pathway responses
        strong_responders = [s for s in aging_scores if abs(s.intervention_response) > 0.5]
        if strong_responders:
            findings.append(
                f"Strong intervention response in {len(strong_responders)} aging pathways"
            )

        # Specific pathway findings
        senescence_scores = [s for s in aging_scores if 'senescence' in s.pathway_name.lower()]
        if senescence_scores and senescence_scores[0].intervention_response < -0.3:
            findings.append("Cellular senescence pathway shows strong downregulation")

        autophagy_scores = [s for s in aging_scores if 'proteostasis' in s.pathway_name.lower()]
        if autophagy_scores and autophagy_scores[0].intervention_response > 0.3:
            findings.append("Proteostasis pathways show beneficial upregulation")

        return findings


def create_pathway_enrichment_demo():
    """Create demonstration of pathway enrichment analysis."""

    # Create sample gene list (aging-related genes)
    gene_list = [
        'TP53', 'CDKN1A', 'CDKN2A', 'SIRT1', 'FOXO3',
        'ATM', 'BRCA1', 'PGC1A', 'MTOR', 'BCL2'
    ]

    # Create sample expression data
    n_samples = 20
    sample_names = [f"Sample_{i}" for i in range(n_samples)]
    n_genes = 100

    # Generate expression data
    np.random.seed(42)
    expression_data = pd.DataFrame(
        np.random.lognormal(5, 1, (n_genes, n_samples)),
        index=[f"Gene_{i}" for i in range(n_genes)],
        columns=sample_names
    )

    # Add some known aging genes
    for i, gene in enumerate(gene_list[:10]):
        if i < len(expression_data):
            expression_data.index.values[i] = gene

    # Define conditions
    condition_comparison = {
        'control': sample_names[:10],
        'treatment': sample_names[10:]
    }

    # Initialize analyzer
    analyzer = PathwayEnrichmentAnalyzer()

    # Perform enrichment analysis
    enrichment_results = analyzer.perform_enrichment_analysis(gene_list)

    # Analyze aging pathways
    aging_scores = analyzer.analyze_aging_pathways(expression_data, condition_comparison)

    # Generate report
    report = analyzer.generate_pathway_report(enrichment_results, aging_scores)

    return analyzer, enrichment_results, aging_scores, report


if __name__ == "__main__":
    # Run demonstration
    analyzer, enrichment_results, aging_scores, report = create_pathway_enrichment_demo()

    print("Pathway enrichment analysis completed!")
    print(f"Enrichment results: {len(enrichment_results)}")
    print(f"Significant pathways: {report['enrichment_summary']['significant_pathways']}")
    print(f"Aging pathways analyzed: {len(aging_scores)}")
    print(f"Strong responders: {len(report['intervention_effects']['strong_responders'])}")
    print("Key findings:")
    for finding in report['key_findings']:
        print(f"  - {finding}")
