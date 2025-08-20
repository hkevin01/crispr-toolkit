"""
Multi-modal target prioritization for aging interventions.

This module implements learning-to-rank models that combine literature evidence,
network centrality, expression patterns, and CRISPR feasibility to prioritize
gene targets for aging research applications.
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


class AgingTargetPrioritizer:
    """Multi-modal target prioritization for aging interventions."""

    def __init__(self, model_config: Optional[Dict] = None):
        """Initialize prioritizer with configuration."""
        self.model_config = model_config or self._default_config()
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = []

    def _default_config(self) -> Dict:
        """Default model configuration."""
        return {
            "model_type": "lightgbm_ltr",
            "params": {
                "objective": "lambdarank",
                "metric": "ndcg",
                "ndcg_eval_at": [5, 10, 20],
                "num_leaves": 50,
                "learning_rate": 0.1,
                "feature_fraction": 0.9,
                "verbosity": -1
            },
            "features": [
                "expression_tissue_specificity",
                "ppi_degree",
                "literature_evidence_score",
                "editability_score",
                "prior_screen_logfc",
                "aging_pathway_membership",
                "safety_score"
            ]
        }

    def compute_features(
        self,
        genes: List[str],
        tissue: str = "liver",
        data_sources: Optional[Dict] = None
    ) -> pd.DataFrame:
        """Compute multi-modal features for gene targets."""
        logger.info(f"Computing features for {len(genes)} genes in {tissue}")

        features = []
        for gene in genes:
            feature_dict = {
                "gene": gene,
                "tissue": tissue,
                # Expression features
                "expression_tissue_specificity": self._compute_expression_specificity(gene, tissue),
                "expression_age_gradient": self._compute_age_gradient(gene, tissue),

                # Network features
                "ppi_degree": self._compute_ppi_degree(gene),
                "aging_pathway_membership": self._compute_pathway_membership(gene),

                # Literature features
                "literature_evidence_score": self._compute_literature_score(gene),
                "publication_count": self._compute_publication_count(gene),

                # CRISPR features
                "editability_score": self._compute_editability(gene),
                "pam_availability": self._compute_pam_availability(gene),

                # Prior evidence
                "prior_screen_logfc": self._compute_prior_screen_evidence(gene, tissue),
                "genage_evidence": self._compute_genage_evidence(gene),

                # Safety features
                "safety_score": self._compute_safety_score(gene),
                "essentiality_score": self._compute_essentiality(gene)
            }
            features.append(feature_dict)

        return pd.DataFrame(features)

    def _compute_expression_specificity(self, gene: str, tissue: str) -> float:
        """Compute tissue-specific expression score."""
        # Placeholder - integrate with GTEx/HPA data
        return np.random.uniform(0.3, 0.9)

    def _compute_age_gradient(self, gene: str, tissue: str) -> float:
        """Compute age-related expression changes."""
        # Placeholder - integrate with age-stratified expression data
        return np.random.uniform(-2.0, 2.0)

    def _compute_ppi_degree(self, gene: str) -> float:
        """Compute protein-protein interaction network degree."""
        # Placeholder - integrate with STRING/BioGRID
        return np.random.uniform(1, 100)

    def _compute_pathway_membership(self, gene: str) -> float:
        """Compute aging pathway membership score."""
        # Key aging pathways
        aging_pathways = {
            "senescence": ["TP53", "CDKN2A", "RB1", "CDKN1A"],
            "dna_repair": ["ATM", "BRCA1", "BRCA2", "TP53"],
            "autophagy": ["ATG5", "ATG7", "BECN1", "MTOR"],
            "mitochondrial": ["PGC1A", "TFAM", "NRF1", "SIRT1"],
            "inflammation": ["TNF", "IL6", "NFKB1", "STAT3"],
            "proteostasis": ["HSF1", "HSPA1A", "HSPB1", "UBE3A"]
        }

        score = 0.0
        for pathway, pathway_genes in aging_pathways.items():
            if gene in pathway_genes:
                score += 1.0

        return min(score / len(aging_pathways), 1.0)

    def _compute_literature_score(self, gene: str) -> float:
        """Compute literature evidence score for aging relevance."""
        # Placeholder - would integrate with PubMed mining
        aging_keywords = ["aging", "senescence", "longevity", "rejuvenation"]

        # Known aging genes get higher scores
        known_aging_genes = [
            "TP53", "CDKN2A", "TERT", "SIRT1", "FOXO3", "IGF1R",
            "MTOR", "ATG7", "PGC1A", "HSPA1A", "TNF", "IL6"
        ]

        if gene in known_aging_genes:
            return np.random.uniform(0.7, 0.95)
        else:
            return np.random.uniform(0.1, 0.6)

    def _compute_publication_count(self, gene: str) -> int:
        """Count publications mentioning gene and aging."""
        # Placeholder
        return np.random.randint(1, 500)

    def _compute_editability(self, gene: str) -> float:
        """Compute CRISPR editability score."""
        # Placeholder - would integrate PAM site analysis
        return np.random.uniform(0.4, 0.95)

    def _compute_pam_availability(self, gene: str) -> float:
        """Compute PAM site availability score."""
        # Placeholder
        return np.random.uniform(0.3, 1.0)

    def _compute_prior_screen_evidence(self, gene: str, tissue: str) -> float:
        """Compute evidence from prior CRISPR screens."""
        # Placeholder - would integrate screen databases
        return np.random.uniform(-3.0, 3.0)

    def _compute_genage_evidence(self, gene: str) -> float:
        """Compute GenAge database evidence."""
        # Known longevity genes
        genage_genes = [
            "FOXO3", "APOE", "SIRT1", "TP53", "IGF1R", "TERT",
            "CDKN2A", "WRN", "LMNA", "GHR"
        ]
        return 1.0 if gene in genage_genes else 0.0

    def _compute_safety_score(self, gene: str) -> float:
        """Compute safety score (higher = safer)."""
        # Essential genes and oncogenes get lower safety scores
        essential_genes = ["TP53", "RB1", "BRCA1", "BRCA2", "ATM"]
        oncogenes = ["MYC", "RAS", "AKT1", "PIK3CA", "EGFR"]

        if gene in essential_genes or gene in oncogenes:
            return np.random.uniform(0.2, 0.6)
        else:
            return np.random.uniform(0.6, 0.9)

    def _compute_essentiality(self, gene: str) -> float:
        """Compute gene essentiality score."""
        # Placeholder - would integrate DepMap data
        return np.random.uniform(0.0, 1.0)

    def train(
        self,
        training_data: pd.DataFrame,
        labels: np.ndarray,
        groups: Optional[np.ndarray] = None
    ) -> Dict[str, Any]:
        """Train the prioritization model."""
        logger.info("Training target prioritization model")

        # Prepare features
        feature_cols = [col for col in training_data.columns
                       if col not in ['gene', 'tissue']]
        X = training_data[feature_cols]
        self.feature_names = feature_cols

        # Scale features
        X_scaled = self.scaler.fit_transform(X)

        if self.model_config["model_type"] == "lightgbm_ltr":
            # LightGBM Learning-to-Rank
            train_data = lgb.Dataset(X_scaled, label=labels, group=groups)
            self.model = lgb.train(
                self.model_config["params"],
                train_data,
                num_boost_round=100,
                valid_sets=[train_data]
            )
        else:
            # Fallback to RandomForest
            self.model = RandomForestRegressor(
                n_estimators=100,
                random_state=42
            )
            self.model.fit(X_scaled, labels)

        # Compute feature importance
        if hasattr(self.model, 'feature_importance'):
            importance = self.model.feature_importance()
        else:
            importance = self.model.feature_importances_

        feature_importance = dict(zip(self.feature_names, importance))

        return {
            "model_type": self.model_config["model_type"],
            "feature_importance": feature_importance,
            "n_features": len(feature_cols),
            "n_samples": len(training_data)
        }

    def predict(self, features: pd.DataFrame) -> pd.DataFrame:
        """Predict target priorities."""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")

        # Prepare features
        feature_cols = [col for col in features.columns
                       if col not in ['gene', 'tissue']]
        X = features[feature_cols]
        X_scaled = self.scaler.transform(X)

        # Predict scores
        scores = self.model.predict(X_scaled)

        # Create results dataframe
        results = features[['gene', 'tissue']].copy()
        results['priority_score'] = scores
        results['rank'] = results['priority_score'].rank(ascending=False)

        # Add confidence intervals (placeholder)
        results['confidence_lower'] = scores - 0.1
        results['confidence_upper'] = scores + 0.1

        return results.sort_values('priority_score', ascending=False)

    def rank_targets(
        self,
        genes: List[str],
        tissue: str = "liver",
        phenotype: str = "senescence",
        n_targets: int = 20
    ) -> List[Dict]:
        """Main interface for target ranking."""
        logger.info(f"Ranking {len(genes)} targets for {tissue} {phenotype}")

        # Compute features
        features_df = self.compute_features(genes, tissue)

        # If no trained model, use heuristic scoring
        if self.model is None:
            logger.warning("No trained model found, using heuristic scoring")
            return self._heuristic_ranking(features_df, n_targets)

        # Predict with trained model
        results = self.predict(features_df)

        # Format output
        ranked_targets = []
        for _, row in results.head(n_targets).iterrows():
            target = {
                "gene": row["gene"],
                "tissue": row["tissue"],
                "priority_score": float(row["priority_score"]),
                "rank": int(row["rank"]),
                "confidence_interval": [
                    float(row["confidence_lower"]),
                    float(row["confidence_upper"])
                ],
                "features": {
                    col: float(features_df[features_df["gene"] == row["gene"]][col].iloc[0])
                    for col in features_df.columns
                    if col not in ["gene", "tissue"]
                },
                "rationale": self._generate_rationale(row, features_df)
            }
            ranked_targets.append(target)

        return ranked_targets

    def _heuristic_ranking(self, features_df: pd.DataFrame, n_targets: int) -> List[Dict]:
        """Heuristic ranking when no trained model is available."""
        # Simple weighted combination
        weights = {
            "literature_evidence_score": 0.3,
            "aging_pathway_membership": 0.25,
            "safety_score": 0.2,
            "editability_score": 0.15,
            "genage_evidence": 0.1
        }

        scores = np.zeros(len(features_df))
        for feature, weight in weights.items():
            if feature in features_df.columns:
                scores += weight * features_df[feature].values

        features_df["heuristic_score"] = scores
        features_df["rank"] = features_df["heuristic_score"].rank(ascending=False)

        ranked_targets = []
        for _, row in features_df.nlargest(n_targets, "heuristic_score").iterrows():
            target = {
                "gene": row["gene"],
                "tissue": row["tissue"],
                "priority_score": float(row["heuristic_score"]),
                "rank": int(row["rank"]),
                "confidence_interval": [
                    float(row["heuristic_score"] - 0.1),
                    float(row["heuristic_score"] + 0.1)
                ],
                "rationale": "High aging pathway relevance and literature evidence"
            }
            ranked_targets.append(target)

        return ranked_targets

    def _generate_rationale(self, row: pd.Series, features_df: pd.DataFrame) -> str:
        """Generate explanation for ranking."""
        gene = row["gene"]
        gene_features = features_df[features_df["gene"] == gene].iloc[0]

        reasons = []
        if gene_features["literature_evidence_score"] > 0.7:
            reasons.append("strong literature evidence")
        if gene_features["aging_pathway_membership"] > 0.5:
            reasons.append("key aging pathway involvement")
        if gene_features["safety_score"] > 0.7:
            reasons.append("favorable safety profile")
        if gene_features["editability_score"] > 0.8:
            reasons.append("high CRISPR editability")

        if not reasons:
            reasons = ["moderate evidence across multiple factors"]

        return f"Prioritized due to {', '.join(reasons)}"

    def save_model(self, filepath: Path) -> None:
        """Save trained model."""
        if self.model is None:
            raise ValueError("No model to save")

        model_data = {
            "model_config": self.model_config,
            "feature_names": self.feature_names,
            "scaler_mean": self.scaler.mean_.tolist(),
            "scaler_scale": self.scaler.scale_.tolist()
        }

        # Save model-specific data
        if self.model_config["model_type"] == "lightgbm_ltr":
            self.model.save_model(str(filepath.with_suffix('.lgb')))
        else:
            import joblib
            joblib.dump(self.model, filepath.with_suffix('.joblib'))

        # Save metadata
        with open(filepath.with_suffix('.json'), 'w') as f:
            json.dump(model_data, f, indent=2)

    def load_model(self, filepath: Path) -> None:
        """Load trained model."""
        # Load metadata
        with open(filepath.with_suffix('.json'), 'r') as f:
            model_data = json.load(f)

        self.model_config = model_data["model_config"]
        self.feature_names = model_data["feature_names"]

        # Restore scaler
        self.scaler.mean_ = np.array(model_data["scaler_mean"])
        self.scaler.scale_ = np.array(model_data["scaler_scale"])

        # Load model
        if self.model_config["model_type"] == "lightgbm_ltr":
            self.model = lgb.Booster(model_file=str(filepath.with_suffix('.lgb')))
        else:
            import joblib
            self.model = joblib.load(filepath.with_suffix('.joblib'))


# Convenience function for direct use
def rank_targets(
    genes: List[str],
    tissue: str = "liver",
    phenotype: str = "senescence",
    n_targets: int = 20,
    model_path: Optional[Path] = None
) -> List[Dict]:
    """Convenience function for target ranking."""
    prioritizer = AgingTargetPrioritizer()

    if model_path and model_path.exists():
        prioritizer.load_model(model_path)

    return prioritizer.rank_targets(genes, tissue, phenotype, n_targets)
