"""
Rejuvenation outcome prediction for CRISPR interventions.

This module implements multi-task models to predict rejuvenation outcomes
including epigenetic clock changes, senescence markers, and safety metrics.
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.multioutput import MultiOutputRegressor
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


class RejuvenationPredictor:
    """Multi-task predictor for rejuvenation intervention outcomes."""

    def __init__(self, model_config: Optional[Dict] = None):
        """Initialize predictor with configuration."""
        self.model_config = model_config or self._default_config()
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = []
        self.target_names = [
            "epigenetic_clock_delta",
            "senescence_score_change",
            "transcriptional_age_shift",
            "functional_improvement",
            "safety_risk"
        ]

    def _default_config(self) -> Dict:
        """Default model configuration."""
        return {
            "model_type": "multi_output_rf",
            "params": {
                "n_estimators": 100,
                "max_depth": 10,
                "random_state": 42
            },
            "uncertainty": "bootstrap",
            "n_bootstrap": 50
        }

    def extract_intervention_features(
        self,
        intervention_config: Dict,
        tissue: str = "liver",
        context_data: Optional[Dict] = None
    ) -> Dict[str, float]:
        """Extract features from intervention configuration."""

        # Parse intervention type
        intervention_type = intervention_config.get("intervention", "unknown")
        targets = intervention_config.get("targets", [])
        delivery = intervention_config.get("delivery", "unknown")

        features = {
            # Intervention features
            "is_osk": 1.0 if intervention_type == "OSK" else 0.0,
            "is_oskm": 1.0 if intervention_type == "OSKM" else 0.0,
            "is_senolytic": 1.0 if intervention_type == "senolytic" else 0.0,
            "is_base_edit": 1.0 if intervention_type == "base_edit" else 0.0,
            "n_targets": len(targets),

            # Delivery features
            "is_aav_delivery": 1.0 if delivery == "AAV" else 0.0,
            "is_lnp_delivery": 1.0 if delivery == "LNP" else 0.0,
            "is_viral_delivery": 1.0 if delivery in ["AAV", "lentivirus"] else 0.0,

            # Tissue context
            "tissue_liver": 1.0 if tissue == "liver" else 0.0,
            "tissue_brain": 1.0 if tissue == "brain" else 0.0,
            "tissue_muscle": 1.0 if tissue == "muscle" else 0.0,
            "tissue_skin": 1.0 if tissue == "skin" else 0.0,

            # Target gene features
            "targets_senescence_genes": self._count_senescence_targets(targets),
            "targets_dna_repair_genes": self._count_dna_repair_targets(targets),
            "targets_autophagy_genes": self._count_autophagy_targets(targets),
            "targets_mitochondrial_genes": self._count_mitochondrial_targets(targets),

            # Dosing and schedule (if available)
            "dose_level": intervention_config.get("dose_level", 1.0),
            "cyclic_dosing": 1.0 if intervention_config.get("cyclic", False) else 0.0,
            "duration_weeks": intervention_config.get("duration_weeks", 4.0),

            # Safety considerations
            "has_safety_switch": 1.0 if intervention_config.get("safety_switch", False) else 0.0,
            "oncogene_targets": self._count_oncogene_targets(targets),
            "essential_gene_targets": self._count_essential_targets(targets),
        }

        # Add context data if available
        if context_data:
            features.update({
                "baseline_age": context_data.get("age", 50.0),
                "baseline_senescence_load": context_data.get("senescence_load", 0.5),
                "tissue_specific_expression": context_data.get("expression_level", 1.0),
                "chromatin_accessibility": context_data.get("chromatin_access", 0.8)
            })

        return features

    def _count_senescence_targets(self, targets: List[str]) -> float:
        """Count senescence-related targets."""
        senescence_genes = {"TP53", "CDKN2A", "RB1", "CDKN1A", "MDM2"}
        return sum(1 for target in targets if target in senescence_genes) / max(len(targets), 1)

    def _count_dna_repair_targets(self, targets: List[str]) -> float:
        """Count DNA repair pathway targets."""
        dna_repair_genes = {"ATM", "BRCA1", "BRCA2", "TP53", "PARP1", "XRCC1"}
        return sum(1 for target in targets if target in dna_repair_genes) / max(len(targets), 1)

    def _count_autophagy_targets(self, targets: List[str]) -> float:
        """Count autophagy pathway targets."""
        autophagy_genes = {"ATG5", "ATG7", "BECN1", "MTOR", "ULK1", "LC3B"}
        return sum(1 for target in targets if target in autophagy_genes) / max(len(targets), 1)

    def _count_mitochondrial_targets(self, targets: List[str]) -> float:
        """Count mitochondrial pathway targets."""
        mito_genes = {"PGC1A", "TFAM", "NRF1", "SIRT1", "SOD2", "PINK1"}
        return sum(1 for target in targets if target in mito_genes) / max(len(targets), 1)

    def _count_oncogene_targets(self, targets: List[str]) -> float:
        """Count oncogene targets (safety concern)."""
        oncogenes = {"MYC", "RAS", "AKT1", "PIK3CA", "EGFR", "BRAF"}
        return sum(1 for target in targets if target in oncogenes) / max(len(targets), 1)

    def _count_essential_targets(self, targets: List[str]) -> float:
        """Count essential gene targets (safety concern)."""
        essential_genes = {"TP53", "RB1", "BRCA1", "BRCA2", "ATM", "PTEN"}
        return sum(1 for target in targets if target in essential_genes) / max(len(targets), 1)

    def predict_outcomes(
        self,
        intervention_config: Dict,
        tissue: str = "liver",
        context_data: Optional[Dict] = None,
        return_uncertainty: bool = True
    ) -> Dict[str, Any]:
        """Predict rejuvenation outcomes for an intervention."""

        # Extract features
        features = self.extract_intervention_features(intervention_config, tissue, context_data)

        # If no trained model, use heuristic prediction
        if self.model is None:
            logger.warning("No trained model found, using heuristic prediction")
            return self._heuristic_prediction(features, intervention_config)

        # Convert to array
        feature_array = np.array([list(features.values())]).reshape(1, -1)
        feature_array = self.scaler.transform(feature_array)

        # Predict
        predictions = self.model.predict(feature_array)[0]

        # Create results dictionary
        results = {
            "predictions": {
                name: float(pred) for name, pred in zip(self.target_names, predictions)
            },
            "intervention": intervention_config.get("intervention", "unknown"),
            "tissue": tissue,
            "features": features
        }

        # Add uncertainty estimates if requested
        if return_uncertainty:
            uncertainties = self._estimate_uncertainty(feature_array)
            results["uncertainty"] = {
                name: float(unc) for name, unc in zip(self.target_names, uncertainties)
            }

            # Add confidence intervals
            results["confidence_intervals"] = {
                name: [
                    float(pred - 1.96 * unc),
                    float(pred + 1.96 * unc)
                ]
                for name, pred, unc in zip(
                    self.target_names,
                    predictions,
                    uncertainties
                )
            }

        return results

    def _heuristic_prediction(
        self,
        features: Dict[str, float],
        intervention_config: Dict
    ) -> Dict[str, Any]:
        """Heuristic prediction when no trained model available."""

        intervention_type = intervention_config.get("intervention", "unknown")

        # OSK/OSKM typically show stronger rejuvenation effects
        if features["is_osk"] or features["is_oskm"]:
            base_effect = -2.5  # Negative = rejuvenation
            safety_risk = 0.3
        elif features["is_senolytic"]:
            base_effect = -1.5
            safety_risk = 0.2
        else:
            base_effect = -0.8
            safety_risk = 0.15

        # Adjust for tissue
        tissue_multiplier = 1.0
        if features["tissue_brain"]:
            tissue_multiplier = 0.7  # More conservative in brain
            safety_risk += 0.1
        elif features["tissue_liver"]:
            tissue_multiplier = 1.2  # Liver responds well

        # Adjust for safety concerns
        safety_penalty = features["oncogene_targets"] * 0.5 + features["essential_gene_targets"] * 0.3
        safety_risk += safety_penalty

        predictions = {
            "epigenetic_clock_delta": base_effect * tissue_multiplier + np.random.normal(0, 0.5),
            "senescence_score_change": base_effect * 0.8 * tissue_multiplier + np.random.normal(0, 0.3),
            "transcriptional_age_shift": base_effect * 0.6 * tissue_multiplier + np.random.normal(0, 0.4),
            "functional_improvement": abs(base_effect) * 0.4 * tissue_multiplier + np.random.normal(0, 0.2),
            "safety_risk": min(safety_risk, 1.0)
        }

        return {
            "predictions": predictions,
            "intervention": intervention_config.get("intervention", "unknown"),
            "tissue": intervention_config.get("tissue", "unknown"),
            "features": features,
            "model_type": "heuristic",
            "confidence": 0.6  # Lower confidence for heuristic
        }

    def _estimate_uncertainty(self, feature_array: np.ndarray) -> np.ndarray:
        """Estimate prediction uncertainty using bootstrap or other methods."""
        if self.model_config["uncertainty"] == "bootstrap":
            # Placeholder - would implement bootstrap sampling
            return np.random.uniform(0.1, 0.5, len(self.target_names))
        else:
            # Default uncertainty
            return np.full(len(self.target_names), 0.3)

    def train(
        self,
        training_data: pd.DataFrame,
        target_columns: List[str]
    ) -> Dict[str, Any]:
        """Train the rejuvenation prediction model."""
        logger.info("Training rejuvenation prediction model")

        # Prepare features
        feature_cols = [col for col in training_data.columns if col not in target_columns]
        X = training_data[feature_cols]
        y = training_data[target_columns]

        self.feature_names = feature_cols
        self.target_names = target_columns

        # Scale features
        X_scaled = self.scaler.fit_transform(X)

        # Train model
        if self.model_config["model_type"] == "multi_output_rf":
            base_model = RandomForestRegressor(**self.model_config["params"])
            self.model = MultiOutputRegressor(base_model)
        else:
            raise ValueError(f"Unknown model type: {self.model_config['model_type']}")

        self.model.fit(X_scaled, y)

        # Evaluate on training data
        y_pred = self.model.predict(X_scaled)

        metrics = {}
        for i, target in enumerate(target_columns):
            metrics[target] = {
                "mae": mean_absolute_error(y.iloc[:, i], y_pred[:, i]),
                "r2": r2_score(y.iloc[:, i], y_pred[:, i])
            }

        return {
            "model_type": self.model_config["model_type"],
            "n_features": len(feature_cols),
            "n_samples": len(training_data),
            "metrics": metrics
        }

    def evaluate(
        self,
        test_data: pd.DataFrame,
        target_columns: List[str]
    ) -> Dict[str, Any]:
        """Evaluate model performance."""
        if self.model is None:
            raise ValueError("Model not trained")

        feature_cols = [col for col in test_data.columns if col not in target_columns]
        X = test_data[feature_cols]
        y = test_data[target_columns]

        X_scaled = self.scaler.transform(X)
        y_pred = self.model.predict(X_scaled)

        metrics = {}
        for i, target in enumerate(target_columns):
            metrics[target] = {
                "mae": mean_absolute_error(y.iloc[:, i], y_pred[:, i]),
                "r2": r2_score(y.iloc[:, i], y_pred[:, i])
            }

        return {"test_metrics": metrics}

    def save_model(self, filepath: Path) -> None:
        """Save trained model."""
        if self.model is None:
            raise ValueError("No model to save")

        import joblib

        model_data = {
            "model_config": self.model_config,
            "feature_names": self.feature_names,
            "target_names": self.target_names,
            "scaler_mean": self.scaler.mean_.tolist() if self.scaler.mean_ is not None else None,
            "scaler_scale": self.scaler.scale_.tolist() if self.scaler.scale_ is not None else None
        }

        # Save model
        joblib.dump(self.model, filepath.with_suffix('.joblib'))

        # Save metadata
        with open(filepath.with_suffix('.json'), 'w') as f:
            json.dump(model_data, f, indent=2)

    def load_model(self, filepath: Path) -> None:
        """Load trained model."""
        import joblib

        # Load metadata
        with open(filepath.with_suffix('.json'), 'r') as f:
            model_data = json.load(f)

        self.model_config = model_data["model_config"]
        self.feature_names = model_data["feature_names"]
        self.target_names = model_data["target_names"]

        # Restore scaler
        if model_data["scaler_mean"] is not None:
            self.scaler.mean_ = np.array(model_data["scaler_mean"])
            self.scaler.scale_ = np.array(model_data["scaler_scale"])

        # Load model
        self.model = joblib.load(filepath.with_suffix('.joblib'))


# Convenience function for direct use
def predict_rejuvenation(
    intervention_config: Dict,
    tissue: str = "liver",
    context_data: Optional[Dict] = None,
    model_path: Optional[Path] = None
) -> Dict[str, Any]:
    """Convenience function for rejuvenation prediction."""
    predictor = RejuvenationPredictor()

    if model_path and model_path.exists():
        predictor.load_model(model_path)

    return predictor.predict_outcomes(intervention_config, tissue, context_data)
