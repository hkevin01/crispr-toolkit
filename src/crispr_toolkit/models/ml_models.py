"""
Machine Learning models for aging intervention analysis.

This module contains ML models for predicting rejuvenation effectiveness
and intervention outcomes.
"""

import logging
import pickle
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


class RejuvenationModel:
    """ML model for predicting rejuvenation effectiveness."""

    def __init__(self, model_type: str = 'random_forest'):
        """Initialize rejuvenation model."""
        self.model_type = model_type
        self.model = None
        self.scaler = StandardScaler()
        self.is_trained = False
        self.feature_names = None

        # Initialize model based on type
        if model_type == 'random_forest':
            self.model = RandomForestRegressor(
                n_estimators=100,
                max_depth=10,
                random_state=42,
                n_jobs=-1
            )
        elif model_type == 'gradient_boosting':
            self.model = GradientBoostingRegressor(
                n_estimators=100,
                max_depth=6,
                learning_rate=0.1,
                random_state=42
            )
        else:
            raise ValueError(f"Unsupported model type: {model_type}")

    def train(self, X: np.ndarray, y: np.ndarray,
             feature_names: List[str] = None) -> Dict[str, float]:
        """Train the rejuvenation model."""
        logger.info(f"Training {self.model_type} rejuvenation model")

        # Store feature names
        if feature_names:
            self.feature_names = feature_names
        elif hasattr(X, 'columns'):
            self.feature_names = list(X.columns)
        else:
            self.feature_names = [f"feature_{i}" for i in range(X.shape[1])]

        # Convert to numpy if needed
        if hasattr(X, 'values'):
            X = X.values
        if hasattr(y, 'values'):
            y = y.values

        # Scale features
        X_scaled = self.scaler.fit_transform(X)

        # Split for training and validation
        X_train, X_val, y_train, y_val = train_test_split(
            X_scaled, y, test_size=0.2, random_state=42
        )

        # Train model
        self.model.fit(X_train, y_train)
        self.is_trained = True

        # Evaluate on validation set
        y_pred = self.model.predict(X_val)

        metrics = {
            'mse': mean_squared_error(y_val, y_pred),
            'rmse': np.sqrt(mean_squared_error(y_val, y_pred)),
            'mae': mean_absolute_error(y_val, y_pred),
            'r2': r2_score(y_val, y_pred)
        }

        # Cross-validation score
        cv_scores = cross_val_score(
            self.model, X_train, y_train, cv=5, scoring='r2'
        )
        metrics['cv_r2_mean'] = cv_scores.mean()
        metrics['cv_r2_std'] = cv_scores.std()

        logger.info(f"Model training completed. R² score: {metrics['r2']:.4f}")

        return metrics

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions with the trained model."""
        if not self.is_trained:
            raise ValueError("Model must be trained before making predictions")

        # Convert to numpy if needed
        if hasattr(X, 'values'):
            X = X.values

        # Scale features
        X_scaled = self.scaler.transform(X)

        return self.model.predict(X_scaled)

    def predict_proba(self, X: np.ndarray,
                     uncertainty_method: str = 'std') -> Tuple[np.ndarray, np.ndarray]:
        """
        Predict with uncertainty estimation.

        For tree-based models, use individual trees for uncertainty.
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before making predictions")

        # Convert to numpy if needed
        if hasattr(X, 'values'):
            X = X.values

        # Scale features
        X_scaled = self.scaler.transform(X)

        if hasattr(self.model, 'estimators_'):
            # For ensemble models, get predictions from each estimator
            predictions = np.array([
                tree.predict(X_scaled) for tree in self.model.estimators_
            ])

            mean_pred = predictions.mean(axis=0)
            std_pred = predictions.std(axis=0)

            return mean_pred, std_pred
        else:
            # For single models, return prediction with zero uncertainty
            pred = self.model.predict(X_scaled)
            return pred, np.zeros_like(pred)

    def get_feature_importance(self) -> Dict[str, float]:
        """Get feature importance scores."""
        if not self.is_trained:
            raise ValueError("Model must be trained before getting importance")

        if hasattr(self.model, 'feature_importances_'):
            importance_dict = {}
            for name, importance in zip(
                self.feature_names, self.model.feature_importances_
            ):
                importance_dict[name] = importance

            # Sort by importance
            return dict(sorted(
                importance_dict.items(),
                key=lambda x: x[1],
                reverse=True
            ))
        else:
            logger.warning("Model does not support feature importance")
            return {}

    def evaluate_performance(self, X: np.ndarray, y: np.ndarray) -> Dict[str, float]:
        """Evaluate model performance on given data."""
        if not self.is_trained:
            raise ValueError("Model must be trained before evaluation")

        # Make predictions
        y_pred = self.predict(X)

        # Calculate metrics
        metrics = {
            'mse': mean_squared_error(y, y_pred),
            'rmse': np.sqrt(mean_squared_error(y, y_pred)),
            'mae': mean_absolute_error(y, y_pred),
            'r2': r2_score(y, y_pred)
        }

        return metrics

    def save_model(self, filepath: str):
        """Save the trained model to disk."""
        if not self.is_trained:
            raise ValueError("Model must be trained before saving")

        model_data = {
            'model': self.model,
            'scaler': self.scaler,
            'model_type': self.model_type,
            'feature_names': self.feature_names,
            'is_trained': self.is_trained
        }

        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)

        logger.info(f"Model saved to {filepath}")

    def load_model(self, filepath: str):
        """Load a trained model from disk."""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)

        self.model = model_data['model']
        self.scaler = model_data['scaler']
        self.model_type = model_data['model_type']
        self.feature_names = model_data['feature_names']
        self.is_trained = model_data['is_trained']

        logger.info(f"Model loaded from {filepath}")


class InterventionPrioritizer:
    """Model for prioritizing aging interventions."""

    def __init__(self):
        """Initialize intervention prioritizer."""
        self.effectiveness_model = RejuvenationModel('random_forest')
        self.safety_model = RejuvenationModel('gradient_boosting')
        self.is_trained = False

    def train(self, training_data: pd.DataFrame) -> Dict[str, Any]:
        """Train both effectiveness and safety models."""
        logger.info("Training intervention prioritization models")

        # Prepare features
        feature_cols = [
            'n_targets', 'dose_level', 'duration_weeks',
            'tissue_encoded', 'intervention_encoded'
        ]

        # Encode categorical variables
        data = training_data.copy()

        # Simple encoding for demonstration
        tissue_map = {tissue: i for i, tissue in enumerate(data['tissue'].unique())}
        intervention_map = {
            intv: i for i, intv in enumerate(data['intervention_type'].unique())
        }

        data['tissue_encoded'] = data['tissue'].map(tissue_map)
        data['intervention_encoded'] = data['intervention_type'].map(intervention_map)

        X = data[feature_cols]

        # Train effectiveness model
        y_effectiveness = data['functional_improvement']
        effectiveness_metrics = self.effectiveness_model.train(
            X, y_effectiveness, feature_cols
        )

        # Train safety model
        y_safety = 1 - data['safety_risk']  # Convert risk to safety score
        safety_metrics = self.safety_model.train(X, y_safety, feature_cols)

        self.is_trained = True

        return {
            'effectiveness_metrics': effectiveness_metrics,
            'safety_metrics': safety_metrics,
            'feature_importance_effectiveness': (
                self.effectiveness_model.get_feature_importance()
            ),
            'feature_importance_safety': self.safety_model.get_feature_importance()
        }

    def prioritize_interventions(self, interventions: List[Dict]) -> List[Dict]:
        """Prioritize a list of interventions by predicted effectiveness and safety."""
        if not self.is_trained:
            raise ValueError("Models must be trained before prioritization")

        results = []

        for intervention in interventions:
            # Create feature vector (simplified)
            features = np.array([[
                intervention.get('n_targets', 1),
                intervention.get('dose_level', 0.5),
                intervention.get('duration_weeks', 4),
                intervention.get('tissue_encoded', 0),
                intervention.get('intervention_encoded', 0)
            ]])

            # Predict effectiveness and safety
            effectiveness = self.effectiveness_model.predict(features)[0]
            safety = self.safety_model.predict(features)[0]

            # Calculate combined priority score
            priority_score = (effectiveness * 0.6 + safety * 0.4)

            results.append({
                'intervention': intervention,
                'predicted_effectiveness': effectiveness,
                'predicted_safety': safety,
                'priority_score': priority_score
            })

        # Sort by priority score
        results.sort(key=lambda x: x['priority_score'], reverse=True)

        return results


def create_rejuvenation_model(model_type: str = 'random_forest') -> RejuvenationModel:
    """Convenience function to create a rejuvenation model."""
    return RejuvenationModel(model_type)
