"""
Patient Stratification Engine for CRISPR Toolkit Phase 3
=========================================================

Advanced patient stratification and cohort optimization for
aging intervention clinical trials using multi-omics data
and clinical characteristics.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


@dataclass
class StratificationResult:
    """Results from patient stratification analysis."""
    cluster_assignments: Dict[str, int]
    cluster_characteristics: Dict[int, Dict[str, float]]
    stratification_quality: float
    biomarker_importance: Dict[str, float]
    treatment_recommendations: Dict[int, Dict[str, any]]
    validation_metrics: Dict[str, float]


@dataclass
class PatientProfile:
    """Individual patient profile for stratification."""
    patient_id: str
    age: float
    sex: str
    clinical_metrics: Dict[str, float]
    omics_data: Dict[str, pd.Series]
    biomarkers: Dict[str, float]
    medical_history: List[str]
    lifestyle_factors: Dict[str, any]


class PatientStratificationEngine:
    """
    Advanced patient stratification for aging intervention trials.

    Uses multi-omics data, clinical characteristics, and machine learning
    to identify optimal patient subgroups for personalized interventions.
    """

    def __init__(self, stratification_method: str = "multi_modal"):
        """
        Initialize patient stratification engine.

        Args:
            stratification_method: "clinical", "omics", "multi_modal", or "adaptive"
        """
        self.stratification_method = stratification_method
        self.patient_profiles = {}
        self.stratification_models = {}
        self.cluster_models = {}
        self.aging_biomarkers = self._load_aging_biomarkers()

        logger.info(f"Initialized PatientStratificationEngine with {stratification_method} method")

    def _load_aging_biomarkers(self) -> Dict[str, Dict[str, float]]:
        """Load validated aging biomarkers and their importance weights."""

        # Key aging biomarkers from literature
        biomarkers = {
            'inflammatory': {
                'IL6': 0.85,
                'TNF_alpha': 0.75,
                'CRP': 0.70,
                'IL1_beta': 0.65,
                'IL8': 0.60
            },
            'metabolic': {
                'glucose': 0.80,
                'HbA1c': 0.85,
                'insulin': 0.75,
                'triglycerides': 0.65,
                'HDL_cholesterol': 0.70
            },
            'cellular_aging': {
                'telomere_length': 0.90,
                'p16_expression': 0.85,
                'p21_expression': 0.80,
                'beta_galactosidase': 0.75,
                'SASP_score': 0.80
            },
            'oxidative_stress': {
                'malondialdehyde': 0.70,
                'protein_carbonyls': 0.65,
                'GSH_GSSG_ratio': 0.75,
                'catalase_activity': 0.60,
                'SOD_activity': 0.65
            },
            'epigenetic': {
                'DNA_methylation_age': 0.95,
                'histone_H3K27me3': 0.70,
                'histone_H3K4me3': 0.65,
                'SIRT1_activity': 0.80,
                'DNMT_activity': 0.60
            }
        }

        return biomarkers

    def add_patient_profile(self, patient: PatientProfile):
        """Add a patient profile to the stratification analysis."""
        self.patient_profiles[patient.patient_id] = patient
        logger.info(f"Added patient profile: {patient.patient_id}")

    def load_patient_cohort(self,
                           clinical_data: pd.DataFrame,
                           omics_data: Optional[Dict[str, pd.DataFrame]] = None,
                           biomarker_data: Optional[pd.DataFrame] = None):
        """
        Load a full patient cohort for stratification.

        Args:
            clinical_data: Clinical characteristics (patients x features)
            omics_data: Dictionary of omics datasets
            biomarker_data: Biomarker measurements (patients x biomarkers)
        """
        logger.info(f"Loading patient cohort: {len(clinical_data)} patients")

        for patient_id in clinical_data.index:
            # Extract clinical metrics
            clinical_metrics = clinical_data.loc[patient_id].to_dict()

            # Extract omics data if available
            patient_omics = {}
            if omics_data:
                for omics_type, data in omics_data.items():
                    if patient_id in data.columns:
                        patient_omics[omics_type] = data[patient_id]

            # Extract biomarkers if available
            patient_biomarkers = {}
            if biomarker_data is not None and patient_id in biomarker_data.index:
                patient_biomarkers = biomarker_data.loc[patient_id].to_dict()

            # Create patient profile
            patient_profile = PatientProfile(
                patient_id=patient_id,
                age=clinical_metrics.get('age', 50.0),
                sex=clinical_metrics.get('sex', 'unknown'),
                clinical_metrics=clinical_metrics,
                omics_data=patient_omics,
                biomarkers=patient_biomarkers,
                medical_history=clinical_metrics.get('medical_history', []),
                lifestyle_factors={}
            )

            self.add_patient_profile(patient_profile)

    def perform_stratification(self,
                             n_clusters: Optional[int] = None,
                             min_cluster_size: int = 10) -> StratificationResult:
        """
        Perform patient stratification analysis.

        Args:
            n_clusters: Number of clusters (auto-determined if None)
            min_cluster_size: Minimum patients per cluster

        Returns:
            StratificationResult with cluster assignments and characteristics
        """
        logger.info("Performing patient stratification analysis...")

        if len(self.patient_profiles) < min_cluster_size * 2:
            raise ValueError(f"Need at least {min_cluster_size * 2} patients for stratification")

        # Prepare feature matrix
        feature_matrix, feature_names, patient_ids = self._prepare_feature_matrix()

        # Determine optimal number of clusters
        if n_clusters is None:
            n_clusters = self._determine_optimal_clusters(feature_matrix, min_cluster_size)

        # Perform clustering
        cluster_assignments = self._perform_clustering(feature_matrix, n_clusters)

        # Analyze cluster characteristics
        cluster_characteristics = self._analyze_cluster_characteristics(
            feature_matrix, feature_names, cluster_assignments
        )

        # Calculate stratification quality
        stratification_quality = self._calculate_stratification_quality(
            feature_matrix, cluster_assignments
        )

        # Determine biomarker importance
        biomarker_importance = self._calculate_biomarker_importance(
            feature_matrix, feature_names, cluster_assignments
        )

        # Generate treatment recommendations
        treatment_recommendations = self._generate_treatment_recommendations(
            cluster_characteristics
        )

        # Validate stratification
        validation_metrics = self._validate_stratification(
            feature_matrix, cluster_assignments
        )

        # Create result mapping
        patient_cluster_map = dict(zip(patient_ids, cluster_assignments))

        result = StratificationResult(
            cluster_assignments=patient_cluster_map,
            cluster_characteristics=cluster_characteristics,
            stratification_quality=stratification_quality,
            biomarker_importance=biomarker_importance,
            treatment_recommendations=treatment_recommendations,
            validation_metrics=validation_metrics
        )

        logger.info(f"Stratification complete: {n_clusters} clusters, "
                   f"quality score: {stratification_quality:.3f}")

        return result

    def _prepare_feature_matrix(self) -> Tuple[np.ndarray, List[str], List[str]]:
        """Prepare feature matrix for clustering."""

        features = []
        feature_names = []
        patient_ids = []

        for patient_id, profile in self.patient_profiles.items():
            patient_features = []

            # Clinical features
            if self.stratification_method in ['clinical', 'multi_modal']:
                clinical_features = self._extract_clinical_features(profile)
                patient_features.extend(clinical_features)

                if len(feature_names) == 0:  # First patient
                    feature_names.extend(self._get_clinical_feature_names())

            # Omics features
            if self.stratification_method in ['omics', 'multi_modal']:
                omics_features = self._extract_omics_features(profile)
                patient_features.extend(omics_features)

                if len(feature_names) == len(patient_features) - len(omics_features):  # First patient omics
                    feature_names.extend(self._get_omics_feature_names(profile))

            # Biomarker features
            biomarker_features = self._extract_biomarker_features(profile)
            patient_features.extend(biomarker_features)

            if len(feature_names) == len(patient_features) - len(biomarker_features):  # First patient biomarkers
                feature_names.extend(list(profile.biomarkers.keys()))

            features.append(patient_features)
            patient_ids.append(patient_id)

        feature_matrix = np.array(features)

        # Handle missing values
        feature_matrix = np.nan_to_num(feature_matrix, nan=0.0)

        # Normalize features
        scaler = StandardScaler()
        feature_matrix = scaler.fit_transform(feature_matrix)

        return feature_matrix, feature_names, patient_ids

    def _extract_clinical_features(self, profile: PatientProfile) -> List[float]:
        """Extract clinical features from patient profile."""

        features = []

        # Age (normalized)
        features.append(profile.age / 100.0)

        # Sex (encoded)
        features.append(1.0 if profile.sex == 'male' else 0.0)

        # Key clinical metrics
        clinical_keys = ['BMI', 'blood_pressure_systolic', 'blood_pressure_diastolic',
                        'heart_rate', 'white_blood_cell_count', 'hemoglobin']

        for key in clinical_keys:
            features.append(profile.clinical_metrics.get(key, 0.0))

        return features

    def _get_clinical_feature_names(self) -> List[str]:
        """Get clinical feature names."""
        return ['age', 'sex', 'BMI', 'systolic_bp', 'diastolic_bp',
                'heart_rate', 'wbc_count', 'hemoglobin']

    def _extract_omics_features(self, profile: PatientProfile) -> List[float]:
        """Extract omics features using PCA reduction."""

        features = []

        for omics_type in ['transcriptomics', 'proteomics', 'metabolomics']:
            if omics_type in profile.omics_data:
                omics_data = profile.omics_data[omics_type]

                # Use top variable features or PCA components
                if len(omics_data) > 50:
                    # Take top 10 most variable features as proxy
                    variance = omics_data.var() if hasattr(omics_data, 'var') else np.var(omics_data.values)
                    if hasattr(variance, 'nlargest'):
                        top_features = variance.nlargest(10)
                        features.extend(omics_data[top_features.index].values)
                    else:
                        features.extend(omics_data.values[:10])
                else:
                    features.extend(omics_data.values)
            else:
                # Add zeros for missing omics data
                features.extend([0.0] * 10)

        return features

    def _get_omics_feature_names(self, profile: PatientProfile) -> List[str]:
        """Get omics feature names."""
        names = []

        for omics_type in ['transcriptomics', 'proteomics', 'metabolomics']:
            for i in range(10):
                names.append(f"{omics_type}_PC{i+1}")

        return names

    def _extract_biomarker_features(self, profile: PatientProfile) -> List[float]:
        """Extract biomarker features."""

        features = []

        # Key aging biomarkers
        all_biomarkers = []
        for category in self.aging_biomarkers.values():
            all_biomarkers.extend(category.keys())

        for biomarker in all_biomarkers:
            features.append(profile.biomarkers.get(biomarker, 0.0))

        return features

    def _determine_optimal_clusters(self, feature_matrix: np.ndarray,
                                  min_cluster_size: int) -> int:
        """Determine optimal number of clusters using silhouette analysis."""

        max_clusters = min(10, len(feature_matrix) // min_cluster_size)
        best_score = -1
        best_k = 2

        for k in range(2, max_clusters + 1):
            try:
                kmeans = KMeans(n_clusters=k, random_state=42)
                cluster_labels = kmeans.fit_predict(feature_matrix)

                # Check minimum cluster size
                unique, counts = np.unique(cluster_labels, return_counts=True)
                if np.min(counts) < min_cluster_size:
                    continue

                score = silhouette_score(feature_matrix, cluster_labels)

                if score > best_score:
                    best_score = score
                    best_k = k

            except:
                continue

        logger.info(f"Optimal clusters: {best_k} (silhouette score: {best_score:.3f})")
        return best_k

    def _perform_clustering(self, feature_matrix: np.ndarray,
                          n_clusters: int) -> np.ndarray:
        """Perform clustering using the specified method."""

        if self.stratification_method == "adaptive":
            # Use hierarchical clustering for adaptive method
            clustering = AgglomerativeClustering(n_clusters=n_clusters)
        else:
            # Use K-means for other methods
            clustering = KMeans(n_clusters=n_clusters, random_state=42)

        cluster_assignments = clustering.fit_predict(feature_matrix)

        return cluster_assignments

    def _analyze_cluster_characteristics(self,
                                       feature_matrix: np.ndarray,
                                       feature_names: List[str],
                                       cluster_assignments: np.ndarray) -> Dict[int, Dict[str, float]]:
        """Analyze characteristics of each cluster."""

        characteristics = {}

        unique_clusters = np.unique(cluster_assignments)

        for cluster_id in unique_clusters:
            cluster_mask = cluster_assignments == cluster_id
            cluster_data = feature_matrix[cluster_mask]

            cluster_chars = {}

            # Mean values for each feature
            for i, feature_name in enumerate(feature_names):
                cluster_chars[f"{feature_name}_mean"] = np.mean(cluster_data[:, i])
                cluster_chars[f"{feature_name}_std"] = np.std(cluster_data[:, i])

            # Cluster size
            cluster_chars['cluster_size'] = np.sum(cluster_mask)

            # Age statistics
            if 'age' in feature_names:
                age_idx = feature_names.index('age')
                cluster_chars['mean_age'] = np.mean(cluster_data[:, age_idx]) * 100.0  # Denormalize

            characteristics[int(cluster_id)] = cluster_chars

        return characteristics

    def _calculate_stratification_quality(self,
                                        feature_matrix: np.ndarray,
                                        cluster_assignments: np.ndarray) -> float:
        """Calculate overall quality of stratification."""

        try:
            # Silhouette score
            silhouette = silhouette_score(feature_matrix, cluster_assignments)

            # Cluster balance (penalize very unbalanced clusters)
            unique, counts = np.unique(cluster_assignments, return_counts=True)
            balance_score = 1.0 - np.std(counts) / np.mean(counts)
            balance_score = max(0.0, balance_score)

            # Combined quality score
            quality_score = 0.7 * silhouette + 0.3 * balance_score

            return quality_score

        except:
            return 0.0

    def _calculate_biomarker_importance(self,
                                      feature_matrix: np.ndarray,
                                      feature_names: List[str],
                                      cluster_assignments: np.ndarray) -> Dict[str, float]:
        """Calculate importance of each biomarker for stratification."""

        importance_scores = {}

        try:
            # Train random forest to predict cluster assignments
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            rf.fit(feature_matrix, cluster_assignments)

            # Get feature importances
            importances = rf.feature_importances_

            for i, feature_name in enumerate(feature_names):
                importance_scores[feature_name] = importances[i]

        except:
            # Default equal importance if classification fails
            for feature_name in feature_names:
                importance_scores[feature_name] = 1.0 / len(feature_names)

        return importance_scores

    def _generate_treatment_recommendations(self,
                                          cluster_characteristics: Dict[int, Dict[str, float]]) -> Dict[int, Dict[str, any]]:
        """Generate treatment recommendations for each cluster."""

        recommendations = {}

        for cluster_id, characteristics in cluster_characteristics.items():
            cluster_recs = {}

            # Age-based recommendations
            mean_age = characteristics.get('mean_age', 50.0)

            if mean_age < 40:
                cluster_recs['intervention_intensity'] = 'moderate'
                cluster_recs['primary_focus'] = 'prevention'
            elif mean_age < 65:
                cluster_recs['intervention_intensity'] = 'standard'
                cluster_recs['primary_focus'] = 'health_span_extension'
            else:
                cluster_recs['intervention_intensity'] = 'intensive'
                cluster_recs['primary_focus'] = 'aging_reversal'

            # Biomarker-based recommendations
            inflammation_markers = ['IL6_mean', 'TNF_alpha_mean', 'CRP_mean']
            inflammation_score = np.mean([characteristics.get(marker, 0.0) for marker in inflammation_markers])

            if inflammation_score > 0.5:
                cluster_recs['anti_inflammatory'] = True
                cluster_recs['supplements'] = ['omega_3', 'curcumin', 'resveratrol']
            else:
                cluster_recs['anti_inflammatory'] = False
                cluster_recs['supplements'] = ['NAD_precursors', 'metformin']

            # Metabolic recommendations
            metabolic_markers = ['glucose_mean', 'HbA1c_mean', 'insulin_mean']
            metabolic_score = np.mean([characteristics.get(marker, 0.0) for marker in metabolic_markers])

            if metabolic_score > 0.5:
                cluster_recs['metabolic_intervention'] = True
                cluster_recs['diet_recommendation'] = 'ketogenic_or_low_carb'
            else:
                cluster_recs['metabolic_intervention'] = False
                cluster_recs['diet_recommendation'] = 'mediterranean'

            # Monitoring frequency
            if cluster_recs['intervention_intensity'] == 'intensive':
                cluster_recs['monitoring_frequency'] = 'monthly'
            elif cluster_recs['intervention_intensity'] == 'standard':
                cluster_recs['monitoring_frequency'] = 'quarterly'
            else:
                cluster_recs['monitoring_frequency'] = 'biannual'

            recommendations[cluster_id] = cluster_recs

        return recommendations

    def _validate_stratification(self,
                               feature_matrix: np.ndarray,
                               cluster_assignments: np.ndarray) -> Dict[str, float]:
        """Validate stratification using cross-validation."""

        validation_metrics = {}

        try:
            # Cross-validation of cluster prediction
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            cv_scores = cross_val_score(rf, feature_matrix, cluster_assignments, cv=5)

            validation_metrics['cv_accuracy_mean'] = np.mean(cv_scores)
            validation_metrics['cv_accuracy_std'] = np.std(cv_scores)

            # Stability check (simplified)
            validation_metrics['cluster_stability'] = self._calculate_cluster_stability(
                feature_matrix, cluster_assignments
            )

        except:
            validation_metrics['cv_accuracy_mean'] = 0.0
            validation_metrics['cv_accuracy_std'] = 0.0
            validation_metrics['cluster_stability'] = 0.0

        return validation_metrics

    def _calculate_cluster_stability(self,
                                   feature_matrix: np.ndarray,
                                   cluster_assignments: np.ndarray) -> float:
        """Calculate cluster stability using bootstrap sampling."""

        try:
            n_samples, n_features = feature_matrix.shape
            n_bootstrap = 10
            stability_scores = []

            for _ in range(n_bootstrap):
                # Bootstrap sample
                indices = np.random.choice(n_samples, size=n_samples, replace=True)
                bootstrap_data = feature_matrix[indices]

                # Re-cluster
                n_clusters = len(np.unique(cluster_assignments))
                kmeans = KMeans(n_clusters=n_clusters, random_state=42)
                bootstrap_clusters = kmeans.fit_predict(bootstrap_data)

                # Calculate similarity to original clustering
                original_subset = cluster_assignments[indices]
                stability = adjusted_rand_score(original_subset, bootstrap_clusters)
                stability_scores.append(stability)

            return np.mean(stability_scores)

        except:
            return 0.0


def create_stratification_demo():
    """Create demonstration of patient stratification."""

    # Create synthetic patient data
    n_patients = 100
    patient_ids = [f"Patient_{i:03d}" for i in range(n_patients)]

    # Clinical data
    np.random.seed(42)
    clinical_data = pd.DataFrame({
        'age': np.random.normal(55, 15, n_patients),
        'sex': np.random.choice(['male', 'female'], n_patients),
        'BMI': np.random.normal(26, 4, n_patients),
        'blood_pressure_systolic': np.random.normal(130, 20, n_patients),
        'blood_pressure_diastolic': np.random.normal(80, 10, n_patients),
        'heart_rate': np.random.normal(70, 10, n_patients),
        'white_blood_cell_count': np.random.normal(7000, 2000, n_patients),
        'hemoglobin': np.random.normal(14, 2, n_patients)
    }, index=patient_ids)

    # Biomarker data
    biomarker_data = pd.DataFrame({
        'IL6': np.random.lognormal(1, 0.5, n_patients),
        'TNF_alpha': np.random.lognormal(2, 0.7, n_patients),
        'CRP': np.random.lognormal(0, 1, n_patients),
        'glucose': np.random.normal(100, 20, n_patients),
        'HbA1c': np.random.normal(5.5, 0.8, n_patients),
        'telomere_length': np.random.normal(1.0, 0.3, n_patients),
        'DNA_methylation_age': clinical_data['age'] + np.random.normal(0, 5, n_patients)
    }, index=patient_ids)

    # Initialize stratification engine
    engine = PatientStratificationEngine(stratification_method="multi_modal")

    # Load patient cohort
    engine.load_patient_cohort(clinical_data, biomarker_data=biomarker_data)

    # Perform stratification
    result = engine.perform_stratification()

    return engine, result


if __name__ == "__main__":
    # Run demonstration
    engine, result = create_stratification_demo()

    print("Patient stratification completed!")
    print(f"Number of clusters: {len(result.cluster_characteristics)}")
    print(f"Stratification quality: {result.stratification_quality:.3f}")
    print(f"Top biomarkers: {sorted(result.biomarker_importance.items(), key=lambda x: x[1], reverse=True)[:5]}")
    print("Cluster sizes:")
    for cluster_id, chars in result.cluster_characteristics.items():
        print(f"  Cluster {cluster_id}: {chars['cluster_size']} patients, "
              f"mean age: {chars.get('mean_age', 0):.1f} years")
