"""
Temporal Multi-Omics Modeling for CRISPR Toolkit Phase 3
========================================================

Advanced temporal modeling for tracking intervention responses
across multiple omics layers over time, with longitudinal
analysis and intervention trajectory prediction.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import savgol_filter

logger = logging.getLogger(__name__)


@dataclass
class TemporalOmicsData:
    """Container for temporal omics data across timepoints."""
    omics_type: str
    timepoints: List[str]
    expression_matrices: Dict[str, pd.DataFrame]  # timepoint -> expression matrix
    sample_metadata: pd.DataFrame
    intervention_timeline: Dict[str, str]  # timepoint -> intervention status


@dataclass
class InterventionTrajectory:
    """Results from intervention trajectory analysis."""
    patient_id: str
    omics_trajectories: Dict[str, pd.DataFrame]  # omics_type -> features x timepoints
    response_metrics: Dict[str, float]
    trajectory_classification: str  # "responder", "non_responder", "delayed_response"
    predicted_endpoints: Dict[str, float]
    confidence_scores: Dict[str, float]


@dataclass
class TemporalModelResult:
    """Results from temporal modeling analysis."""
    model_performance: Dict[str, float]
    feature_trajectories: Dict[str, pd.DataFrame]
    intervention_effects: Dict[str, Dict[str, float]]
    temporal_signatures: Dict[str, List[str]]
    prediction_intervals: Dict[str, Tuple[float, float]]


class TemporalMultiOmicsModel:
    """
    Advanced temporal modeling for multi-omics intervention tracking.

    This class implements sophisticated time-series analysis for tracking
    aging intervention responses across multiple omics layers, including
    trajectory prediction and intervention response classification.
    """

    def __init__(self, modeling_approach: str = "longitudinal_mixed"):
        """
        Initialize temporal multi-omics model.

        Args:
            modeling_approach: "longitudinal_mixed", "time_series", "trajectory_clustering"
        """
        self.modeling_approach = modeling_approach
        self.temporal_data = {}
        self.fitted_models = {}
        self.baseline_profiles = {}
        self.intervention_timepoints = []
        self.aging_trajectories = self._load_aging_trajectory_patterns()

        logger.info(f"Initialized TemporalMultiOmicsModel with {modeling_approach} approach")

    def _load_aging_trajectory_patterns(self) -> Dict[str, Dict[str, float]]:
        """Load known aging trajectory patterns from literature."""

        patterns = {
            'inflammatory_aging': {
                'IL6': 0.8,      # Strong upward trajectory with age
                'TNF_alpha': 0.7,
                'CRP': 0.6,
                'IL1_beta': 0.5,
                'NFkB_activity': 0.8
            },
            'cellular_senescence': {
                'p16_expression': 0.9,    # Strong senescence markers
                'p21_expression': 0.8,
                'SASP_IL6': 0.7,
                'SASP_IL8': 0.6,
                'telomere_length': -0.9   # Negative trajectory (shortening)
            },
            'metabolic_decline': {
                'glucose': 0.4,           # Gradual metabolic changes
                'insulin': 0.6,
                'HbA1c': 0.5,
                'adiponectin': -0.3,      # Beneficial hormone decline
                'HOMA_IR': 0.7
            },
            'mitochondrial_dysfunction': {
                'ATP_production': -0.8,    # Decline in mitochondrial function
                'COX_activity': -0.7,
                'citrate_synthase': -0.6,
                'ROS_production': 0.8,
                'mtDNA_copy_number': -0.5
            },
            'epigenetic_drift': {
                'DNA_methylation_age': 0.95,  # Strong epigenetic aging
                'H3K27me3_global': 0.6,
                'H3K4me3_promoters': -0.4,
                'SIRT1_activity': -0.7,
                'DNMT_activity': 0.5
            }
        }

        return patterns

    def add_temporal_omics_data(self, temporal_data: TemporalOmicsData):
        """Add temporal omics data for analysis."""
        self.temporal_data[temporal_data.omics_type] = temporal_data

        # Update intervention timepoints
        for timepoint in temporal_data.timepoints:
            if timepoint not in self.intervention_timepoints:
                self.intervention_timepoints.append(timepoint)

        logger.info(f"Added temporal {temporal_data.omics_type} data: "
                   f"{len(temporal_data.timepoints)} timepoints")

    def load_longitudinal_cohort(self,
                                cohort_data: Dict[str, Dict[str, pd.DataFrame]],
                                timepoint_metadata: pd.DataFrame):
        """
        Load longitudinal cohort data across multiple omics and timepoints.

        Args:
            cohort_data: {omics_type: {timepoint: expression_matrix}}
            timepoint_metadata: Metadata about timepoints and interventions
        """
        logger.info("Loading longitudinal cohort data...")

        for omics_type, timepoint_data in cohort_data.items():
            # Sort timepoints chronologically
            sorted_timepoints = sorted(timepoint_data.keys())

            # Create intervention timeline
            intervention_timeline = {}
            for timepoint in sorted_timepoints:
                if timepoint in timepoint_metadata.index:
                    intervention_timeline[timepoint] = timepoint_metadata.loc[timepoint, 'intervention_status']
                else:
                    intervention_timeline[timepoint] = 'unknown'

            # Create temporal omics data
            temporal_omics = TemporalOmicsData(
                omics_type=omics_type,
                timepoints=sorted_timepoints,
                expression_matrices=timepoint_data,
                sample_metadata=timepoint_metadata,
                intervention_timeline=intervention_timeline
            )

            self.add_temporal_omics_data(temporal_omics)

        # Sort all intervention timepoints
        self.intervention_timepoints = sorted(list(set(self.intervention_timepoints)))

        logger.info(f"Loaded {len(cohort_data)} omics types across "
                   f"{len(self.intervention_timepoints)} timepoints")

    def calculate_baseline_profiles(self):
        """Calculate baseline omics profiles before intervention."""

        logger.info("Calculating baseline omics profiles...")

        for omics_type, temporal_data in self.temporal_data.items():
            # Find baseline timepoint (usually first timepoint)
            baseline_timepoint = temporal_data.timepoints[0]

            if baseline_timepoint in temporal_data.expression_matrices:
                baseline_matrix = temporal_data.expression_matrices[baseline_timepoint]

                # Calculate baseline statistics
                baseline_profile = {
                    'mean_expression': baseline_matrix.mean(axis=1),
                    'std_expression': baseline_matrix.std(axis=1),
                    'cv_expression': baseline_matrix.std(axis=1) / baseline_matrix.mean(axis=1),
                    'n_samples': baseline_matrix.shape[1]
                }

                self.baseline_profiles[omics_type] = baseline_profile

                logger.info(f"Calculated baseline profile for {omics_type}: "
                           f"{len(baseline_profile['mean_expression'])} features")

    def fit_temporal_models(self,
                          features_to_model: Optional[Dict[str, List[str]]] = None) -> Dict[str, TemporalModelResult]:
        """
        Fit temporal models for intervention response prediction.

        Args:
            features_to_model: {omics_type: [feature_names]} to focus modeling

        Returns:
            Dictionary of temporal model results per omics type
        """
        logger.info("Fitting temporal models for intervention tracking...")

        if not self.baseline_profiles:
            self.calculate_baseline_profiles()

        model_results = {}

        for omics_type, temporal_data in self.temporal_data.items():
            logger.info(f"Modeling temporal patterns for {omics_type}...")

            # Select features to model
            if features_to_model and omics_type in features_to_model:
                selected_features = features_to_model[omics_type]
            else:
                # Use most variable features
                baseline_cv = self.baseline_profiles[omics_type]['cv_expression']
                selected_features = baseline_cv.nlargest(100).index.tolist()

            # Prepare temporal feature matrix
            temporal_matrix, time_vector = self._prepare_temporal_matrix(
                temporal_data, selected_features
            )

            # Fit models based on approach
            if self.modeling_approach == "longitudinal_mixed":
                model_result = self._fit_longitudinal_mixed_model(
                    temporal_matrix, time_vector, selected_features, omics_type
                )
            elif self.modeling_approach == "time_series":
                model_result = self._fit_time_series_model(
                    temporal_matrix, time_vector, selected_features, omics_type
                )
            else:  # trajectory_clustering
                model_result = self._fit_trajectory_clustering_model(
                    temporal_matrix, time_vector, selected_features, omics_type
                )

            model_results[omics_type] = model_result

            # Store fitted models
            self.fitted_models[omics_type] = model_result

        logger.info(f"Fitted temporal models for {len(model_results)} omics types")
        return model_results

    def _prepare_temporal_matrix(self,
                               temporal_data: TemporalOmicsData,
                               selected_features: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare temporal feature matrix for modeling."""

        # Find common samples across timepoints
        common_samples = None
        for timepoint in temporal_data.timepoints:
            if timepoint in temporal_data.expression_matrices:
                timepoint_samples = set(temporal_data.expression_matrices[timepoint].columns)
                if common_samples is None:
                    common_samples = timepoint_samples
                else:
                    common_samples = common_samples.intersection(timepoint_samples)

        common_samples = list(common_samples)

        if len(common_samples) == 0:
            raise ValueError(f"No common samples found across timepoints for {temporal_data.omics_type}")

        # Create time vector (days or numerical timepoints)
        time_vector = np.array([self._timepoint_to_numeric(tp) for tp in temporal_data.timepoints])

        # Build temporal matrix: [timepoints x features x samples]
        temporal_matrices = []

        for timepoint in temporal_data.timepoints:
            if timepoint in temporal_data.expression_matrices:
                timepoint_matrix = temporal_data.expression_matrices[timepoint]

                # Select features and common samples
                available_features = [f for f in selected_features if f in timepoint_matrix.index]
                timepoint_data = timepoint_matrix.loc[available_features, common_samples]

                temporal_matrices.append(timepoint_data.values)
            else:
                # Fill missing timepoint with NaN
                n_features = len(selected_features)
                n_samples = len(common_samples)
                temporal_matrices.append(np.full((n_features, n_samples), np.nan))

        # Stack into 3D array
        temporal_matrix = np.stack(temporal_matrices, axis=0)  # [timepoints x features x samples]

        return temporal_matrix, time_vector

    def _timepoint_to_numeric(self, timepoint: str) -> float:
        """Convert timepoint string to numeric value (days)."""

        # Handle common timepoint formats
        timepoint_lower = timepoint.lower()

        if 'baseline' in timepoint_lower or 'day_0' in timepoint_lower:
            return 0.0
        elif 'week' in timepoint_lower:
            # Extract week number
            week_num = float(timepoint_lower.split('week_')[1].split('_')[0])
            return week_num * 7.0
        elif 'day' in timepoint_lower:
            # Extract day number
            day_num = float(timepoint_lower.split('day_')[1].split('_')[0])
            return day_num
        elif 'month' in timepoint_lower:
            # Extract month number
            month_num = float(timepoint_lower.split('month_')[1].split('_')[0])
            return month_num * 30.0
        else:
            # Try to extract any number
            import re
            numbers = re.findall(r'\d+', timepoint)
            if numbers:
                return float(numbers[0])
            else:
                return 0.0

    def _fit_longitudinal_mixed_model(self,
                                    temporal_matrix: np.ndarray,
                                    time_vector: np.ndarray,
                                    feature_names: List[str],
                                    omics_type: str) -> TemporalModelResult:
        """Fit longitudinal mixed effects model."""

        n_timepoints, n_features, n_samples = temporal_matrix.shape

        model_performance = {}
        feature_trajectories = {}
        intervention_effects = {}
        temporal_signatures = {}
        prediction_intervals = {}

        # Analyze each feature separately
        for i, feature_name in enumerate(feature_names):
            feature_data = temporal_matrix[:, i, :]  # [timepoints x samples]

            # Skip features with too much missing data
            if np.isnan(feature_data).sum() > (n_timepoints * n_samples * 0.5):
                continue

            # Fit linear trend for each sample
            sample_slopes = []
            sample_r2_scores = []

            for j in range(n_samples):
                sample_trajectory = feature_data[:, j]

                # Remove NaN timepoints
                valid_mask = ~np.isnan(sample_trajectory)
                if np.sum(valid_mask) < 3:  # Need at least 3 timepoints
                    continue

                valid_times = time_vector[valid_mask]
                valid_values = sample_trajectory[valid_mask]

                # Fit linear regression
                try:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(valid_times, valid_values)
                    sample_slopes.append(slope)
                    sample_r2_scores.append(r_value**2)
                except:
                    continue

            if len(sample_slopes) > 0:
                # Calculate feature trajectory statistics
                feature_trajectories[feature_name] = pd.DataFrame({
                    'timepoint': time_vector,
                    'mean_expression': np.nanmean(feature_data, axis=1),
                    'std_expression': np.nanstd(feature_data, axis=1),
                    'n_samples': np.sum(~np.isnan(feature_data), axis=1)
                })

                # Calculate intervention effects
                baseline_mean = np.nanmean(feature_data[0, :])  # First timepoint
                endpoint_mean = np.nanmean(feature_data[-1, :])  # Last timepoint

                intervention_effects[feature_name] = {
                    'baseline_mean': baseline_mean,
                    'endpoint_mean': endpoint_mean,
                    'absolute_change': endpoint_mean - baseline_mean,
                    'relative_change': (endpoint_mean - baseline_mean) / baseline_mean if baseline_mean != 0 else 0,
                    'mean_slope': np.mean(sample_slopes),
                    'slope_std': np.std(sample_slopes),
                    'mean_r2': np.mean(sample_r2_scores)
                }

        # Identify temporal signatures (features with strong temporal patterns)
        strong_responders = []
        for feature_name, effects in intervention_effects.items():
            if abs(effects['mean_slope']) > 0.1 and effects['mean_r2'] > 0.5:
                strong_responders.append(feature_name)

        temporal_signatures['strong_temporal_response'] = strong_responders

        # Calculate overall model performance
        model_performance['n_features_modeled'] = len(intervention_effects)
        model_performance['mean_r2'] = np.mean([effects['mean_r2'] for effects in intervention_effects.values()])
        model_performance['temporal_signature_size'] = len(strong_responders)

        return TemporalModelResult(
            model_performance=model_performance,
            feature_trajectories=pd.DataFrame(feature_trajectories),
            intervention_effects=intervention_effects,
            temporal_signatures=temporal_signatures,
            prediction_intervals=prediction_intervals
        )

    def _fit_time_series_model(self,
                             temporal_matrix: np.ndarray,
                             time_vector: np.ndarray,
                             feature_names: List[str],
                             omics_type: str) -> TemporalModelResult:
        """Fit time series model with smoothing."""

        n_timepoints, n_features, n_samples = temporal_matrix.shape

        model_performance = {}
        feature_trajectories = {}
        intervention_effects = {}
        temporal_signatures = {}
        prediction_intervals = {}

        # Apply smoothing to trajectories
        for i, feature_name in enumerate(feature_names):
            feature_data = temporal_matrix[:, i, :]  # [timepoints x samples]

            # Calculate mean trajectory
            mean_trajectory = np.nanmean(feature_data, axis=1)

            # Apply smoothing if enough timepoints
            if len(time_vector) >= 5 and not np.any(np.isnan(mean_trajectory)):
                try:
                    # Savitzky-Golay smoothing
                    window_length = min(5, len(time_vector) if len(time_vector) % 2 == 1 else len(time_vector) - 1)
                    if window_length >= 3:
                        smoothed_trajectory = savgol_filter(mean_trajectory, window_length, 2)
                    else:
                        smoothed_trajectory = mean_trajectory
                except:
                    smoothed_trajectory = mean_trajectory
            else:
                smoothed_trajectory = mean_trajectory

            # Store trajectory
            feature_trajectories[feature_name] = pd.DataFrame({
                'timepoint': time_vector,
                'raw_mean': mean_trajectory,
                'smoothed_mean': smoothed_trajectory,
                'std_expression': np.nanstd(feature_data, axis=1)
            })

            # Calculate derivative (rate of change)
            if len(time_vector) > 1:
                time_diff = np.diff(time_vector)
                trajectory_diff = np.diff(smoothed_trajectory)
                rates = trajectory_diff / time_diff

                intervention_effects[feature_name] = {
                    'mean_rate': np.mean(rates),
                    'max_rate': np.max(rates),
                    'min_rate': np.min(rates),
                    'rate_variance': np.var(rates),
                    'total_change': smoothed_trajectory[-1] - smoothed_trajectory[0]
                }

        # Identify features with rapid changes (temporal signatures)
        rapid_changers = []
        for feature_name, effects in intervention_effects.items():
            if abs(effects['mean_rate']) > np.percentile([abs(e['mean_rate']) for e in intervention_effects.values()], 75):
                rapid_changers.append(feature_name)

        temporal_signatures['rapid_response'] = rapid_changers

        # Model performance
        model_performance['n_features_modeled'] = len(feature_trajectories)
        model_performance['rapid_responders'] = len(rapid_changers)

        return TemporalModelResult(
            model_performance=model_performance,
            feature_trajectories=pd.DataFrame(feature_trajectories),
            intervention_effects=intervention_effects,
            temporal_signatures=temporal_signatures,
            prediction_intervals=prediction_intervals
        )

    def _fit_trajectory_clustering_model(self,
                                       temporal_matrix: np.ndarray,
                                       time_vector: np.ndarray,
                                       feature_names: List[str],
                                       omics_type: str) -> TemporalModelResult:
        """Fit trajectory clustering model."""

        # This is a simplified implementation
        # In practice, would use more sophisticated trajectory clustering

        from sklearn.cluster import KMeans

        n_timepoints, n_features, n_samples = temporal_matrix.shape

        # Reshape for clustering: [samples x (features * timepoints)]
        reshaped_data = temporal_matrix.transpose(2, 1, 0).reshape(n_samples, -1)

        # Remove samples with too much missing data
        valid_samples = ~np.isnan(reshaped_data).any(axis=1)
        clean_data = reshaped_data[valid_samples]

        model_performance = {}
        feature_trajectories = {}
        intervention_effects = {}
        temporal_signatures = {}
        prediction_intervals = {}

        if len(clean_data) > 10:  # Need sufficient samples for clustering
            # Cluster trajectories
            n_clusters = min(5, len(clean_data) // 3)
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            cluster_labels = kmeans.fit_predict(clean_data)

            # Analyze clusters
            for cluster_id in range(n_clusters):
                cluster_mask = cluster_labels == cluster_id
                cluster_data = clean_data[cluster_mask]

                # Calculate cluster centroid trajectory
                cluster_centroid = cluster_data.mean(axis=0).reshape(n_features, n_timepoints)

                # Store cluster characteristics
                for i, feature_name in enumerate(feature_names):
                    feature_trajectory = cluster_centroid[i, :]

                    if feature_name not in feature_trajectories:
                        feature_trajectories[feature_name] = pd.DataFrame()

                    feature_trajectories[feature_name][f'cluster_{cluster_id}'] = feature_trajectory

            temporal_signatures['trajectory_clusters'] = list(range(n_clusters))
            model_performance['n_clusters'] = n_clusters
            model_performance['silhouette_score'] = 0.5  # Placeholder

        model_performance['n_features_modeled'] = len(feature_names)

        return TemporalModelResult(
            model_performance=model_performance,
            feature_trajectories=pd.DataFrame(feature_trajectories),
            intervention_effects=intervention_effects,
            temporal_signatures=temporal_signatures,
            prediction_intervals=prediction_intervals
        )

    def predict_intervention_trajectory(self,
                                      patient_baseline: Dict[str, pd.Series],
                                      prediction_timepoints: List[str]) -> InterventionTrajectory:
        """
        Predict intervention trajectory for a new patient.

        Args:
            patient_baseline: {omics_type: baseline_expression}
            prediction_timepoints: List of timepoints to predict

        Returns:
            InterventionTrajectory with predicted responses
        """
        logger.info("Predicting intervention trajectory for new patient...")

        if not self.fitted_models:
            raise ValueError("No fitted models available. Run fit_temporal_models() first.")

        patient_id = f"patient_{np.random.randint(1000, 9999)}"
        omics_trajectories = {}
        response_metrics = {}
        confidence_scores = {}

        # Convert prediction timepoints to numeric
        prediction_times = np.array([self._timepoint_to_numeric(tp) for tp in prediction_timepoints])

        for omics_type, baseline_data in patient_baseline.items():
            if omics_type in self.fitted_models:
                model_result = self.fitted_models[omics_type]

                # Predict trajectory for each feature
                feature_predictions = {}

                for feature_name in baseline_data.index:
                    if feature_name in model_result.intervention_effects:
                        baseline_value = baseline_data[feature_name]
                        effects = model_result.intervention_effects[feature_name]

                        # Simple linear prediction based on mean slope
                        predicted_values = baseline_value + effects['mean_slope'] * prediction_times
                        feature_predictions[feature_name] = predicted_values

                # Store trajectory predictions
                if feature_predictions:
                    trajectory_df = pd.DataFrame(feature_predictions, index=prediction_timepoints)
                    omics_trajectories[omics_type] = trajectory_df

                    # Calculate response metrics
                    baseline_mean = baseline_data.mean()
                    predicted_endpoint = trajectory_df.iloc[-1].mean()

                    response_metrics[f'{omics_type}_response'] = (predicted_endpoint - baseline_mean) / baseline_mean
                    confidence_scores[f'{omics_type}_confidence'] = 0.7  # Placeholder

        # Classify overall response
        overall_response = np.mean(list(response_metrics.values()))

        if overall_response > 0.1:
            trajectory_classification = "responder"
        elif overall_response < -0.1:
            trajectory_classification = "non_responder"
        else:
            trajectory_classification = "delayed_response"

        # Predict endpoints
        predicted_endpoints = {}
        for omics_type, trajectory in omics_trajectories.items():
            predicted_endpoints[f'{omics_type}_endpoint'] = trajectory.iloc[-1].mean()

        return InterventionTrajectory(
            patient_id=patient_id,
            omics_trajectories=omics_trajectories,
            response_metrics=response_metrics,
            trajectory_classification=trajectory_classification,
            predicted_endpoints=predicted_endpoints,
            confidence_scores=confidence_scores
        )

    def generate_temporal_report(self) -> Dict[str, any]:
        """Generate comprehensive temporal modeling report."""

        report = {
            'temporal_summary': {},
            'model_performance': {},
            'intervention_signatures': {},
            'trajectory_predictions': {},
            'recommendations': []
        }

        # Temporal data summary
        report['temporal_summary'] = {
            'omics_types': list(self.temporal_data.keys()),
            'timepoints': len(self.intervention_timepoints),
            'timespan_days': max([self._timepoint_to_numeric(tp) for tp in self.intervention_timepoints]),
            'modeling_approach': self.modeling_approach
        }

        # Model performance summary
        if self.fitted_models:
            performance_metrics = {}
            for omics_type, model_result in self.fitted_models.items():
                performance_metrics[omics_type] = model_result.model_performance

            report['model_performance'] = performance_metrics

            # Intervention signatures
            signatures = {}
            for omics_type, model_result in self.fitted_models.items():
                signatures[omics_type] = model_result.temporal_signatures

            report['intervention_signatures'] = signatures

        # Generate recommendations
        recommendations = self._generate_temporal_recommendations()
        report['recommendations'] = recommendations

        return report

    def _generate_temporal_recommendations(self) -> List[str]:
        """Generate recommendations based on temporal analysis."""

        recommendations = []

        if self.fitted_models:
            # Check model quality
            mean_performance = np.mean([
                model.model_performance.get('mean_r2', 0)
                for model in self.fitted_models.values()
            ])

            if mean_performance > 0.7:
                recommendations.append("High-quality temporal models - proceed with trajectory predictions")
            elif mean_performance > 0.4:
                recommendations.append("Moderate model quality - collect additional timepoints for better predictions")
            else:
                recommendations.append("Low model quality - consider alternative modeling approaches")

            # Check temporal signatures
            total_signatures = sum([
                len(model.temporal_signatures.get('strong_temporal_response', []))
                for model in self.fitted_models.values()
            ])

            if total_signatures > 50:
                recommendations.append("Strong temporal signatures detected - intervention showing clear effects")
            elif total_signatures > 10:
                recommendations.append("Moderate temporal response - monitor for sustained effects")
            else:
                recommendations.append("Limited temporal signatures - consider intervention optimization")

        # General recommendations
        recommendations.extend([
            "Continue longitudinal monitoring to improve model accuracy",
            "Consider adding more frequent early timepoints to capture initial response",
            "Validate temporal patterns with independent cohorts"
        ])

        return recommendations


def create_temporal_modeling_demo():
    """Create demonstration of temporal multi-omics modeling."""

    # Create synthetic longitudinal data
    n_patients = 50
    n_genes = 200
    timepoints = ['baseline', 'week_2', 'week_4', 'week_8', 'week_12']

    patient_ids = [f"Patient_{i:03d}" for i in range(n_patients)]
    gene_names = [f"Gene_{i}" for i in range(n_genes)]

    # Create cohort data with intervention effects
    np.random.seed(42)
    cohort_data = {}

    for omics_type in ['transcriptomics', 'proteomics']:
        timepoint_data = {}

        for i, timepoint in enumerate(timepoints):
            # Simulate intervention effects over time
            time_effect = i * 0.2  # Gradual improvement

            # Generate expression data with intervention effect
            expression_matrix = pd.DataFrame(
                np.random.lognormal(5 + time_effect * np.random.randn(n_genes, 1), 0.5, (n_genes, n_patients)),
                index=gene_names,
                columns=patient_ids
            )

            timepoint_data[timepoint] = expression_matrix

        cohort_data[omics_type] = timepoint_data

    # Create timepoint metadata
    timepoint_metadata = pd.DataFrame({
        'intervention_status': ['baseline', 'treatment', 'treatment', 'treatment', 'treatment'],
        'days_from_baseline': [0, 14, 28, 56, 84]
    }, index=timepoints)

    # Initialize temporal model
    temporal_model = TemporalMultiOmicsModel(modeling_approach="longitudinal_mixed")

    # Load longitudinal cohort
    temporal_model.load_longitudinal_cohort(cohort_data, timepoint_metadata)

    # Fit temporal models
    model_results = temporal_model.fit_temporal_models()

    # Generate report
    report = temporal_model.generate_temporal_report()

    # Test trajectory prediction for new patient
    new_patient_baseline = {
        'transcriptomics': pd.Series(np.random.lognormal(5, 0.5, n_genes), index=gene_names),
        'proteomics': pd.Series(np.random.lognormal(5, 0.5, n_genes), index=gene_names)
    }

    predicted_trajectory = temporal_model.predict_intervention_trajectory(
        new_patient_baseline, ['week_4', 'week_8', 'week_12']
    )

    return temporal_model, model_results, report, predicted_trajectory


if __name__ == "__main__":
    # Run demonstration
    temporal_model, model_results, report, trajectory = create_temporal_modeling_demo()

    print("Temporal multi-omics modeling completed!")
    print(f"Model results for {len(model_results)} omics types")
    print(f"Timespan: {report['temporal_summary']['timespan_days']} days")
    print(f"Predicted trajectory classification: {trajectory.trajectory_classification}")
    print(f"Response metrics: {trajectory.response_metrics}")
    print("Recommendations:")
    for rec in report['recommendations']:
        print(f"  - {rec}")
