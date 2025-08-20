"""
Multi-Omics Data Fusion Module for CRISPR Toolkit Phase 3
=========================================================

Advanced multi-omics integration algorithms for comprehensive
aging intervention analysis combining transcriptomics, proteomics,
metabolomics, and epigenomics data.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


@dataclass
class OmicsDataLayer:
    """Container for a single omics data layer."""
    data_type: str  # 'transcriptomics', 'proteomics', 'metabolomics', 'epigenomics'
    expression_matrix: pd.DataFrame
    feature_annotations: pd.DataFrame
    sample_metadata: pd.DataFrame
    quality_metrics: Dict[str, float]


@dataclass
class MultiOmicsIntegrationResult:
    """Results from multi-omics integration analysis."""
    integrated_features: pd.DataFrame
    dimension_reduction: Dict[str, np.ndarray]
    cluster_assignments: Dict[str, np.ndarray]
    pathway_scores: Dict[str, float]
    biomarker_signatures: Dict[str, List[str]]
    temporal_patterns: Optional[Dict[str, np.ndarray]]


class MultiOmicsFusion:
    """
    Advanced multi-omics data fusion for aging intervention research.

    This class implements multiple integration strategies including
    early fusion, late fusion, and intermediate fusion approaches
    for comprehensive multi-omics analysis.
    """

    def __init__(self, integration_strategy: str = "intermediate",
                 normalization_method: str = "standard"):
        """
        Initialize multi-omics fusion analyzer.

        Args:
            integration_strategy: "early", "intermediate", or "late"
            normalization_method: "standard", "robust", or "quantile"
        """
        self.integration_strategy = integration_strategy
        self.normalization_method = normalization_method
        self.omics_layers = {}
        self.integration_weights = {}
        self.scalers = {}

        logger.info(f"Initialized MultiOmicsFusion with {integration_strategy} strategy")

    def add_omics_layer(self, layer: OmicsDataLayer, weight: float = 1.0):
        """Add an omics data layer to the fusion analysis."""
        self.omics_layers[layer.data_type] = layer
        self.integration_weights[layer.data_type] = weight

        logger.info(f"Added {layer.data_type} layer: {layer.expression_matrix.shape}")

    def preprocess_omics_data(self) -> Dict[str, pd.DataFrame]:
        """
        Preprocess and normalize all omics layers.

        Returns:
            Dictionary of preprocessed data matrices
        """
        logger.info("Preprocessing omics data layers...")

        preprocessed_data = {}

        for data_type, layer in self.omics_layers.items():
            # Handle missing values
            data = layer.expression_matrix.copy()
            data = data.fillna(data.median())

            # Log transform if appropriate (for expression data)
            if data_type in ['transcriptomics', 'proteomics']:
                # Add small constant to avoid log(0)
                data = np.log2(data + 1)

            # Normalization
            if self.normalization_method == "standard":
                scaler = StandardScaler()
            elif self.normalization_method == "robust":
                from sklearn.preprocessing import RobustScaler
                scaler = RobustScaler()
            else:  # quantile
                from sklearn.preprocessing import QuantileTransformer
                scaler = QuantileTransformer(output_distribution='normal')

            # Fit scaler and transform data
            data_normalized = pd.DataFrame(
                scaler.fit_transform(data.T).T,
                index=data.index,
                columns=data.columns
            )

            self.scalers[data_type] = scaler
            preprocessed_data[data_type] = data_normalized

            logger.info(f"Preprocessed {data_type}: {data_normalized.shape}")

        return preprocessed_data

    def perform_early_fusion(self, preprocessed_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """
        Perform early fusion by concatenating features from all omics layers.

        Args:
            preprocessed_data: Dictionary of preprocessed omics data

        Returns:
            Concatenated feature matrix
        """
        logger.info("Performing early fusion integration...")

        # Find common samples
        common_samples = None
        for data_type, data in preprocessed_data.items():
            if common_samples is None:
                common_samples = set(data.columns)
            else:
                common_samples = common_samples.intersection(set(data.columns))

        common_samples = list(common_samples)

        if len(common_samples) == 0:
            raise ValueError("No common samples found across omics layers")

        # Concatenate features
        concatenated_features = []
        feature_names = []

        for data_type, data in preprocessed_data.items():
            # Apply integration weight
            weight = self.integration_weights[data_type]
            weighted_data = data[common_samples] * weight

            concatenated_features.append(weighted_data)

            # Create feature names with data type prefix
            prefixed_names = [f"{data_type}_{feature}" for feature in data.index]
            feature_names.extend(prefixed_names)

        # Combine all features
        integrated_matrix = pd.concat(concatenated_features, axis=0)
        integrated_matrix.index = feature_names

        logger.info(f"Early fusion result: {integrated_matrix.shape}")
        return integrated_matrix

    def perform_intermediate_fusion(self, preprocessed_data: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """
        Perform intermediate fusion using dimensionality reduction on each layer.

        Args:
            preprocessed_data: Dictionary of preprocessed omics data

        Returns:
            Integrated feature matrix from reduced dimensions
        """
        logger.info("Performing intermediate fusion integration...")

        # Find common samples
        common_samples = self._find_common_samples(preprocessed_data)

        # Perform PCA on each omics layer
        reduced_features = []
        feature_names = []

        for data_type, data in preprocessed_data.items():
            layer_data = data[common_samples].T  # samples x features

            # Determine number of components (explained variance > 80%)
            pca = PCA()
            pca.fit(layer_data)
            cumsum_var = np.cumsum(pca.explained_variance_ratio_)
            n_components = np.argmax(cumsum_var >= 0.8) + 1
            n_components = min(n_components, min(layer_data.shape) - 1)

            # Apply PCA with selected components
            pca_final = PCA(n_components=n_components)
            reduced_data = pca_final.fit_transform(layer_data)

            # Apply integration weight
            weight = self.integration_weights[data_type]
            weighted_reduced = reduced_data * weight

            reduced_features.append(weighted_reduced.T)  # features x samples

            # Create feature names
            pc_names = [f"{data_type}_PC{i+1}" for i in range(n_components)]
            feature_names.extend(pc_names)

            logger.info(f"{data_type}: reduced to {n_components} components "
                       f"(explaining {cumsum_var[n_components-1]:.2%} variance)")

        # Combine reduced features
        integrated_matrix = np.vstack(reduced_features)
        integrated_df = pd.DataFrame(
            integrated_matrix,
            index=feature_names,
            columns=common_samples
        )

        logger.info(f"Intermediate fusion result: {integrated_df.shape}")
        return integrated_df

    def perform_late_fusion(self, preprocessed_data: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
        """
        Perform late fusion by analyzing each omics layer separately.

        Args:
            preprocessed_data: Dictionary of preprocessed omics data

        Returns:
            Dictionary of analysis results for each omics layer
        """
        logger.info("Performing late fusion integration...")

        layer_results = {}

        for data_type, data in preprocessed_data.items():
            # Perform layer-specific analysis
            layer_analysis = self._analyze_single_layer(data, data_type)
            layer_results[data_type] = layer_analysis

        return layer_results

    def integrate_omics_data(self) -> MultiOmicsIntegrationResult:
        """
        Integrate all omics layers using the specified strategy.

        Returns:
            MultiOmicsIntegrationResult with integrated analysis
        """
        logger.info(f"Starting multi-omics integration with {self.integration_strategy} strategy...")

        # Preprocess data
        preprocessed_data = self.preprocess_omics_data()

        # Perform integration based on strategy
        if self.integration_strategy == "early":
            integrated_features = self.perform_early_fusion(preprocessed_data)
        elif self.integration_strategy == "intermediate":
            integrated_features = self.perform_intermediate_fusion(preprocessed_data)
        else:  # late fusion
            layer_results = self.perform_late_fusion(preprocessed_data)
            # For late fusion, combine results at the end
            integrated_features = self._combine_late_fusion_results(layer_results)

        # Perform downstream analysis
        dimension_reduction = self._perform_dimension_reduction(integrated_features)
        cluster_assignments = self._perform_clustering(integrated_features)
        pathway_scores = self._calculate_pathway_scores(integrated_features)
        biomarker_signatures = self._identify_biomarker_signatures(integrated_features)

        # Temporal analysis if timepoint data available
        temporal_patterns = self._analyze_temporal_patterns(integrated_features)

        result = MultiOmicsIntegrationResult(
            integrated_features=integrated_features,
            dimension_reduction=dimension_reduction,
            cluster_assignments=cluster_assignments,
            pathway_scores=pathway_scores,
            biomarker_signatures=biomarker_signatures,
            temporal_patterns=temporal_patterns
        )

        logger.info("Multi-omics integration completed successfully")
        return result

    def _find_common_samples(self, data_dict: Dict[str, pd.DataFrame]) -> List[str]:
        """Find samples common to all omics layers."""
        common_samples = None
        for data in data_dict.values():
            if common_samples is None:
                common_samples = set(data.columns)
            else:
                common_samples = common_samples.intersection(set(data.columns))
        return list(common_samples)

    def _analyze_single_layer(self, data: pd.DataFrame, data_type: str) -> pd.DataFrame:
        """Analyze a single omics layer."""
        # Perform PCA
        pca = PCA(n_components=min(10, min(data.shape) - 1))
        pca_result = pca.fit_transform(data.T)

        # Create results dataframe
        pc_names = [f"{data_type}_PC{i+1}" for i in range(pca_result.shape[1])]
        results_df = pd.DataFrame(
            pca_result.T,
            index=pc_names,
            columns=data.columns
        )

        return results_df

    def _combine_late_fusion_results(self, layer_results: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Combine results from late fusion analysis."""
        # Simply concatenate the top components from each layer
        combined_features = []

        for data_type, results in layer_results.items():
            # Take top 5 components from each layer
            top_components = results.iloc[:5]
            combined_features.append(top_components)

        return pd.concat(combined_features, axis=0)

    def _perform_dimension_reduction(self, data: pd.DataFrame) -> Dict[str, np.ndarray]:
        """Perform multiple dimension reduction techniques."""
        data_matrix = data.T  # samples x features

        results = {}

        # PCA
        try:
            pca = PCA(n_components=min(10, min(data_matrix.shape) - 1))
            results['pca'] = pca.fit_transform(data_matrix)
        except:
            results['pca'] = np.empty((data_matrix.shape[0], 0))

        # t-SNE (if feasible)
        if data_matrix.shape[0] > 10 and data_matrix.shape[1] > 1:
            try:
                tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, data_matrix.shape[0]//4))
                results['tsne'] = tsne.fit_transform(data_matrix)
            except:
                results['tsne'] = np.empty((data_matrix.shape[0], 2))
        else:
            results['tsne'] = np.empty((data_matrix.shape[0], 2))

        return results

    def _perform_clustering(self, data: pd.DataFrame) -> Dict[str, np.ndarray]:
        """Perform clustering analysis."""
        data_matrix = data.T  # samples x features

        results = {}

        # K-means clustering with different k values
        for k in [2, 3, 4, 5]:
            if k < data_matrix.shape[0]:
                try:
                    kmeans = KMeans(n_clusters=k, random_state=42)
                    cluster_labels = kmeans.fit_predict(data_matrix)
                    results[f'kmeans_k{k}'] = cluster_labels
                except:
                    results[f'kmeans_k{k}'] = np.zeros(data_matrix.shape[0])

        return results

    def _calculate_pathway_scores(self, data: pd.DataFrame) -> Dict[str, float]:
        """Calculate pathway-level scores from integrated data."""
        # This is a simplified version - real implementation would use pathway databases
        pathway_scores = {}

        # Group features by omics type and calculate summary scores
        omics_types = ['transcriptomics', 'proteomics', 'metabolomics', 'epigenomics']

        for omics_type in omics_types:
            # Find features from this omics type
            omics_features = [idx for idx in data.index if idx.startswith(omics_type)]

            if omics_features:
                omics_data = data.loc[omics_features]
                # Calculate variance across samples as a simple pathway activity score
                pathway_scores[f'{omics_type}_activity'] = omics_data.var(axis=1).mean()
            else:
                pathway_scores[f'{omics_type}_activity'] = 0.0

        return pathway_scores

    def _identify_biomarker_signatures(self, data: pd.DataFrame) -> Dict[str, List[str]]:
        """Identify biomarker signatures from integrated data."""
        signatures = {}

        # Find high-variance features as potential biomarkers
        feature_variance = data.var(axis=1)
        top_features = feature_variance.nlargest(20).index.tolist()

        signatures['high_variance'] = top_features

        # Group by omics type
        for omics_type in ['transcriptomics', 'proteomics', 'metabolomics', 'epigenomics']:
            omics_features = [f for f in top_features if f.startswith(omics_type)]
            if omics_features:
                signatures[f'{omics_type}_biomarkers'] = omics_features

        return signatures

    def _analyze_temporal_patterns(self, data: pd.DataFrame) -> Optional[Dict[str, np.ndarray]]:
        """Analyze temporal patterns if timepoint data is available."""
        # Check if any layer has timepoint information
        has_temporal_data = False
        for layer in self.omics_layers.values():
            if 'timepoint' in layer.sample_metadata.columns:
                has_temporal_data = True
                break

        if not has_temporal_data:
            return None

        # Simplified temporal analysis
        temporal_patterns = {}

        # This would be more sophisticated in real implementation
        # For now, just return feature means over time
        temporal_patterns['feature_trends'] = data.mean(axis=1).values

        return temporal_patterns


class TemporalOmicsModel:
    """
    Model for temporal multi-omics analysis of aging interventions.
    """

    def __init__(self, time_resolution: str = "weeks"):
        """
        Initialize temporal omics model.

        Args:
            time_resolution: "hours", "days", "weeks", or "months"
        """
        self.time_resolution = time_resolution
        self.temporal_data = {}
        self.intervention_timepoints = []

    def add_temporal_data(self,
                         timepoint: str,
                         omics_data: Dict[str, pd.DataFrame],
                         intervention_status: str = "baseline"):
        """Add omics data for a specific timepoint."""
        self.temporal_data[timepoint] = {
            'omics_data': omics_data,
            'intervention_status': intervention_status
        }

        if timepoint not in self.intervention_timepoints:
            self.intervention_timepoints.append(timepoint)

        logger.info(f"Added temporal data for timepoint: {timepoint}")

    def analyze_intervention_trajectory(self) -> Dict[str, pd.DataFrame]:
        """Analyze the trajectory of omics changes over time."""

        # Sort timepoints
        sorted_timepoints = sorted(self.intervention_timepoints)

        trajectories = {}

        # For each omics type, track changes over time
        for omics_type in ['transcriptomics', 'proteomics', 'metabolomics']:
            feature_trajectories = []

            for timepoint in sorted_timepoints:
                if (timepoint in self.temporal_data and
                    omics_type in self.temporal_data[timepoint]['omics_data']):

                    data = self.temporal_data[timepoint]['omics_data'][omics_type]
                    # Calculate mean expression across samples
                    mean_expression = data.mean(axis=1)
                    feature_trajectories.append(mean_expression)

            if feature_trajectories:
                trajectory_df = pd.concat(feature_trajectories, axis=1)
                trajectory_df.columns = sorted_timepoints
                trajectories[omics_type] = trajectory_df

        return trajectories

    def identify_temporal_biomarkers(self,
                                   trajectories: Dict[str, pd.DataFrame],
                                   min_fold_change: float = 1.5) -> Dict[str, List[str]]:
        """Identify features with significant temporal changes."""

        temporal_biomarkers = {}

        for omics_type, trajectory_data in trajectories.items():
            biomarkers = []

            for feature in trajectory_data.index:
                feature_trajectory = trajectory_data.loc[feature]

                # Calculate maximum fold change over time
                max_val = feature_trajectory.max()
                min_val = feature_trajectory.min()

                if min_val > 0:
                    fold_change = max_val / min_val

                    if fold_change >= min_fold_change:
                        biomarkers.append(feature)

            temporal_biomarkers[omics_type] = biomarkers

        return temporal_biomarkers


def create_multi_omics_integration_demo():
    """Create a demonstration of multi-omics integration."""

    # Create synthetic multi-omics data
    n_samples = 30
    sample_names = [f"Sample_{i}" for i in range(n_samples)]

    # Transcriptomics data
    transcriptomics_data = pd.DataFrame(
        np.random.lognormal(5, 1, (1000, n_samples)),
        index=[f"Gene_{i}" for i in range(1000)],
        columns=sample_names
    )

    # Proteomics data
    proteomics_data = pd.DataFrame(
        np.random.lognormal(10, 0.8, (200, n_samples)),
        index=[f"Protein_{i}" for i in range(200)],
        columns=sample_names
    )

    # Metabolomics data
    metabolomics_data = pd.DataFrame(
        np.random.lognormal(8, 1.2, (100, n_samples)),
        index=[f"Metabolite_{i}" for i in range(100)],
        columns=sample_names
    )

    # Sample metadata
    metadata = pd.DataFrame({
        'condition': ['control'] * 15 + ['treatment'] * 15,
        'age': np.random.randint(20, 80, n_samples),
        'timepoint': ['baseline'] * 10 + ['week_2'] * 10 + ['week_4'] * 10
    }, index=sample_names)

    # Create omics layers
    transcriptomics_layer = OmicsDataLayer(
        data_type='transcriptomics',
        expression_matrix=transcriptomics_data,
        feature_annotations=pd.DataFrame(),
        sample_metadata=metadata,
        quality_metrics={'missing_data': 0.02}
    )

    proteomics_layer = OmicsDataLayer(
        data_type='proteomics',
        expression_matrix=proteomics_data,
        feature_annotations=pd.DataFrame(),
        sample_metadata=metadata,
        quality_metrics={'missing_data': 0.05}
    )

    metabolomics_layer = OmicsDataLayer(
        data_type='metabolomics',
        expression_matrix=metabolomics_data,
        feature_annotations=pd.DataFrame(),
        sample_metadata=metadata,
        quality_metrics={'missing_data': 0.03}
    )

    # Initialize fusion analyzer
    fusion_analyzer = MultiOmicsFusion(integration_strategy="intermediate")

    # Add omics layers
    fusion_analyzer.add_omics_layer(transcriptomics_layer, weight=1.0)
    fusion_analyzer.add_omics_layer(proteomics_layer, weight=1.2)
    fusion_analyzer.add_omics_layer(metabolomics_layer, weight=0.8)

    # Perform integration
    integration_result = fusion_analyzer.integrate_omics_data()

    return fusion_analyzer, integration_result


if __name__ == "__main__":
    # Run demonstration
    fusion_analyzer, result = create_multi_omics_integration_demo()

    print("Multi-omics integration completed!")
    print(f"Integrated features: {result.integrated_features.shape}")
    print(f"Dimension reduction methods: {list(result.dimension_reduction.keys())}")
    print(f"Clustering results: {list(result.cluster_assignments.keys())}")
    print(f"Pathway scores: {len(result.pathway_scores)}")
    print(f"Biomarker signatures: {len(result.biomarker_signatures)}")
