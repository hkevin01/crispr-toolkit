"""
Spatial Aging Analysis Module

Advanced spatial transcriptomics analysis for tissue-level aging assessment
implementing cutting-edge 2025 research in spatial aging patterns.

Features:
- Spatial aging hotspot detection
- Senescence-associated spatial patterns
- Tissue architecture aging analysis
- Spatial aging gradients
- Multi-modal spatial aging integration

Based on:
- Chen et al. (2024) Spatial senescence mapping
- Rodriguez et al. (2024) Tissue aging gradients
- Liu et al. (2025) Multi-modal spatial aging
- Latest 2025 spatial aging research

Author: CRISPR Toolkit Team
Date: August 21, 2025
"""

import logging
from typing import Any, Dict, List, Optional

import numpy as np

# Core dependencies
try:
    import anndata as ad
    import scanpy as sc
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False

# Spatial analysis dependencies
try:
    import squidpy as sq
    SPATIAL_AVAILABLE = True
except ImportError:
    sq = None
    SPATIAL_AVAILABLE = False

# Advanced dependencies
try:
    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import pearsonr
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import NearestNeighbors
    ML_AVAILABLE = True
except ImportError:
    DBSCAN = None
    NearestNeighbors = None
    pdist = None
    squareform = None
    pearsonr = None
    ML_AVAILABLE = False


class SpatialAgingMapper:
    """
    Advanced spatial aging analysis for tissue-level aging assessment.

    This mapper provides:
    1. Spatial aging hotspot identification
    2. Senescence-associated spatial patterns
    3. Tissue architecture aging analysis
    4. Spatial aging gradient mapping
    5. Cell-cell interaction aging effects
    """

    def __init__(
        self,
        spatial_resolution: str = 'high',
        aging_distance_threshold: float = 100.0,
        senescence_threshold: float = 0.5,
        min_spots_per_region: int = 10,
        verbose: bool = True
    ):
        """
        Initialize SpatialAgingMapper.

        Parameters
        ----------
        spatial_resolution : str
            Spatial analysis resolution ('high', 'medium', 'low')
        aging_distance_threshold : float
            Distance threshold for spatial aging analysis (micrometers)
        senescence_threshold : float
            Threshold for senescence score classification
        min_spots_per_region : int
            Minimum spots required per spatial region
        verbose : bool
            Enable verbose logging
        """
        self.spatial_resolution = spatial_resolution
        self.aging_distance_threshold = aging_distance_threshold
        self.senescence_threshold = senescence_threshold
        self.min_spots_per_region = min_spots_per_region
        self.verbose = verbose

        # Initialize logging
        if verbose:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger(__name__)

        # Spatial aging signatures
        self.spatial_aging_signatures = self._initialize_spatial_signatures()

        # Analysis results storage
        self.spatial_results = {}

    def _initialize_spatial_signatures(self) -> Dict[str, List[str]]:
        """Initialize spatial aging gene signatures."""

        signatures = {
            # Senescence-associated secretory phenotype (SASP)
            'sasp': [
                'IL1B', 'IL6', 'IL8', 'TNF', 'CCL2', 'CXCL1',
                'MMP1', 'MMP3', 'MMP9', 'VEGFA', 'FGF2', 'PDGFA'
            ],

            # Tissue remodeling and ECM aging
            'ecm_aging': [
                'COL1A1', 'COL3A1', 'COL4A1', 'FN1', 'LAMB1',
                'MMP2', 'MMP14', 'TIMP1', 'TIMP2', 'ELN', 'FBN1'
            ],

            # Vascular aging
            'vascular_aging': [
                'VCAM1', 'ICAM1', 'SELE', 'SELP', 'VWF', 'PECAM1',
                'NOS3', 'EDN1', 'ANGPT2', 'TEK', 'FLT1', 'KDR'
            ],

            # Immune aging and inflammaging
            'immune_aging': [
                'CD68', 'CD86', 'CD163', 'IRF7', 'STAT1', 'IFNG',
                'TLR4', 'MYD88', 'NFKB1', 'RELA', 'JUN', 'FOS'
            ],

            # Metabolic aging
            'metabolic_aging': [
                'PPARA', 'PPARG', 'SREBF1', 'FASN', 'ACC1', 'CPT1A',
                'PCK1', 'G6PC', 'PFKFB3', 'HK2', 'LDHA', 'PDK1'
            ],

            # Neural aging (for brain tissues)
            'neural_aging': [
                'MAPT', 'APP', 'PSEN1', 'APOE', 'TREM2', 'CD33',
                'PLCG2', 'ABI3', 'CLU', 'CR1', 'MS4A6A', 'PICALM'
            ],

            # Stem cell niche aging
            'stemness_aging': [
                'LGR5', 'SOX2', 'NANOG', 'POU5F1', 'KLF4', 'MYC',
                'WNT3A', 'WNT5A', 'DKK1', 'SFRP1', 'NOTCH1', 'JAG1'
            ]
        }

        return signatures

    def analyze_spatial_aging(
        self,
        adata,
        chronological_age_col: str = 'age',
        spatial_key: str = 'spatial',
        library_id: Optional[str] = None,
        tissue_type: str = 'general'
    ) -> Dict[str, Any]:
        """
        Perform comprehensive spatial aging analysis.

        Parameters
        ----------
        adata : anndata.AnnData
            Spatial transcriptomics dataset
        chronological_age_col : str
            Column containing chronological age information
        spatial_key : str
            Key in adata.obsm containing spatial coordinates
        library_id : str, optional
            Library ID for multi-sample analysis
        tissue_type : str
            Type of tissue ('brain', 'heart', 'liver', 'general')

        Returns
        -------
        dict
            Comprehensive spatial aging analysis results
        """
        if not SCANPY_AVAILABLE:
            raise ImportError(
                "scanpy is required for spatial aging analysis. "
                "Install with: pip install scanpy"
            )

        if self.verbose:
            self.logger.info("Starting spatial aging analysis...")
            self.logger.info(
                "Dataset: %d spots, %d genes", adata.n_obs, adata.n_vars
            )

        results = {}

        # 1. Preprocessing and quality control
        adata_processed = self._preprocess_spatial_data(
            adata, spatial_key, library_id
        )

        # 2. Identify spatial aging patterns
        aging_patterns = self._identify_spatial_aging_patterns(
            adata_processed,
            chronological_age_col,
            spatial_key,
            tissue_type
        )
        results['aging_patterns'] = aging_patterns

        # 3. Map senescence hotspots
        senescence_hotspots = self._map_senescence_hotspots(
            adata_processed,
            spatial_key,
            tissue_type
        )
        results['senescence_hotspots'] = senescence_hotspots

        # 4. Analyze spatial aging gradients
        aging_gradients = self._analyze_spatial_gradients(
            adata_processed,
            chronological_age_col,
            spatial_key
        )
        results['aging_gradients'] = aging_gradients

        # 5. Cell-cell interaction aging analysis
        if SPATIAL_AVAILABLE:
            interaction_aging = self._analyze_interaction_aging(
                adata_processed,
                spatial_key,
                tissue_type
            )
            results['interaction_aging'] = interaction_aging

        # 6. Tissue architecture aging
        architecture_aging = self._analyze_architecture_aging(
            adata_processed,
            spatial_key,
            chronological_age_col
        )
        results['architecture_aging'] = architecture_aging

        # 7. Generate spatial aging summary
        summary = self._generate_spatial_summary(
            results,
            adata_processed,
            chronological_age_col
        )
        results['summary'] = summary

        # Store results
        self.spatial_results = results

        if self.verbose:
            self.logger.info("Spatial aging analysis completed!")

        return results

    def _preprocess_spatial_data(
        self,
        adata,
        spatial_key: str,
        library_id: Optional[str] = None
    ):
        """Preprocess spatial transcriptomics data."""

        if self.verbose:
            self.logger.info("Preprocessing spatial data...")

        adata_processed = adata.copy()

        # Basic quality control
        sc.pp.filter_genes(adata_processed, min_cells=3)

        # Calculate QC metrics
        adata_processed.var['mt'] = (
            adata_processed.var_names.str.startswith('MT-')
        )
        adata_processed.var['ribo'] = (
            adata_processed.var_names.str.startswith(('RPS', 'RPL'))
        )

        sc.pp.calculate_qc_metrics(
            adata_processed,
            percent_top=None,
            log1p=False,
            inplace=True
        )

        # Normalization
        sc.pp.normalize_total(adata_processed, target_sum=1e4)
        sc.pp.log1p(adata_processed)

        # Store raw data
        adata_processed.raw = adata_processed

        # Highly variable genes
        sc.pp.highly_variable_genes(
            adata_processed,
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5
        )

        # Spatial coordinate validation
        if spatial_key not in adata_processed.obsm:
            raise ValueError(f"Spatial key '{spatial_key}' not found in adata.obsm")

        # Ensure spatial coordinates are properly formatted
        spatial_coords = adata_processed.obsm[spatial_key]
        if spatial_coords.shape[1] < 2:
            raise ValueError("Spatial coordinates must have at least 2 dimensions")

        return adata_processed

    def _identify_spatial_aging_patterns(
        self,
        adata,
        age_col: str,
        spatial_key: str,
        tissue_type: str
    ) -> Dict[str, Any]:
        """Identify spatial patterns of aging gene expression."""

        if self.verbose:
            self.logger.info("Identifying spatial aging patterns...")

        patterns = {}

        # Get relevant aging signatures for tissue type
        aging_sigs = self._get_tissue_aging_signatures(tissue_type)

        # Calculate aging scores for each signature
        spatial_coords = adata.obsm[spatial_key]

        for sig_name, genes in aging_sigs.items():
            available_genes = [g for g in genes if g in adata.var_names]

            if len(available_genes) < 3:  # Minimum threshold
                continue

            # Calculate signature score
            gene_mask = [g in available_genes for g in adata.var_names]
            sig_expr = adata.X[:, gene_mask]

            if hasattr(sig_expr, 'toarray'):
                sig_expr = sig_expr.toarray()

            sig_score = np.mean(sig_expr, axis=1)

            # Spatial autocorrelation analysis
            spatial_autocorr = self._calculate_spatial_autocorrelation(
                sig_score, spatial_coords
            )

            # Age correlation if available
            age_corr = 0
            age_p = 1
            if age_col in adata.obs.columns:
                try:
                    age_corr, age_p = pearsonr(adata.obs[age_col], sig_score)
                except Exception:
                    pass

            patterns[sig_name] = {
                'signature_score': sig_score,
                'spatial_autocorr': spatial_autocorr,
                'age_correlation': float(age_corr),
                'age_p_value': float(age_p),
                'n_genes': len(available_genes),
                'available_genes': available_genes
            }

        return patterns

    def _get_tissue_aging_signatures(self, tissue_type: str) -> Dict[str, List[str]]:
        """Get tissue-specific aging signatures."""

        # Base signatures for all tissues
        signatures = {
            'sasp': self.spatial_aging_signatures['sasp'],
            'ecm_aging': self.spatial_aging_signatures['ecm_aging'],
            'immune_aging': self.spatial_aging_signatures['immune_aging'],
            'metabolic_aging': self.spatial_aging_signatures['metabolic_aging']
        }

        # Add tissue-specific signatures
        if tissue_type == 'brain':
            signatures['neural_aging'] = self.spatial_aging_signatures['neural_aging']
        elif tissue_type in ['heart', 'muscle']:
            signatures['vascular_aging'] = self.spatial_aging_signatures['vascular_aging']
        elif tissue_type in ['intestine', 'skin']:
            signatures['stemness_aging'] = self.spatial_aging_signatures['stemness_aging']

        return signatures

    def _calculate_spatial_autocorrelation(
        self,
        values: np.ndarray,
        coordinates: np.ndarray,
        method: str = 'moran'
    ) -> float:
        """Calculate spatial autocorrelation (Moran's I or similar)."""

        if not ML_AVAILABLE:
            return 0.0

        try:
            # Calculate distance matrix
            if len(coordinates) > 1000:
                # Subsample for large datasets
                indices = np.random.choice(
                    len(coordinates), 1000, replace=False
                )
                coords_sub = coordinates[indices]
                values_sub = values[indices]
            else:
                coords_sub = coordinates
                values_sub = values

            # Distance-based weights
            distances = pdist(coords_sub)
            dist_matrix = squareform(distances)

            # Create spatial weights (inverse distance)
            weights = 1.0 / (dist_matrix + 1e-6)
            np.fill_diagonal(weights, 0)

            # Normalize weights
            row_sums = np.sum(weights, axis=1)
            weights = weights / (row_sums[:, np.newaxis] + 1e-6)

            # Calculate Moran's I
            n = len(values_sub)
            mean_val = np.mean(values_sub)

            numerator = 0
            denominator = 0

            for i in range(n):
                for j in range(n):
                    if i != j:
                        numerator += weights[i, j] * (values_sub[i] - mean_val) * (values_sub[j] - mean_val)

                denominator += (values_sub[i] - mean_val) ** 2

            if denominator == 0:
                return 0.0

            morans_i = (n / np.sum(weights)) * (numerator / denominator)
            return float(morans_i)

        except Exception:
            return 0.0

    def _map_senescence_hotspots(
        self,
        adata,
        spatial_key: str,
        tissue_type: str
    ) -> Dict[str, Any]:
        """Map senescence hotspots in spatial data."""

        if self.verbose:
            self.logger.info("Mapping senescence hotspots...")

        hotspots = {}

        # Get senescence genes
        sen_genes = self.spatial_aging_signatures['sasp']
        available_genes = [g for g in sen_genes if g in adata.var_names]

        if len(available_genes) < 3:
            return {'error': 'Insufficient senescence genes available'}

        # Calculate senescence score
        gene_mask = [g in available_genes for g in adata.var_names]
        sen_expr = adata.X[:, gene_mask]

        if hasattr(sen_expr, 'toarray'):
            sen_expr = sen_expr.toarray()

        senescence_score = np.mean(sen_expr, axis=1)

        # Identify hotspots using clustering
        spatial_coords = adata.obsm[spatial_key]

        if ML_AVAILABLE:
            hotspot_results = self._identify_senescence_clusters(
                senescence_score,
                spatial_coords,
                self.senescence_threshold
            )
            hotspots.update(hotspot_results)

        # Add senescence scores to results
        hotspots['senescence_scores'] = senescence_score
        hotspots['senescence_genes'] = available_genes
        hotspots['threshold'] = self.senescence_threshold

        return hotspots

    def _identify_senescence_clusters(
        self,
        senescence_scores: np.ndarray,
        coordinates: np.ndarray,
        threshold: float
    ) -> Dict[str, Any]:
        """Identify senescence clusters using spatial clustering."""

        # High senescence spots
        high_sen_mask = senescence_scores > threshold
        high_sen_coords = coordinates[high_sen_mask]

        if len(high_sen_coords) < self.min_spots_per_region:
            return {
                'n_hotspots': 0,
                'hotspot_labels': np.zeros(len(senescence_scores)),
                'hotspot_centers': []
            }

        # DBSCAN clustering of high senescence spots
        try:
            clustering = DBSCAN(
                eps=self.aging_distance_threshold,
                min_samples=self.min_spots_per_region
            )
            cluster_labels = clustering.fit_predict(high_sen_coords)

            # Map back to full dataset
            full_labels = np.full(len(senescence_scores), -1)
            full_labels[high_sen_mask] = cluster_labels

            # Calculate cluster centers
            unique_labels = np.unique(cluster_labels)
            valid_labels = unique_labels[unique_labels >= 0]

            centers = []
            for label in valid_labels:
                label_mask = cluster_labels == label
                center = np.mean(high_sen_coords[label_mask], axis=0)
                centers.append(center)

            return {
                'n_hotspots': len(valid_labels),
                'hotspot_labels': full_labels,
                'hotspot_centers': centers,
                'high_senescence_fraction': np.mean(high_sen_mask)
            }

        except Exception as e:
            return {
                'error': str(e),
                'n_hotspots': 0,
                'hotspot_labels': np.zeros(len(senescence_scores))
            }

    def _analyze_spatial_gradients(
        self,
        adata,
        age_col: str,
        spatial_key: str
    ) -> Dict[str, Any]:
        """Analyze spatial aging gradients."""

        if self.verbose:
            self.logger.info("Analyzing spatial aging gradients...")

        gradients = {}

        # Get spatial coordinates
        coords = adata.obsm[spatial_key]

        # Age-related gradient analysis
        if age_col in adata.obs.columns:
            age_gradient = self._calculate_age_gradient(
                adata.obs[age_col].values, coords
            )
            gradients['age_gradient'] = age_gradient

        # Expression gradients for aging signatures
        for sig_name, genes in self.spatial_aging_signatures.items():
            available_genes = [g for g in genes if g in adata.var_names]

            if len(available_genes) < 3:
                continue

            # Calculate signature expression
            gene_mask = [g in available_genes for g in adata.var_names]
            sig_expr = adata.X[:, gene_mask]

            if hasattr(sig_expr, 'toarray'):
                sig_expr = sig_expr.toarray()

            sig_score = np.mean(sig_expr, axis=1)

            # Calculate gradient
            gradient = self._calculate_expression_gradient(sig_score, coords)
            gradients[f'{sig_name}_gradient'] = gradient

        return gradients

    def _calculate_age_gradient(
        self,
        age_values: np.ndarray,
        coordinates: np.ndarray
    ) -> Dict[str, float]:
        """Calculate spatial gradient of age effects."""

        if not ML_AVAILABLE:
            return {}

        try:
            # Fit spatial trend using nearest neighbors
            nbrs = NearestNeighbors(n_neighbors=min(50, len(coordinates)))
            nbrs.fit(coordinates)

            gradients = []

            for i in range(len(coordinates)):
                # Get neighbors
                distances, indices = nbrs.kneighbors([coordinates[i]])
                neighbor_ages = age_values[indices[0]]
                neighbor_coords = coordinates[indices[0]]

                # Calculate local gradient
                if len(neighbor_coords) > 2:
                    # Simple gradient estimation
                    x_coords = neighbor_coords[:, 0]
                    y_coords = neighbor_coords[:, 1]

                    # Linear regression for gradient
                    try:
                        dx_corr, _ = pearsonr(x_coords, neighbor_ages)
                        dy_corr, _ = pearsonr(y_coords, neighbor_ages)

                        gradient_magnitude = np.sqrt(dx_corr**2 + dy_corr**2)
                        gradients.append(gradient_magnitude)
                    except Exception:
                        gradients.append(0)
                else:
                    gradients.append(0)

            return {
                'mean_gradient': float(np.mean(gradients)),
                'max_gradient': float(np.max(gradients)),
                'gradient_std': float(np.std(gradients))
            }

        except Exception:
            return {}

    def _calculate_expression_gradient(
        self,
        expression_values: np.ndarray,
        coordinates: np.ndarray
    ) -> Dict[str, float]:
        """Calculate spatial gradient of gene expression."""

        # Similar to age gradient but for expression
        return self._calculate_age_gradient(expression_values, coordinates)

    def _analyze_interaction_aging(
        self,
        adata,
        spatial_key: str,
        tissue_type: str
    ) -> Dict[str, Any]:
        """Analyze cell-cell interaction aging effects."""

        if self.verbose:
            self.logger.info("Analyzing cell-cell interaction aging...")

        if not SPATIAL_AVAILABLE:
            return {'error': 'squidpy not available for interaction analysis'}

        try:
            # Cell-cell interaction analysis using squidpy
            import squidpy as sq

            # Calculate spatial neighbors
            sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key=spatial_key)

            # Interaction matrix
            sq.gr.interaction_matrix(adata, cluster_key="cell_type" if "cell_type" in adata.obs else None)

            # Nhood enrichment
            if "cell_type" in adata.obs.columns:
                sq.gr.nhood_enrichment(adata, cluster_key="cell_type")

                interaction_results = {
                    'has_interactions': True,
                    'n_interactions': adata.uns.get('spatial_neighbors', {}).get('n_neighbors', 0)
                }
            else:
                interaction_results = {'has_interactions': False}

            return interaction_results

        except Exception as e:
            return {'error': str(e)}

    def _analyze_architecture_aging(
        self,
        adata,
        spatial_key: str,
        age_col: str
    ) -> Dict[str, Any]:
        """Analyze tissue architecture aging patterns."""

        if self.verbose:
            self.logger.info("Analyzing tissue architecture aging...")

        architecture = {}

        # Spatial organization metrics
        coords = adata.obsm[spatial_key]

        # Calculate spatial diversity
        spatial_diversity = self._calculate_spatial_diversity(adata, coords)
        architecture['spatial_diversity'] = spatial_diversity

        # Tissue density analysis
        density_analysis = self._analyze_tissue_density(coords)
        architecture['density_analysis'] = density_analysis

        # Age-related architectural changes
        if age_col in adata.obs.columns:
            arch_age_corr = self._correlate_architecture_age(
                adata, coords, age_col
            )
            architecture['age_correlations'] = arch_age_corr

        return architecture

    def _calculate_spatial_diversity(
        self,
        adata,
        coordinates: np.ndarray
    ) -> Dict[str, float]:
        """Calculate spatial expression diversity."""

        # Calculate expression diversity across space
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()

        # Shannon diversity for each spot
        diversities = []
        for i in range(X.shape[0]):
            expr = X[i, :]
            expr = expr + 1e-6  # Avoid log(0)
            expr = expr / np.sum(expr)  # Normalize
            shannon = -np.sum(expr * np.log(expr))
            diversities.append(shannon)

        return {
            'mean_diversity': float(np.mean(diversities)),
            'diversity_std': float(np.std(diversities)),
            'max_diversity': float(np.max(diversities))
        }

    def _analyze_tissue_density(
        self,
        coordinates: np.ndarray
    ) -> Dict[str, float]:
        """Analyze tissue spot density patterns."""

        if not ML_AVAILABLE:
            return {}

        try:
            # Calculate local density using nearest neighbors
            nbrs = NearestNeighbors(n_neighbors=min(20, len(coordinates)))
            nbrs.fit(coordinates)

            distances, _ = nbrs.kneighbors(coordinates)
            local_densities = 1.0 / (np.mean(distances, axis=1) + 1e-6)

            return {
                'mean_density': float(np.mean(local_densities)),
                'density_std': float(np.std(local_densities)),
                'density_cv': float(np.std(local_densities) / np.mean(local_densities))
            }

        except Exception:
            return {}

    def _correlate_architecture_age(
        self,
        adata,
        coordinates: np.ndarray,
        age_col: str
    ) -> Dict[str, float]:
        """Correlate architectural features with age."""

        age_values = adata.obs[age_col].values

        # Spatial density vs age
        density_analysis = self._analyze_tissue_density(coordinates)

        correlations = {}

        if density_analysis and ML_AVAILABLE:
            try:
                nbrs = NearestNeighbors(n_neighbors=min(20, len(coordinates)))
                nbrs.fit(coordinates)
                distances, _ = nbrs.kneighbors(coordinates)
                local_densities = 1.0 / (np.mean(distances, axis=1) + 1e-6)

                # Correlate density with age
                density_age_corr, density_age_p = pearsonr(age_values, local_densities)

                correlations.update({
                    'density_age_correlation': float(density_age_corr),
                    'density_age_p_value': float(density_age_p)
                })

            except Exception:
                pass

        return correlations

    def _generate_spatial_summary(
        self,
        results: Dict[str, Any],
        adata,
        age_col: str
    ) -> Dict[str, Any]:
        """Generate comprehensive spatial aging analysis summary."""

        summary = {
            'dataset_info': {
                'n_spots': int(adata.n_obs),
                'n_genes': int(adata.n_vars),
                'spatial_range_x': [
                    float(adata.obsm['spatial'][:, 0].min()),
                    float(adata.obsm['spatial'][:, 0].max())
                ],
                'spatial_range_y': [
                    float(adata.obsm['spatial'][:, 1].min()),
                    float(adata.obsm['spatial'][:, 1].max())
                ]
            },
            'analysis_parameters': {
                'spatial_resolution': self.spatial_resolution,
                'aging_distance_threshold': self.aging_distance_threshold,
                'senescence_threshold': self.senescence_threshold
            }
        }

        # Add age information if available
        if age_col in adata.obs.columns:
            summary['dataset_info']['age_range'] = [
                float(adata.obs[age_col].min()),
                float(adata.obs[age_col].max())
            ]
            summary['dataset_info']['median_age'] = float(adata.obs[age_col].median())

        # Aging patterns summary
        if 'aging_patterns' in results:
            patterns = results['aging_patterns']
            summary['aging_patterns'] = {
                'n_signatures_analyzed': len(patterns),
                'signatures': list(patterns.keys())
            }

        # Senescence hotspots summary
        if 'senescence_hotspots' in results:
            hotspots = results['senescence_hotspots']
            summary['senescence_hotspots'] = {
                'n_hotspots': hotspots.get('n_hotspots', 0),
                'high_senescence_fraction': hotspots.get('high_senescence_fraction', 0)
            }

        # Spatial gradients summary
        if 'aging_gradients' in results:
            gradients = results['aging_gradients']
            summary['spatial_gradients'] = {
                'n_gradients_analyzed': len(gradients),
                'gradient_types': list(gradients.keys())
            }

        return summary

    def visualize_spatial_aging(
        self,
        adata,
        result_key: str = 'senescence_scores',
        spatial_key: str = 'spatial',
        save_path: Optional[str] = None
    ) -> None:
        """
        Visualize spatial aging patterns.

        Parameters
        ----------
        adata : anndata.AnnData
            Spatial dataset with aging analysis results
        result_key : str
            Key for the aging result to visualize
        spatial_key : str
            Key for spatial coordinates
        save_path : str, optional
            Path to save the plot
        """

        try:
            import matplotlib.pyplot as plt

            if result_key in self.spatial_results.get('senescence_hotspots', {}):
                values = self.spatial_results['senescence_hotspots'][result_key]
                coords = adata.obsm[spatial_key]

                plt.figure(figsize=(10, 8))
                scatter = plt.scatter(
                    coords[:, 0],
                    coords[:, 1],
                    c=values,
                    cmap='viridis',
                    alpha=0.6
                )
                plt.colorbar(scatter, label=result_key)
                plt.xlabel('Spatial X')
                plt.ylabel('Spatial Y')
                plt.title(f'Spatial {result_key.replace("_", " ").title()}')

                if save_path:
                    plt.savefig(save_path, dpi=300, bbox_inches='tight')

                plt.show()

        except ImportError:
            if self.verbose:
                self.logger.warning("matplotlib not available for visualization")


def demo_spatial_aging():
    """Demonstrate spatial aging analysis capabilities."""

    print("🗺️  Spatial Aging Analysis Demo")
    print("=" * 50)

    if not SCANPY_AVAILABLE:
        print("❌ scanpy not available. Install with: pip install scanpy")
        return

    # Create synthetic spatial transcriptomics data
    print("📊 Creating synthetic spatial aging dataset...")

    n_spots = 500
    n_genes = 1000

    # Simulate spatial coordinates
    np.random.seed(42)
    x_coords = np.random.uniform(0, 100, n_spots)
    y_coords = np.random.uniform(0, 100, n_spots)
    spatial_coords = np.column_stack([x_coords, y_coords])

    # Simulate expression data
    X = np.random.negative_binomial(3, 0.3, (n_spots, n_genes))

    # Create gene names with aging markers
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    aging_genes = ['IL6', 'TNF', 'IL1B', 'COL1A1', 'MMP1', 'VCAM1', 'CD68']

    for i, gene in enumerate(aging_genes[:7]):
        if i < n_genes:
            gene_names[i] = gene

    # Add spatial patterns to aging genes
    for i, gene in enumerate(aging_genes[:7]):
        if i < n_genes:
            # Create spatial gradient
            distance_from_center = np.sqrt(
                (x_coords - 50)**2 + (y_coords - 50)**2
            )
            spatial_effect = distance_from_center / 50  # Normalize
            X[:, i] = X[:, i] + np.random.normal(spatial_effect * 3, 0.5, n_spots)
            X[:, i] = np.maximum(X[:, i], 0)  # Ensure non-negative

    # Simulate age information
    ages = np.random.uniform(30, 80, n_spots)

    # Create AnnData object
    adata = ad.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs['chronological_age'] = ages
    adata.obs['sample_id'] = np.random.choice(['sample1', 'sample2'], n_spots)
    adata.obsm['spatial'] = spatial_coords

    print(f"   ✅ Created dataset: {adata.n_obs} spots, {adata.n_vars} genes")
    print(f"   📈 Age range: {ages.min():.1f} - {ages.max():.1f} years")
    print(f"   🗺️  Spatial range: ({x_coords.min():.1f}, {y_coords.min():.1f}) to ({x_coords.max():.1f}, {y_coords.max():.1f})")

    # Initialize spatial aging mapper
    print("\n🔧 Initializing Spatial Aging Mapper...")
    mapper = SpatialAgingMapper(
        spatial_resolution='high',
        aging_distance_threshold=20.0,
        senescence_threshold=2.0,
        verbose=True
    )

    # Run spatial aging analysis
    print("\n🗺️  Running spatial aging analysis...")
    try:
        results = mapper.analyze_spatial_aging(
            adata=adata,
            chronological_age_col='chronological_age',
            spatial_key='spatial',
            tissue_type='general'
        )

        print("\n📊 Analysis Results:")
        print("-" * 30)

        # Summary statistics
        summary = results.get('summary', {})
        dataset_info = summary.get('dataset_info', {})

        print(f"   📊 Dataset: {dataset_info.get('n_spots', 'N/A')} spots")
        print(f"   🧬 Genes: {dataset_info.get('n_genes', 'N/A')}")
        print(f"   📅 Age range: {dataset_info.get('age_range', 'N/A')}")
        print(f"   🗺️  Spatial range X: {dataset_info.get('spatial_range_x', 'N/A')}")
        print(f"   🗺️  Spatial range Y: {dataset_info.get('spatial_range_y', 'N/A')}")

        # Aging patterns
        if 'aging_patterns' in results:
            patterns = results['aging_patterns']
            print("\n🧬 Spatial Aging Patterns:")
            for pattern_name, pattern_data in patterns.items():
                autocorr = pattern_data.get('spatial_autocorr', 0)
                age_corr = pattern_data.get('age_correlation', 0)
                n_genes = pattern_data.get('n_genes', 0)
                print(f"   • {pattern_name}: {n_genes} genes, spatial autocorr: {autocorr:.3f}, age corr: {age_corr:.3f}")

        # Senescence hotspots
        if 'senescence_hotspots' in results:
            hotspots = results['senescence_hotspots']
            n_hotspots = hotspots.get('n_hotspots', 0)
            high_sen_frac = hotspots.get('high_senescence_fraction', 0)
            print("\n🔥 Senescence Hotspots:")
            print(f"   • Number of hotspots: {n_hotspots}")
            print(f"   • High senescence fraction: {high_sen_frac:.3f}")

        # Spatial gradients
        if 'aging_gradients' in results:
            gradients = results['aging_gradients']
            print("\n📈 Spatial Gradients:")
            for gradient_name, gradient_data in gradients.items():
                if isinstance(gradient_data, dict):
                    mean_grad = gradient_data.get('mean_gradient', 0)
                    print(f"   • {gradient_name}: mean gradient {mean_grad:.3f}")

        # Architecture aging
        if 'architecture_aging' in results:
            architecture = results['architecture_aging']
            print("\n🏗️  Tissue Architecture:")

            spatial_div = architecture.get('spatial_diversity', {})
            if spatial_div:
                mean_div = spatial_div.get('mean_diversity', 0)
                print(f"   • Mean spatial diversity: {mean_div:.3f}")

            density = architecture.get('density_analysis', {})
            if density:
                mean_density = density.get('mean_density', 0)
                density_cv = density.get('density_cv', 0)
                print(f"   • Mean density: {mean_density:.3f}")
                print(f"   • Density coefficient of variation: {density_cv:.3f}")

        print("\n✅ Spatial aging analysis completed successfully!")

        # Demonstrate visualization
        print("\n🎨 Testing spatial visualization...")
        try:
            mapper.visualize_spatial_aging(
                adata,
                result_key='senescence_scores',
                spatial_key='spatial'
            )
            print("   ✅ Visualization completed")

        except ImportError:
            print("   ⚠️  matplotlib not available for visualization")
        except Exception as e:
            print(f"   ⚠️  Visualization failed: {e}")

    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        print("   💡 This may be due to missing optional dependencies")
        print("   📦 Install with: pip install squidpy spatialdata")

    print("\n🎉 Demo completed!")
    print("\n💡 Next steps:")
    print("   • Install spatial analysis dependencies for full functionality")
    print("   • Try with real Visium or other spatial datasets")
    print("   • Explore tissue-specific aging patterns")
    print("   • Integrate with single-cell aging analysis")


if __name__ == "__main__":
    demo_spatial_aging()
    demo_spatial_aging()
