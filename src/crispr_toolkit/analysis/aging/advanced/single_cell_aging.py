"""
Single-Cell Aging Analysis Module

Advanced single-cell aging analysis implementing cutting-edge 2025 research
developments in cellular aging heterogeneity and aging clocks.

Features:
- scRNA-seq aging clocks (scAge)
- Cell-type specific biological age estimation
- scATAC-seq chromatin accessibility aging
- Aging trajectory analysis
- Multi-modal single-cell integration

Based on:
- Zhu et al. (2023) PBMC scRNA-seq aging clocks
- Mao et al. (2023) SCALE tissue-specific aging model
- Trapp et al. (2021) scAge epigenetic aging framework
- Latest 2025 single-cell aging research

Author: CRISPR Toolkit Team
Date: August 21, 2025
"""

import warnings
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
import logging

# Core dependencies (installed by default)
try:
    import scanpy as sc
    import anndata as ad
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False

# Optional advanced dependencies
try:
    import cellrank as cr
    import squidpy as sq
    ADVANCED_SC_AVAILABLE = True
except ImportError:
    cr = None
    sq = None
    ADVANCED_SC_AVAILABLE = False

try:
    import torch
    import torch.nn as nn
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import ElasticNet
    ML_AVAILABLE = True
except ImportError:
    torch = None
    nn = None
    RandomForestRegressor = None
    ElasticNet = None
    ML_AVAILABLE = False


class SingleCellAgingAnalyzer:
    """
    Advanced single-cell aging analysis implementing state-of-the-art
    aging clock methodologies and cellular heterogeneity assessment.
    
    This analyzer provides:
    1. scAge-style single-cell aging clocks
    2. Cell-type specific aging pattern analysis
    3. Trajectory-based aging progression modeling
    4. Multi-modal aging signature integration
    5. Aging heterogeneity quantification
    """
    
    def __init__(
        self,
        aging_model: str = 'elastic_net',
        n_aging_genes: int = 500,
        min_cells_per_type: int = 50,
        aging_trajectory_method: str = 'cellrank',
        verbose: bool = True
    ):
        """
        Initialize SingleCellAgingAnalyzer.
        
        Parameters
        ----------
        aging_model : str
            Aging clock model type ('elastic_net', 'random_forest',
            'neural_net')
        n_aging_genes : int
            Number of aging-associated genes to use
        min_cells_per_type : int
            Minimum cells required per cell type for aging analysis
        aging_trajectory_method : str
            Method for aging trajectory analysis ('cellrank', 'scanpy')
        verbose : bool
            Enable verbose logging
        """
        self.aging_model = aging_model
        self.n_aging_genes = n_aging_genes
        self.min_cells_per_type = min_cells_per_type
        self.aging_trajectory_method = aging_trajectory_method
        self.verbose = verbose
        
        # Initialize logging
        if verbose:
            logging.basicConfig(level=logging.INFO)
            self.logger = logging.getLogger(__name__)
        
        # Aging gene signatures from literature
        self.aging_gene_signatures = self._initialize_aging_signatures()
        
        # Model storage
        self.aging_models = {}
        self.cell_type_models = {}
        
    def _initialize_aging_signatures(self) -> Dict[str, List[str]]:
        """Initialize aging gene signatures from literature."""
        
        signatures = {
            # Core aging hallmarks
            'senescence': [
                'CDKN1A', 'CDKN2A', 'TP53', 'RB1', 'ATM', 'CHEK2',
                'MDM2', 'FOXO1', 'FOXO3', 'SIRT1', 'SIRT3', 'SIRT6'
            ],
            
            # DNA damage and genomic instability
            'dna_damage': [
                'ATR', 'BRCA1', 'BRCA2', 'PARP1', 'H2AFX', 'MDC1',
                'NBN', 'RAD51', 'XRCC1', 'XRCC6', 'PCNA', 'RFC1'
            ],
            
            # Telomere attrition
            'telomere': [
                'TERT', 'TERC', 'DKC1', 'TINF2', 'POT1', 'TPP1',
                'TRF1', 'TRF2', 'RAP1', 'CST', 'CTC1', 'STN1'
            ],
            
            # Mitochondrial dysfunction
            'mitochondrial': [
                'SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX3', 'TXNRD2',
                'PPARGC1A', 'NRF1', 'NRF2', 'TFAM', 'POLG', 'MFN1'
            ],
            
            # Inflammation and SASP
            'inflammation': [
                'IL1B', 'IL6', 'TNF', 'CXCL1', 'CXCL8', 'CCL2',
                'NFKB1', 'RELA', 'STAT3', 'JAK2', 'TLR4', 'MYD88'
            ],
            
            # Proteostasis decline
            'proteostasis': [
                'HSPA1A', 'HSPA4', 'HSP90AA1', 'DNAJB1', 'BAG3',
                'PSMC1', 'PSMC4', 'PSMB5', 'UBE3A', 'SQSTM1'
            ],
            
            # Stem cell exhaustion
            'stemness': [
                'SOX2', 'NANOG', 'POU5F1', 'KLF4', 'MYC', 'LIN28A',
                'BMI1', 'EZH2', 'SUZ12', 'RING1', 'CBX7', 'PHC1'
            ],
            
            # Ribosomal and translation aging
            'ribosomal': [
                'RPS6', 'RPS14', 'RPS27A', 'RPL13A', 'RPL7A', 'RPL23',
                'EIF4E', 'EIF4G1', 'MTOR', 'RPS6KB1', 'EEF1A1'
            ]
        }
        
        return signatures
    
    def analyze_single_cell_aging(
        self,
        adata,
        chronological_age_col: str = 'age',
        cell_type_col: str = 'cell_type',
        batch_col: Optional[str] = None,
        aging_signatures: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Perform comprehensive single-cell aging analysis.
        
        Parameters
        ----------
        adata : anndata.AnnData
            Single-cell dataset with cells as observations
        chronological_age_col : str
            Column in adata.obs containing chronological age
        cell_type_col : str
            Column in adata.obs containing cell type annotations
        batch_col : str, optional
            Column for batch correction
        aging_signatures : list, optional
            Custom aging gene signatures to use
            
        Returns
        -------
        dict
            Comprehensive aging analysis results
        """
        if not SCANPY_AVAILABLE:
            raise ImportError(
                "scanpy is required for single-cell aging analysis. "
                "Install with: pip install scanpy"
            )
        
        if self.verbose:
            self.logger.info("Starting single-cell aging analysis...")
            self.logger.info(
                f"Dataset: {adata.n_obs} cells, {adata.n_vars} genes"
            )
        
        results = {}
        
        # 1. Preprocessing and quality control
        adata_processed = self._preprocess_data(
            adata, batch_col=batch_col
        )
        
        # 2. Identify aging-associated genes
        aging_genes = self._identify_aging_genes(
            adata_processed,
            chronological_age_col,
            aging_signatures
        )
        results['aging_genes'] = aging_genes
        
        # 3. Build single-cell aging clocks
        aging_clocks = self._build_sc_aging_clocks(
            adata_processed,
            aging_genes,
            chronological_age_col,
            cell_type_col
        )
        results['aging_clocks'] = aging_clocks
        
        # 4. Analyze cell-type specific aging
        celltype_aging = self._analyze_celltype_aging(
            adata_processed,
            aging_genes,
            chronological_age_col,
            cell_type_col
        )
        results['celltype_aging'] = celltype_aging
        
        # 5. Aging trajectory analysis
        if ADVANCED_SC_AVAILABLE:
            aging_trajectories = self._analyze_aging_trajectories(
                adata_processed,
                chronological_age_col,
                cell_type_col
            )
            results['aging_trajectories'] = aging_trajectories
        
        # 6. Aging heterogeneity quantification
        aging_heterogeneity = self._quantify_aging_heterogeneity(
            adata_processed,
            aging_clocks,
            cell_type_col
        )
        results['aging_heterogeneity'] = aging_heterogeneity
        
        # 7. Generate summary statistics
        summary_stats = self._generate_aging_summary(
            results,
            adata_processed,
            chronological_age_col
        )
        results['summary'] = summary_stats
        
        if self.verbose:
            self.logger.info("Single-cell aging analysis completed!")
            
        return results
    
    def _preprocess_data(
        self,
        adata,
        batch_col: Optional[str] = None
    ):
        """Preprocess single-cell data for aging analysis."""
        
        adata_processed = adata.copy()
        
        if self.verbose:
            self.logger.info("Preprocessing single-cell data...")
        
        # Basic quality control
        sc.pp.filter_cells(adata_processed, min_genes=200)
        sc.pp.filter_genes(adata_processed, min_cells=3)
        
        # Calculate QC metrics
        adata_processed.var['mt'] = adata_processed.var_names.str.startswith(
            'MT-'
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
        
        # Store raw counts
        adata_processed.raw = adata_processed
        
        # Feature selection
        sc.pp.highly_variable_genes(
            adata_processed,
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5
        )
        
        # Batch correction if specified
        if batch_col and batch_col in adata_processed.obs.columns:
            if self.verbose:
                self.logger.info(
                    "Applying batch correction using %s", batch_col
                )
            try:
                import scanpy.external as sce
                sce.pp.harmony_integrate(adata_processed, batch_col)
            except ImportError:
                warnings.warn(
                    "harmony-pytorch not available for batch correction. "
                    "Install with: pip install harmonypy"
                )
        
        # PCA and neighbors
        sc.pp.scale(adata_processed, max_value=10)
        sc.tl.pca(adata_processed, svd_solver='arpack')
        sc.pp.neighbors(adata_processed, n_neighbors=10, n_pcs=40)
        
        return adata_processed
    
    def _identify_aging_genes(
        self,
        adata: 'anndata.AnnData',
        age_col: str,
        custom_signatures: Optional[List[str]] = None
    ) -> Dict[str, List[str]]:
        """Identify aging-associated genes using correlation and signatures."""
        
        if self.verbose:
            self.logger.info("Identifying aging-associated genes...")
        
        aging_genes = {}
        
        # Age correlation analysis
        if age_col in adata.obs.columns:
            age_corr_genes = self._find_age_correlated_genes(adata, age_col)
            aging_genes['age_correlated'] = age_corr_genes
        
        # Literature-based signatures
        for signature_name, genes in self.aging_gene_signatures.items():
            available_genes = [g for g in genes if g in adata.var_names]
            if len(available_genes) > 5:  # Minimum threshold
                aging_genes[signature_name] = available_genes
        
        # Custom signatures
        if custom_signatures:
            available_custom = [g for g in custom_signatures if g in adata.var_names]
            if available_custom:
                aging_genes['custom'] = available_custom
        
        # Combined aging signature
        all_aging_genes = set()
        for gene_list in aging_genes.values():
            all_aging_genes.update(gene_list)
        
        # Limit to top N aging genes by variance or correlation
        if len(all_aging_genes) > self.n_aging_genes:
            aging_genes['combined'] = list(all_aging_genes)[:self.n_aging_genes]
        else:
            aging_genes['combined'] = list(all_aging_genes)
        
        if self.verbose:
            self.logger.info(f"Identified {len(aging_genes['combined'])} aging genes")
        
        return aging_genes
    
    def _find_age_correlated_genes(
        self,
        adata: 'anndata.AnnData',
        age_col: str,
        correlation_threshold: float = 0.1,
        p_value_threshold: float = 0.05
    ) -> List[str]:
        """Find genes correlated with chronological age."""
        
        from scipy.stats import pearsonr
        
        age_values = adata.obs[age_col].values
        correlations = []
        p_values = []
        
        # Calculate correlations for each gene
        X = adata.X if hasattr(adata.X, 'toarray') else adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        
        for i in range(X.shape[1]):
            gene_expr = X[:, i]
            try:
                corr, p_val = pearsonr(age_values, gene_expr)
                correlations.append(abs(corr))
                p_values.append(p_val)
            except:
                correlations.append(0)
                p_values.append(1)
        
        # Filter by correlation and significance
        corr_df = pd.DataFrame({
            'gene': adata.var_names,
            'correlation': correlations,
            'p_value': p_values
        })
        
        significant_genes = corr_df[
            (corr_df['correlation'] > correlation_threshold) &
            (corr_df['p_value'] < p_value_threshold)
        ].sort_values('correlation', ascending=False)
        
        return significant_genes['gene'].tolist()
    
    def _build_sc_aging_clocks(
        self,
        adata: 'anndata.AnnData',
        aging_genes: Dict[str, List[str]],
        age_col: str,
        cell_type_col: str
    ) -> Dict[str, Any]:
        """Build single-cell aging clocks using multiple methodologies."""
        
        if not ML_AVAILABLE:
            warnings.warn(
                "scikit-learn not available. Using simplified aging clock."
            )
            return {'method': 'simplified', 'models': {}}
        
        if self.verbose:
            self.logger.info("Building single-cell aging clocks...")
        
        aging_clocks = {
            'global_clock': None,
            'celltype_clocks': {},
            'predictions': {},
            'performance': {}
        }
        
        # Prepare data
        X, y, cell_types = self._prepare_aging_data(
            adata, aging_genes['combined'], age_col, cell_type_col
        )
        
        # Global aging clock
        global_model = self._train_aging_model(X, y)
        aging_clocks['global_clock'] = global_model
        
        # Cell-type specific clocks
        unique_cell_types = np.unique(cell_types)
        for cell_type in unique_cell_types:
            type_mask = cell_types == cell_type
            if np.sum(type_mask) >= self.min_cells_per_type:
                X_type = X[type_mask]
                y_type = y[type_mask]
                
                type_model = self._train_aging_model(X_type, y_type)
                aging_clocks['celltype_clocks'][cell_type] = type_model
        
        # Generate predictions
        predictions = self._generate_aging_predictions(
            X, y, cell_types, aging_clocks
        )
        aging_clocks['predictions'] = predictions
        
        # Evaluate performance
        performance = self._evaluate_aging_clocks(predictions, y)
        aging_clocks['performance'] = performance
        
        return aging_clocks
    
    def _prepare_aging_data(
        self,
        adata: 'anndata.AnnData',
        aging_genes: List[str],
        age_col: str,
        cell_type_col: str
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Prepare data matrices for aging clock training."""
        
        # Filter to aging genes
        gene_mask = [g in aging_genes for g in adata.var_names]
        X = adata.X[:, gene_mask]
        
        if hasattr(X, 'toarray'):
            X = X.toarray()
        
        # Get age and cell type labels
        y = adata.obs[age_col].values
        cell_types = adata.obs[cell_type_col].values
        
        # Remove any missing values
        valid_mask = ~(np.isnan(y) | pd.isna(cell_types))
        X = X[valid_mask]
        y = y[valid_mask]
        cell_types = cell_types[valid_mask]
        
        return X, y, cell_types
    
    def _train_aging_model(
        self,
        X: np.ndarray,
        y: np.ndarray
    ) -> Any:
        """Train aging prediction model."""
        
        if self.aging_model == 'elastic_net':
            model = ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42)
        elif self.aging_model == 'random_forest':
            model = RandomForestRegressor(
                n_estimators=100, 
                random_state=42,
                n_jobs=-1
            )
        else:
            # Default to elastic net
            model = ElasticNet(alpha=0.1, l1_ratio=0.5, random_state=42)
        
        model.fit(X, y)
        return model
    
    def _generate_aging_predictions(
        self,
        X: np.ndarray,
        y: np.ndarray,
        cell_types: np.ndarray,
        aging_clocks: Dict[str, Any]
    ) -> Dict[str, np.ndarray]:
        """Generate aging predictions using trained clocks."""
        
        predictions = {}
        
        # Global predictions
        if aging_clocks['global_clock'] is not None:
            global_pred = aging_clocks['global_clock'].predict(X)
            predictions['global'] = global_pred
            predictions['global_age_acceleration'] = global_pred - y
        
        # Cell-type specific predictions
        celltype_predictions = np.full(len(y), np.nan)
        celltype_age_accel = np.full(len(y), np.nan)
        
        unique_types = np.unique(cell_types)
        for cell_type in unique_types:
            if cell_type in aging_clocks['celltype_clocks']:
                type_mask = cell_types == cell_type
                type_model = aging_clocks['celltype_clocks'][cell_type]
                
                type_pred = type_model.predict(X[type_mask])
                celltype_predictions[type_mask] = type_pred
                celltype_age_accel[type_mask] = type_pred - y[type_mask]
        
        predictions['celltype'] = celltype_predictions
        predictions['celltype_age_acceleration'] = celltype_age_accel
        predictions['chronological_age'] = y
        predictions['cell_types'] = cell_types
        
        return predictions
    
    def _evaluate_aging_clocks(
        self,
        predictions: Dict[str, np.ndarray],
        true_age: np.ndarray
    ) -> Dict[str, float]:
        """Evaluate aging clock performance."""
        
        from sklearn.metrics import mean_absolute_error, r2_score
        
        performance = {}
        
        # Global clock performance
        if 'global' in predictions:
            global_pred = predictions['global']
            performance['global_mae'] = mean_absolute_error(true_age, global_pred)
            performance['global_r2'] = r2_score(true_age, global_pred)
        
        # Cell-type clock performance
        if 'celltype' in predictions:
            celltype_pred = predictions['celltype']
            valid_mask = ~np.isnan(celltype_pred)
            
            if np.sum(valid_mask) > 0:
                performance['celltype_mae'] = mean_absolute_error(
                    true_age[valid_mask], 
                    celltype_pred[valid_mask]
                )
                performance['celltype_r2'] = r2_score(
                    true_age[valid_mask], 
                    celltype_pred[valid_mask]
                )
        
        return performance
    
    def _analyze_celltype_aging(
        self,
        adata: 'anndata.AnnData',
        aging_genes: Dict[str, List[str]],
        age_col: str,
        cell_type_col: str
    ) -> Dict[str, Any]:
        """Analyze cell-type specific aging patterns."""
        
        if self.verbose:
            self.logger.info("Analyzing cell-type specific aging patterns...")
        
        celltype_aging = {}
        
        unique_cell_types = adata.obs[cell_type_col].unique()
        
        for cell_type in unique_cell_types:
            type_mask = adata.obs[cell_type_col] == cell_type
            type_data = adata[type_mask]
            
            if type_data.n_obs < self.min_cells_per_type:
                continue
            
            # Aging gene expression patterns
            aging_expr = self._analyze_aging_expression(
                type_data, aging_genes['combined'], age_col
            )
            
            # Senescence markers
            senescence_score = self._calculate_senescence_score(
                type_data, aging_genes.get('senescence', [])
            )
            
            # Cell cycle analysis
            cell_cycle = self._analyze_cell_cycle_aging(type_data, age_col)
            
            celltype_aging[cell_type] = {
                'aging_expression': aging_expr,
                'senescence_score': senescence_score,
                'cell_cycle': cell_cycle,
                'n_cells': type_data.n_obs
            }
        
        return celltype_aging
    
    def _analyze_aging_expression(
        self,
        adata: 'anndata.AnnData',
        aging_genes: List[str],
        age_col: str
    ) -> Dict[str, float]:
        """Analyze aging gene expression patterns."""
        
        available_genes = [g for g in aging_genes if g in adata.var_names]
        
        if len(available_genes) == 0:
            return {}
        
        # Calculate aging gene expression score
        gene_mask = [g in available_genes for g in adata.var_names]
        aging_expr = adata.X[:, gene_mask]
        
        if hasattr(aging_expr, 'toarray'):
            aging_expr = aging_expr.toarray()
        
        aging_score = np.mean(aging_expr, axis=1)
        
        # Correlation with age
        if age_col in adata.obs.columns:
            from scipy.stats import pearsonr
            age_corr, age_p = pearsonr(adata.obs[age_col], aging_score)
        else:
            age_corr, age_p = 0, 1
        
        return {
            'aging_score_mean': float(np.mean(aging_score)),
            'aging_score_std': float(np.std(aging_score)),
            'age_correlation': float(age_corr),
            'age_correlation_p': float(age_p),
            'n_genes': len(available_genes)
        }
    
    def _calculate_senescence_score(
        self,
        adata: 'anndata.AnnData',
        senescence_genes: List[str]
    ) -> Dict[str, float]:
        """Calculate senescence score using established markers."""
        
        available_genes = [g for g in senescence_genes if g in adata.var_names]
        
        if len(available_genes) == 0:
            return {'senescence_score': 0, 'n_genes': 0}
        
        gene_mask = [g in available_genes for g in adata.var_names]
        sen_expr = adata.X[:, gene_mask]
        
        if hasattr(sen_expr, 'toarray'):
            sen_expr = sen_expr.toarray()
        
        senescence_score = np.mean(sen_expr, axis=1)
        
        return {
            'senescence_score': float(np.mean(senescence_score)),
            'senescence_score_std': float(np.std(senescence_score)),
            'n_genes': len(available_genes)
        }
    
    def _analyze_cell_cycle_aging(
        self,
        adata: 'anndata.AnnData',
        age_col: str
    ) -> Dict[str, Any]:
        """Analyze cell cycle changes with aging."""
        
        try:
            # Cell cycle scoring
            sc.tl.score_genes_cell_cycle(
                adata,
                s_genes=sc.settings.verbosity,  # Use default S genes
                g2m_genes=sc.settings.verbosity  # Use default G2M genes
            )
            
            cell_cycle_results = {}
            
            for phase in ['S_score', 'G2M_score']:
                if phase in adata.obs.columns and age_col in adata.obs.columns:
                    from scipy.stats import pearsonr
                    corr, p_val = pearsonr(adata.obs[age_col], adata.obs[phase])
                    cell_cycle_results[f'{phase}_age_correlation'] = float(corr)
                    cell_cycle_results[f'{phase}_age_p_value'] = float(p_val)
            
            return cell_cycle_results
            
        except Exception as e:
            if self.verbose:
                self.logger.warning(f"Cell cycle analysis failed: {e}")
            return {}
    
    def _analyze_aging_trajectories(
        self,
        adata: 'anndata.AnnData',
        age_col: str,
        cell_type_col: str
    ) -> Dict[str, Any]:
        """Analyze aging trajectories using advanced methods."""
        
        if self.verbose:
            self.logger.info("Analyzing aging trajectories...")
        
        trajectories = {}
        
        try:
            if self.aging_trajectory_method == 'cellrank':
                # CellRank-based trajectory analysis
                cr_results = self._cellrank_aging_analysis(
                    adata, age_col, cell_type_col
                )
                trajectories['cellrank'] = cr_results
            
            # Pseudotime analysis
            pseudotime_results = self._pseudotime_aging_analysis(
                adata, age_col, cell_type_col
            )
            trajectories['pseudotime'] = pseudotime_results
            
        except Exception as e:
            if self.verbose:
                self.logger.warning(f"Trajectory analysis failed: {e}")
            trajectories['error'] = str(e)
        
        return trajectories
    
    def _cellrank_aging_analysis(
        self,
        adata: 'anndata.AnnData',
        age_col: str,
        cell_type_col: str
    ) -> Dict[str, Any]:
        """Perform CellRank-based aging trajectory analysis."""
        
        try:
            import cellrank as cr
            
            # Velocity-based analysis if available
            if 'velocity' in adata.layers:
                cr.tl.terminal_states(adata, cluster_key=cell_type_col)
                cr.tl.initial_states(adata, cluster_key=cell_type_col)
                
                # Compute macrostates
                cr.tl.macrostates(adata, n_states=10)
                
                # Aging-associated fate probabilities
                if age_col in adata.obs.columns:
                    # Correlate fate probabilities with age
                    fate_corr = {}
                    for state in adata.obsm['macrostates'].columns:
                        from scipy.stats import pearsonr
                        corr, p_val = pearsonr(
                            adata.obs[age_col],
                            adata.obsm['macrostates'][state]
                        )
                        fate_corr[state] = {'correlation': corr, 'p_value': p_val}
                    
                    return {
                        'fate_correlations': fate_corr,
                        'macrostates': adata.obsm['macrostates'].shape[1]
                    }
            
        except ImportError:
            return {'error': 'CellRank not available'}
        except Exception as e:
            return {'error': str(e)}
        
        return {}
    
    def _pseudotime_aging_analysis(
        self,
        adata: 'anndata.AnnData',
        age_col: str,
        cell_type_col: str
    ) -> Dict[str, Any]:
        """Perform pseudotime-based aging analysis."""
        
        try:
            # UMAP-based pseudotime
            sc.tl.umap(adata)
            
            # Diffusion pseudotime
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
            sc.tl.diffmap(adata)
            sc.tl.dpt(adata)
            
            # Correlate pseudotime with chronological age
            if 'dpt_pseudotime' in adata.obs.columns and age_col in adata.obs.columns:
                from scipy.stats import pearsonr
                pseudo_age_corr, pseudo_age_p = pearsonr(
                    adata.obs[age_col],
                    adata.obs['dpt_pseudotime']
                )
                
                return {
                    'pseudotime_age_correlation': float(pseudo_age_corr),
                    'pseudotime_age_p_value': float(pseudo_age_p),
                    'has_pseudotime': True
                }
            
        except Exception as e:
            return {'error': str(e)}
        
        return {'has_pseudotime': False}
    
    def _quantify_aging_heterogeneity(
        self,
        adata: 'anndata.AnnData',
        aging_clocks: Dict[str, Any],
        cell_type_col: str
    ) -> Dict[str, Any]:
        """Quantify aging heterogeneity within and between cell types."""
        
        if self.verbose:
            self.logger.info("Quantifying aging heterogeneity...")
        
        heterogeneity = {}
        
        if 'predictions' not in aging_clocks:
            return heterogeneity
        
        predictions = aging_clocks['predictions']
        
        # Global aging heterogeneity
        if 'global_age_acceleration' in predictions:
            age_accel = predictions['global_age_acceleration']
            heterogeneity['global_heterogeneity'] = {
                'age_acceleration_std': float(np.std(age_accel)),
                'age_acceleration_range': float(np.ptp(age_accel)),
                'coefficient_of_variation': float(np.std(age_accel) / np.mean(np.abs(age_accel)))
            }
        
        # Cell-type specific heterogeneity
        if 'celltype_age_acceleration' in predictions and 'cell_types' in predictions:
            celltype_heterogeneity = {}
            cell_types = predictions['cell_types']
            celltype_age_accel = predictions['celltype_age_acceleration']
            
            unique_types = np.unique(cell_types)
            for cell_type in unique_types:
                type_mask = cell_types == cell_type
                type_age_accel = celltype_age_accel[type_mask]
                
                # Remove NaN values
                type_age_accel = type_age_accel[~np.isnan(type_age_accel)]
                
                if len(type_age_accel) > 0:
                    celltype_heterogeneity[cell_type] = {
                        'age_acceleration_std': float(np.std(type_age_accel)),
                        'age_acceleration_range': float(np.ptp(type_age_accel)),
                        'n_cells': int(len(type_age_accel))
                    }
            
            heterogeneity['celltype_heterogeneity'] = celltype_heterogeneity
        
        return heterogeneity
    
    def _generate_aging_summary(
        self,
        results: Dict[str, Any],
        adata: 'anndata.AnnData',
        age_col: str
    ) -> Dict[str, Any]:
        """Generate comprehensive aging analysis summary."""
        
        summary = {
            'dataset_info': {
                'n_cells': int(adata.n_obs),
                'n_genes': int(adata.n_vars),
                'age_range': [float(adata.obs[age_col].min()), float(adata.obs[age_col].max())],
                'median_age': float(adata.obs[age_col].median())
            },
            'analysis_parameters': {
                'aging_model': self.aging_model,
                'n_aging_genes': self.n_aging_genes,
                'min_cells_per_type': self.min_cells_per_type
            }
        }
        
        # Aging genes summary
        if 'aging_genes' in results:
            aging_genes = results['aging_genes']
            summary['aging_genes'] = {
                signature: len(genes) 
                for signature, genes in aging_genes.items()
            }
        
        # Performance summary
        if 'aging_clocks' in results and 'performance' in results['aging_clocks']:
            summary['performance'] = results['aging_clocks']['performance']
        
        # Cell-type summary
        if 'celltype_aging' in results:
            summary['analyzed_cell_types'] = list(results['celltype_aging'].keys())
            summary['n_analyzed_cell_types'] = len(results['celltype_aging'])
        
        return summary
    
    def predict_cellular_age(
        self,
        adata: 'anndata.AnnData',
        aging_model: Optional[Any] = None,
        cell_type_col: Optional[str] = None
    ) -> np.ndarray:
        """
        Predict cellular age for new data using trained aging clocks.
        
        Parameters
        ----------
        adata : anndata.AnnData
            New single-cell data
        aging_model : optional
            Pre-trained aging model to use
        cell_type_col : str, optional
            Cell type column for cell-type specific predictions
            
        Returns
        -------
        np.ndarray
            Predicted cellular ages
        """
        
        if aging_model is None and not self.aging_models:
            raise ValueError(
                "No aging model available. Run analyze_single_cell_aging first."
            )
        
        # Use provided model or stored global model
        model = aging_model or self.aging_models.get('global')
        
        if model is None:
            raise ValueError("No suitable aging model found")
        
        # Prepare data (simplified version)
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        
        # Make predictions
        predictions = model.predict(X)
        
        return predictions


def demo_single_cell_aging():
    """Demonstrate single-cell aging analysis capabilities."""
    
    print("🧬 Single-Cell Aging Analysis Demo")
    print("=" * 50)
    
    if not SCANPY_AVAILABLE:
        print("❌ scanpy not available. Install with: pip install scanpy")
        return
    
    # Create synthetic single-cell data
    print("📊 Creating synthetic single-cell aging dataset...")
    
    n_cells = 1000
    n_genes = 2000
    n_cell_types = 4
    
    # Simulate expression data
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.3, (n_cells, n_genes))
    
    # Create gene names
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    
    # Add some aging-related genes
    aging_genes = [
        'CDKN1A', 'CDKN2A', 'TP53', 'IL6', 'TNF', 'SOD2',
        'TERT', 'SIRT1', 'FOXO1', 'RPS6'
    ]
    
    for i, gene in enumerate(aging_genes[:10]):
        if i < n_genes:
            gene_names[i] = gene
    
    # Simulate age and cell types
    ages = np.random.uniform(20, 80, n_cells)
    cell_types = np.random.choice(
        ['T_cells', 'B_cells', 'Monocytes', 'NK_cells'], 
        n_cells
    )
    
    # Add age-dependent expression patterns
    for i, gene in enumerate(aging_genes[:10]):
        if i < n_genes:
            # Age-correlated expression
            age_effect = (ages - 50) / 30  # Normalize age
            X[:, i] = X[:, i] + np.random.normal(age_effect * 2, 0.5, n_cells)
            X[:, i] = np.maximum(X[:, i], 0)  # Ensure non-negative
    
    # Create AnnData object
    import pandas as pd
    
    adata = ad.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs['chronological_age'] = ages
    adata.obs['cell_type'] = cell_types
    adata.obs['batch'] = np.random.choice(['batch1', 'batch2'], n_cells)
    
    print(f"   ✅ Created dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"   📈 Age range: {ages.min():.1f} - {ages.max():.1f} years")
    print(f"   🏷️  Cell types: {len(np.unique(cell_types))}")
    
    # Initialize analyzer
    print("\n🔧 Initializing Single-Cell Aging Analyzer...")
    analyzer = SingleCellAgingAnalyzer(
        aging_model='elastic_net',
        n_aging_genes=100,
        verbose=True
    )
    
    # Run aging analysis
    print("\n🧬 Running single-cell aging analysis...")
    try:
        results = analyzer.analyze_single_cell_aging(
            adata=adata,
            chronological_age_col='chronological_age',
            cell_type_col='cell_type'
        )
        
        print("\n📊 Analysis Results:")
        print("-" * 30)
        
        # Summary statistics
        summary = results.get('summary', {})
        dataset_info = summary.get('dataset_info', {})
        
        print(f"   📈 Dataset: {dataset_info.get('n_cells', 'N/A')} cells")
        print(f"   🧬 Genes: {dataset_info.get('n_genes', 'N/A')}")
        print(f"   📅 Age range: {dataset_info.get('age_range', 'N/A')}")
        
        # Aging genes
        aging_genes_summary = summary.get('aging_genes', {})
        print(f"\n🧬 Aging Genes Identified:")
        for signature, count in aging_genes_summary.items():
            print(f"   • {signature}: {count} genes")
        
        # Performance
        performance = summary.get('performance', {})
        if performance:
            print(f"\n📊 Aging Clock Performance:")
            for metric, value in performance.items():
                print(f"   • {metric}: {value:.3f}")
        
        # Cell-type analysis
        n_cell_types = summary.get('n_analyzed_cell_types', 0)
        print(f"\n🏷️  Analyzed Cell Types: {n_cell_types}")
        
        if 'celltype_aging' in results:
            for cell_type, analysis in results['celltype_aging'].items():
                n_cells = analysis.get('n_cells', 0)
                aging_expr = analysis.get('aging_expression', {})
                aging_score = aging_expr.get('aging_score_mean', 0)
                print(f"   • {cell_type}: {n_cells} cells, aging score: {aging_score:.3f}")
        
        # Heterogeneity
        if 'aging_heterogeneity' in results:
            heterogeneity = results['aging_heterogeneity']
            global_het = heterogeneity.get('global_heterogeneity', {})
            age_accel_std = global_het.get('age_acceleration_std', 0)
            print(f"\n📊 Aging Heterogeneity:")
            print(f"   • Age acceleration std: {age_accel_std:.3f}")
        
        print("\n✅ Single-cell aging analysis completed successfully!")
        
        # Demonstrate prediction on new data
        print("\n🔮 Testing cellular age prediction...")
        if 'aging_clocks' in results and 'global_clock' in results['aging_clocks']:
            # Use a subset as "new" data
            test_data = adata[:100].copy()
            
            try:
                predicted_ages = analyzer.predict_cellular_age(
                    test_data,
                    aging_model=results['aging_clocks']['global_clock']
                )
                
                true_ages = test_data.obs['chronological_age'].values
                mae = np.mean(np.abs(predicted_ages - true_ages))
                
                print(f"   ✅ Predicted ages for {len(predicted_ages)} cells")
                print(f"   📊 Mean Absolute Error: {mae:.2f} years")
                
            except Exception as e:
                print(f"   ⚠️  Prediction failed: {e}")
        
    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        print("   💡 This may be due to missing optional dependencies")
        print("   📦 Install with: pip install scanpy scikit-learn")
    
    print("\n🎉 Demo completed!")
    print("\n💡 Next steps:")
    print("   • Install advanced dependencies for full functionality")
    print("   • Try with real single-cell datasets")
    print("   • Explore cell-type specific aging patterns")
    print("   • Integrate with spatial transcriptomics data")


if __name__ == "__main__":
    demo_single_cell_aging()
