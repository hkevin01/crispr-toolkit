"""
Aging Biomarker Visualization Module

Provides comprehensive visualization tools for aging biomarker analysis
with specific focus on ReHMGB1/RAGE signaling pathways and senescence
research applications.

Visualization Features:
- Aging clock comparison plots
- Age acceleration analysis
- Pathway-specific aging patterns
- Clinical risk stratification plots
- ReHMGB1 pathway network visualization
- Longitudinal aging trajectory plots

Author: CRISPR Toolkit Development Team
"""

import warnings
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

# Handle optional plotting dependencies
try:
    import matplotlib.pyplot as plt
    import seaborn as sns
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    warnings.warn("matplotlib/seaborn not available - install for visualizations")

try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    warnings.warn("plotly not available - install for interactive plots")


class AgingBiomarkerVisualizer:
    """
    Comprehensive visualization toolkit for aging biomarker analysis.

    Provides both static (matplotlib/seaborn) and interactive (plotly)
    visualizations with ReHMGB1/senescence pathway focus.
    """

    def __init__(self, style: str = 'scientific', verbose: bool = True):
        """
        Initialize aging biomarker visualizer.

        Args:
            style: Plotting style ('scientific', 'clinical', 'presentation')
            verbose: Enable verbose logging
        """
        self.verbose = verbose
        self.style = style

        # Set up plotting parameters
        self._setup_plotting_style()

        # Define color schemes
        self.color_schemes = self._load_color_schemes()

        if verbose:
            available_backends = []
            if MATPLOTLIB_AVAILABLE:
                available_backends.append('matplotlib')
            if PLOTLY_AVAILABLE:
                available_backends.append('plotly')

            print("📊 Visualization tools initialized")
            print(f"   Available backends: {available_backends}")

    def _setup_plotting_style(self):
        """Setup plotting style and parameters."""
        if MATPLOTLIB_AVAILABLE:
            # Set style based on preference
            if self.style == 'scientific':
                plt.style.use('seaborn-v0_8-whitegrid' if hasattr(plt.style, 'available')
                             and 'seaborn-v0_8-whitegrid' in plt.style.available
                             else 'default')
                sns.set_palette("husl")
            elif self.style == 'clinical':
                plt.style.use('seaborn-v0_8-white' if hasattr(plt.style, 'available')
                             and 'seaborn-v0_8-white' in plt.style.available
                             else 'default')
                sns.set_palette("Set2")
            else:  # presentation
                plt.style.use('seaborn-v0_8-talk' if hasattr(plt.style, 'available')
                             and 'seaborn-v0_8-talk' in plt.style.available
                             else 'default')
                sns.set_palette("bright")

    def _load_color_schemes(self) -> Dict[str, Dict]:
        """Load color schemes for different plot types."""
        return {
            'risk_categories': {
                'low': '#2E8B57',      # Sea green
                'moderate': '#FFD700',  # Gold
                'high': '#FF6347',      # Tomato
                'very_high': '#DC143C'  # Crimson
            },
            'aging_clocks': {
                'Horvath2013': '#1f77b4',
                'Hannum2013': '#ff7f0e',
                'PhenoAge': '#2ca02c',
                'GrimAge': '#d62728',
                'DunedinPACE': '#9467bd'
            },
            'pathways': {
                'hmgb1_core': '#8B0000',
                'rage_signaling': '#FF4500',
                'jak_stat_pathway': '#4169E1',
                'nfkb_signaling': '#32CD32',
                'sasp_factors': '#FF1493',
                'oxidative_stress': '#8A2BE2'
            },
            'senescence': {
                'minimal': '#90EE90',
                'low': '#32CD32',
                'moderate': '#FFD700',
                'high': '#FF6347',
                'severe': '#8B0000'
            }
        }

    def plot_aging_clock_comparison(
        self,
        aging_predictions: pd.DataFrame,
        chronological_age: Optional[pd.Series] = None,
        save_path: Optional[str] = None,
        interactive: bool = False
    ) -> Optional[Any]:
        """
        Create aging clock comparison visualization.

        Args:
            aging_predictions: DataFrame with aging clock predictions
            chronological_age: Series with chronological ages
            save_path: Path to save plot
            interactive: Use interactive plotly plot

        Returns:
            Plot object or None
        """
        if interactive and PLOTLY_AVAILABLE:
            return self._plot_clock_comparison_plotly(
                aging_predictions, chronological_age, save_path
            )
        elif MATPLOTLIB_AVAILABLE:
            return self._plot_clock_comparison_matplotlib(
                aging_predictions, chronological_age, save_path
            )
        else:
            if self.verbose:
                print("⚠️  No plotting libraries available")
            return None

    def _plot_clock_comparison_matplotlib(
        self,
        aging_predictions: pd.DataFrame,
        chronological_age: Optional[pd.Series],
        save_path: Optional[str]
    ):
        """Create clock comparison plot with matplotlib."""

        # Determine prediction column
        pred_col = 'predicted_age' if 'predicted_age' in aging_predictions.columns else 'predicted_value'

        # Create figure
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))

        # Plot 1: Clock predictions by sample
        pivot_data = aging_predictions.pivot(
            index='sample_id',
            columns='clock',
            values=pred_col
        )

        # Box plot of predictions by clock
        ax1 = axes[0]
        box_data = [pivot_data[col].dropna() for col in pivot_data.columns]
        box_labels = list(pivot_data.columns)

        bp = ax1.boxplot(box_data, labels=box_labels, patch_artist=True)

        # Color boxes
        colors = [self.color_schemes['aging_clocks'].get(label, '#1f77b4')
                 for label in box_labels]
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        ax1.set_title('Aging Clock Predictions Distribution')
        ax1.set_ylabel('Predicted Age (years)')
        ax1.set_xlabel('Aging Clock')
        ax1.tick_params(axis='x', rotation=45)

        # Plot 2: Age acceleration if chronological age available
        if chronological_age is not None:
            ax2 = axes[1]

            # Calculate age acceleration
            aging_with_age = aging_predictions.copy()
            aging_with_age['chronological_age'] = aging_with_age['sample_id'].map(
                chronological_age.to_dict()
            )
            aging_with_age['age_acceleration'] = (
                aging_with_age[pred_col] - aging_with_age['chronological_age']
            )

            # Box plot of age acceleration by clock
            acceleration_pivot = aging_with_age.pivot(
                index='sample_id',
                columns='clock',
                values='age_acceleration'
            )

            acc_data = [acceleration_pivot[col].dropna() for col in acceleration_pivot.columns]
            acc_labels = list(acceleration_pivot.columns)

            bp2 = ax2.boxplot(acc_data, labels=acc_labels, patch_artist=True)

            # Color boxes
            for patch, color in zip(bp2['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ax2.axhline(y=0, color='red', linestyle='--', alpha=0.7)
            ax2.set_title('Age Acceleration by Clock')
            ax2.set_ylabel('Age Acceleration (years)')
            ax2.set_xlabel('Aging Clock')
            ax2.tick_params(axis='x', rotation=45)
        else:
            # Plot correlation matrix if multiple clocks
            ax2 = axes[1]
            if len(pivot_data.columns) > 1:
                corr_matrix = pivot_data.corr()
                im = ax2.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)

                # Add correlation values
                for i in range(len(corr_matrix)):
                    for j in range(len(corr_matrix)):
                        text = ax2.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                                       ha="center", va="center", color="black")

                ax2.set_xticks(range(len(corr_matrix.columns)))
                ax2.set_yticks(range(len(corr_matrix.columns)))
                ax2.set_xticklabels(corr_matrix.columns, rotation=45)
                ax2.set_yticklabels(corr_matrix.columns)
                ax2.set_title('Clock Correlation Matrix')

                # Add colorbar
                plt.colorbar(im, ax=ax2, shrink=0.8)
            else:
                ax2.text(0.5, 0.5, 'Need multiple clocks\nfor correlation analysis',
                        ha='center', va='center', transform=ax2.transAxes)
                ax2.set_title('Clock Correlations')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            if self.verbose:
                print(f"   💾 Plot saved to: {save_path}")

        return fig

    def _plot_clock_comparison_plotly(
        self,
        aging_predictions: pd.DataFrame,
        chronological_age: Optional[pd.Series],
        save_path: Optional[str]
    ):
        """Create interactive clock comparison plot with plotly."""

        pred_col = 'predicted_age' if 'predicted_age' in aging_predictions.columns else 'predicted_value'

        # Create subplot figure
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=('Clock Predictions', 'Age Acceleration'),
            horizontal_spacing=0.1
        )

        # Plot 1: Box plots of predictions
        clocks = aging_predictions['clock'].unique()
        for i, clock in enumerate(clocks):
            clock_data = aging_predictions[aging_predictions['clock'] == clock]
            color = self.color_schemes['aging_clocks'].get(clock, '#1f77b4')

            fig.add_trace(
                go.Box(
                    y=clock_data[pred_col],
                    name=clock,
                    marker_color=color,
                    boxpoints='outliers'
                ),
                row=1, col=1
            )

        # Plot 2: Age acceleration if available
        if chronological_age is not None:
            aging_with_age = aging_predictions.copy()
            aging_with_age['chronological_age'] = aging_with_age['sample_id'].map(
                chronological_age.to_dict()
            )
            aging_with_age['age_acceleration'] = (
                aging_with_age[pred_col] - aging_with_age['chronological_age']
            )

            for i, clock in enumerate(clocks):
                clock_data = aging_with_age[aging_with_age['clock'] == clock]
                color = self.color_schemes['aging_clocks'].get(clock, '#1f77b4')

                fig.add_trace(
                    go.Box(
                        y=clock_data['age_acceleration'],
                        name=clock,
                        marker_color=color,
                        showlegend=False,
                        boxpoints='outliers'
                    ),
                    row=1, col=2
                )

            # Add horizontal line at y=0
            fig.add_hline(y=0, line_dash="dash", line_color="red",
                         row=1, col=2, opacity=0.7)

        # Update layout
        fig.update_layout(
            title="Aging Clock Analysis",
            height=500,
            showlegend=True
        )

        fig.update_yaxes(title_text="Predicted Age (years)", row=1, col=1)
        fig.update_yaxes(title_text="Age Acceleration (years)", row=1, col=2)
        fig.update_xaxes(title_text="Aging Clock", row=1, col=1)
        fig.update_xaxes(title_text="Aging Clock", row=1, col=2)

        if save_path:
            fig.write_html(save_path)
            if self.verbose:
                print(f"   💾 Interactive plot saved to: {save_path}")

        return fig

    def plot_rehmgb1_pathway_heatmap(
        self,
        pathway_scores: Dict[str, float],
        save_path: Optional[str] = None,
        interactive: bool = False
    ) -> Optional[Any]:
        """
        Create ReHMGB1 pathway scoring heatmap.

        Args:
            pathway_scores: Dictionary of pathway scores
            save_path: Path to save plot
            interactive: Use interactive plotly plot

        Returns:
            Plot object or None
        """
        if not pathway_scores:
            if self.verbose:
                print("⚠️  No pathway scores provided")
            return None

        if interactive and PLOTLY_AVAILABLE:
            return self._plot_pathway_heatmap_plotly(pathway_scores, save_path)
        elif MATPLOTLIB_AVAILABLE:
            return self._plot_pathway_heatmap_matplotlib(pathway_scores, save_path)
        else:
            if self.verbose:
                print("⚠️  No plotting libraries available")
            return None

    def _plot_pathway_heatmap_matplotlib(
        self,
        pathway_scores: Dict[str, float],
        save_path: Optional[str]
    ):
        """Create pathway heatmap with matplotlib."""

        # Prepare data
        pathways = list(pathway_scores.keys())
        scores = list(pathway_scores.values())

        # Create figure
        fig, ax = plt.subplots(figsize=(10, 6))

        # Create heatmap data
        heatmap_data = np.array(scores).reshape(1, -1)

        # Create heatmap
        im = ax.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto', vmin=0, vmax=1)

        # Set ticks and labels
        ax.set_xticks(range(len(pathways)))
        ax.set_xticklabels([p.replace('_', ' ').title() for p in pathways],
                          rotation=45, ha='right')
        ax.set_yticks([0])
        ax.set_yticklabels(['Pathway Score'])

        # Add score values
        for i, score in enumerate(scores):
            text_color = 'white' if score > 0.5 else 'black'
            ax.text(i, 0, f'{score:.2f}', ha="center", va="center",
                   color=text_color, fontweight='bold')

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Pathway Dysregulation Score', rotation=270, labelpad=20)

        ax.set_title('ReHMGB1 Pathway Scoring Profile', fontsize=14, fontweight='bold')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            if self.verbose:
                print(f"   💾 Pathway heatmap saved to: {save_path}")

        return fig

    def _plot_pathway_heatmap_plotly(
        self,
        pathway_scores: Dict[str, float],
        save_path: Optional[str]
    ):
        """Create interactive pathway heatmap with plotly."""

        pathways = list(pathway_scores.keys())
        scores = list(pathway_scores.values())

        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=[scores],
            x=[p.replace('_', ' ').title() for p in pathways],
            y=['Pathway Score'],
            colorscale='RdYlBu_r',
            zmin=0,
            zmax=1,
            text=[[f'{s:.2f}' for s in scores]],
            texttemplate="%{text}",
            textfont={"size": 12},
            colorbar=dict(title="Dysregulation Score")
        ))

        fig.update_layout(
            title="ReHMGB1 Pathway Scoring Profile",
            xaxis_title="Pathway",
            height=400,
            width=800
        )

        if save_path:
            fig.write_html(save_path)
            if self.verbose:
                print(f"   💾 Interactive pathway heatmap saved to: {save_path}")

        return fig

    def plot_clinical_risk_dashboard(
        self,
        aging_profile: Any,  # AgingProfile from clinical_scoring
        risk_scores: Dict[str, Any],  # ClinicalRiskScore objects
        save_path: Optional[str] = None,
        interactive: bool = False
    ) -> Optional[Any]:
        """
        Create comprehensive clinical risk dashboard.

        Args:
            aging_profile: AgingProfile object
            risk_scores: Dictionary of ClinicalRiskScore objects
            save_path: Path to save plot
            interactive: Use interactive plotly plot

        Returns:
            Plot object or None
        """
        if interactive and PLOTLY_AVAILABLE:
            return self._plot_risk_dashboard_plotly(aging_profile, risk_scores, save_path)
        elif MATPLOTLIB_AVAILABLE:
            return self._plot_risk_dashboard_matplotlib(aging_profile, risk_scores, save_path)
        else:
            if self.verbose:
                print("⚠️  No plotting libraries available")
            return None

    def _plot_risk_dashboard_matplotlib(
        self,
        aging_profile: Any,
        risk_scores: Dict[str, Any],
        save_path: Optional[str]
    ):
        """Create clinical risk dashboard with matplotlib."""

        # Create figure with subplots
        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

        # 1. Age comparison (top left)
        ax1 = fig.add_subplot(gs[0, 0])
        ages = ['Chronological', 'Biological']
        age_values = [aging_profile.chronological_age, aging_profile.biological_age]
        colors = ['#1f77b4', '#ff7f0e']

        bars = ax1.bar(ages, age_values, color=colors, alpha=0.7)
        ax1.set_title('Age Comparison')
        ax1.set_ylabel('Age (years)')

        # Add value labels
        for bar, value in zip(bars, age_values):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                    f'{value:.1f}', ha='center', va='bottom')

        # 2. Age acceleration gauge (top middle)
        ax2 = fig.add_subplot(gs[0, 1])
        self._create_gauge_plot(ax2, aging_profile.age_acceleration,
                               'Age Acceleration', 'years', (-15, 15))

        # 3. Senescence burden (top right)
        ax3 = fig.add_subplot(gs[0, 2])
        self._create_gauge_plot(ax3, aging_profile.senescence_burden * 100,
                               'Senescence Burden', '%', (0, 100))

        # 4. Risk scores heatmap (middle row, span all columns)
        ax4 = fig.add_subplot(gs[1, :])
        if risk_scores:
            risk_names = list(risk_scores.keys())
            risk_values = [rs.score for rs in risk_scores.values()]
            risk_categories = [rs.category for rs in risk_scores.values()]

            # Create heatmap
            heatmap_data = np.array(risk_values).reshape(1, -1)
            im = ax4.imshow(heatmap_data, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=1)

            # Set labels
            ax4.set_xticks(range(len(risk_names)))
            ax4.set_xticklabels([name.replace('_', ' ').title() for name in risk_names],
                               rotation=45, ha='right')
            ax4.set_yticks([0])
            ax4.set_yticklabels(['Risk Score'])

            # Add values and categories
            for i, (value, category) in enumerate(zip(risk_values, risk_categories)):
                text_color = 'white' if value > 0.5 else 'black'
                ax4.text(i, 0, f'{value:.1%}\\n({category})', ha="center", va="center",
                        color=text_color, fontweight='bold')

            ax4.set_title('Clinical Risk Assessment')

        # 5. Pathway scores (bottom left)
        ax5 = fig.add_subplot(gs[2, 0])
        if aging_profile.pathway_scores:
            pathway_names = list(aging_profile.pathway_scores.keys())
            pathway_values = list(aging_profile.pathway_scores.values())

            bars = ax5.barh(range(len(pathway_names)), pathway_values,
                           color=plt.cm.viridis(np.linspace(0, 1, len(pathway_names))))
            ax5.set_yticks(range(len(pathway_names)))
            ax5.set_yticklabels([name.replace('_', ' ').title() for name in pathway_names])
            ax5.set_xlabel('Pathway Score')
            ax5.set_title('Pathway Involvement')
            ax5.set_xlim(0, 1)

        # 6. Inflammation score (bottom middle)
        ax6 = fig.add_subplot(gs[2, 1])
        self._create_gauge_plot(ax6, aging_profile.inflammation_score * 100,
                               'Inflammation', '%', (0, 100))

        # 7. Summary text (bottom right)
        ax7 = fig.add_subplot(gs[2, 2])
        ax7.axis('off')

        summary_text = f"""SUMMARY

Biological Age: {aging_profile.biological_age:.1f} years
Age Acceleration: {aging_profile.age_acceleration:+.1f} years

Senescence: {aging_profile.senescence_burden:.1%}
Inflammation: {aging_profile.inflammation_score:.1%}

Overall Risk Level:
{self._get_overall_risk_level(risk_scores)}"""

        ax7.text(0.1, 0.9, summary_text, transform=ax7.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))

        plt.suptitle('Clinical Aging Assessment Dashboard', fontsize=16, fontweight='bold')

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            if self.verbose:
                print(f"   💾 Clinical dashboard saved to: {save_path}")

        return fig

    def _create_gauge_plot(self, ax, value, title, unit, value_range):
        """Create a simple gauge plot."""
        min_val, max_val = value_range

        # Normalize value to 0-1
        normalized = (value - min_val) / (max_val - min_val)
        normalized = max(0, min(1, normalized))

        # Create gauge
        theta = np.linspace(0, np.pi, 100)
        r = np.ones_like(theta)

        # Background arc
        ax.plot(r * np.cos(theta), r * np.sin(theta), 'lightgray', linewidth=8)

        # Value arc
        value_theta = theta[:int(normalized * len(theta))]
        if len(value_theta) > 0:
            color = 'red' if normalized > 0.7 else 'orange' if normalized > 0.4 else 'green'
            ax.plot(np.cos(value_theta), np.sin(value_theta), color, linewidth=8)

        # Add value text
        ax.text(0, -0.3, f'{value:.1f} {unit}', ha='center', va='center',
               fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=10)
        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-0.5, 1.2)
        ax.axis('off')

    def _get_overall_risk_level(self, risk_scores: Dict[str, Any]) -> str:
        """Determine overall risk level from individual risk scores."""
        if not risk_scores:
            return "Unknown"

        categories = [rs.category for rs in risk_scores.values()]

        if any(cat == 'very_high' for cat in categories):
            return "VERY HIGH"
        elif any(cat == 'high' for cat in categories):
            return "HIGH"
        elif any(cat == 'moderate' for cat in categories):
            return "MODERATE"
        else:
            return "LOW"

    def create_aging_report_visualizations(
        self,
        aging_data: Dict[str, Any],
        output_dir: str,
        prefix: str = 'aging_analysis'
    ) -> Dict[str, str]:
        """
        Create comprehensive set of aging analysis visualizations.

        Args:
            aging_data: Dictionary containing all aging analysis results
            output_dir: Output directory for plots
            prefix: File prefix for saved plots

        Returns:
            Dictionary mapping plot types to file paths
        """
        import os
        os.makedirs(output_dir, exist_ok=True)

        created_plots = {}

        # 1. Aging clock comparison
        if aging_data.get('aging_predictions') is not None:
            clock_plot_path = os.path.join(output_dir, f'{prefix}_clock_comparison.png')
            fig = self.plot_aging_clock_comparison(
                aging_data['aging_predictions'],
                aging_data.get('chronological_age'),
                save_path=clock_plot_path
            )
            if fig is not None:
                created_plots['clock_comparison'] = clock_plot_path
                plt.close(fig)

        # 2. ReHMGB1 pathway heatmap
        if aging_data.get('pathway_scores') is not None:
            pathway_plot_path = os.path.join(output_dir, f'{prefix}_pathway_heatmap.png')
            fig = self.plot_rehmgb1_pathway_heatmap(
                aging_data['pathway_scores'],
                save_path=pathway_plot_path
            )
            if fig is not None:
                created_plots['pathway_heatmap'] = pathway_plot_path
                plt.close(fig)

        # 3. Clinical risk dashboard
        if (aging_data.get('aging_profile') is not None and
            aging_data.get('risk_scores') is not None):
            dashboard_plot_path = os.path.join(output_dir, f'{prefix}_clinical_dashboard.png')
            fig = self.plot_clinical_risk_dashboard(
                aging_data['aging_profile'],
                aging_data['risk_scores'],
                save_path=dashboard_plot_path
            )
            if fig is not None:
                created_plots['clinical_dashboard'] = dashboard_plot_path
                plt.close(fig)

        if self.verbose:
            print(f"📊 Created {len(created_plots)} visualization(s)")
            for plot_type, path in created_plots.items():
                print(f"   {plot_type}: {path}")

        return created_plots


def create_demo_visualization() -> str:
    """Create demonstration script for aging biomarker visualization."""

    demo_script = '''#!/usr/bin/env python3
"""
Aging Biomarker Visualization Demo

Demonstrates comprehensive visualization capabilities for aging biomarker
analysis with focus on ReHMGB1/RAGE pathway visualization.
"""

import pandas as pd
import numpy as np
from crispr_toolkit.analysis.aging.biomarkers import AgingBiomarkerVisualizer

def main():
    print("📊 Aging Biomarker Visualization Demo")
    print("=" * 50)

    # Initialize visualizer
    print("\\n1. Initializing aging biomarker visualizer...")
    try:
        visualizer = AgingBiomarkerVisualizer(style='scientific', verbose=True)
    except Exception as e:
        print(f"   ❌ Visualization initialization failed: {e}")
        print("   📦 Install dependencies: pip install matplotlib seaborn plotly")
        return

    # Generate sample data
    print("\\n2. Generating sample aging biomarker data...")

    # Sample aging predictions
    np.random.seed(42)
    n_samples = 30
    sample_ids = [f"Sample_{i:03d}" for i in range(n_samples)]

    aging_predictions = []
    clocks = ['Horvath2013', 'PhenoAge', 'GrimAge', 'DunedinPACE']

    for sample_id in sample_ids:
        base_age = np.random.uniform(40, 75)
        for clock in clocks:
            if clock == 'DunedinPACE':
                # Pace measurement (around 1.0)
                pred_value = np.random.normal(1.0, 0.2)
                pred_age = pred_value
            else:
                # Age prediction with some variation
                pred_age = base_age + np.random.normal(0, 5)
                pred_value = pred_age

            aging_predictions.append({
                'sample_id': sample_id,
                'clock': clock,
                'predicted_age': pred_age,
                'predicted_value': pred_value
            })

    aging_df = pd.DataFrame(aging_predictions)

    # Chronological ages
    chronological_age = pd.Series(
        np.random.uniform(40, 75, n_samples),
        index=sample_ids
    )

    print(f"   📊 Generated data for {n_samples} samples and {len(clocks)} clocks")

    # Create aging clock comparison plot
    print("\\n3. Creating aging clock comparison visualization...")
    try:
        fig1 = visualizer.plot_aging_clock_comparison(
            aging_predictions=aging_df,
            chronological_age=chronological_age,
            save_path='./demo_clock_comparison.png'
        )
        print("   ✅ Clock comparison plot created")
    except Exception as e:
        print(f"   ⚠️  Clock comparison plot failed: {e}")

    # Create ReHMGB1 pathway heatmap
    print("\\n4. Creating ReHMGB1 pathway heatmap...")

    # Sample pathway scores
    pathway_scores = {
        'hmgb1_core': np.random.uniform(0.3, 0.8),
        'rage_signaling': np.random.uniform(0.4, 0.9),
        'jak_stat_pathway': np.random.uniform(0.2, 0.7),
        'nfkb_signaling': np.random.uniform(0.3, 0.8),
        'sasp_factors': np.random.uniform(0.5, 0.9),
        'oxidative_stress': np.random.uniform(0.3, 0.6)
    }

    try:
        fig2 = visualizer.plot_rehmgb1_pathway_heatmap(
            pathway_scores=pathway_scores,
            save_path='./demo_pathway_heatmap.png'
        )
        print("   ✅ Pathway heatmap created")
        print("   🧬 Sample pathway scores:")
        for pathway, score in pathway_scores.items():
            print(f"      {pathway}: {score:.2f}")
    except Exception as e:
        print(f"   ⚠️  Pathway heatmap failed: {e}")

    # Create mock aging profile and risk scores for dashboard
    print("\\n5. Creating clinical risk dashboard...")

    # Mock aging profile (would normally come from ClinicalAgingScorer)
    class MockAgingProfile:
        def __init__(self):
            self.chronological_age = 55.0
            self.biological_age = 62.3
            self.age_acceleration = 7.3
            self.mortality_risk = 65.8
            self.senescence_burden = 0.68
            self.inflammation_score = 0.72
            self.pathway_scores = pathway_scores

    # Mock risk scores
    class MockRiskScore:
        def __init__(self, score, category):
            self.score = score
            self.category = category

    aging_profile = MockAgingProfile()
    risk_scores = {
        'cardiovascular_risk': MockRiskScore(0.23, 'moderate'),
        'all_cause_mortality': MockRiskScore(0.31, 'high'),
        'rehmgb1_related_aging': MockRiskScore(0.67, 'high')
    }

    try:
        fig3 = visualizer.plot_clinical_risk_dashboard(
            aging_profile=aging_profile,
            risk_scores=risk_scores,
            save_path='./demo_clinical_dashboard.png'
        )
        print("   ✅ Clinical dashboard created")
        print(f"   📈 Patient profile:")
        print(f"      Chronological age: {aging_profile.chronological_age:.1f} years")
        print(f"      Biological age: {aging_profile.biological_age:.1f} years")
        print(f"      Age acceleration: {aging_profile.age_acceleration:+.1f} years")
        print(f"      Senescence burden: {aging_profile.senescence_burden:.1%}")
    except Exception as e:
        print(f"   ⚠️  Clinical dashboard failed: {e}")

    # Create interactive plots if plotly available
    print("\\n6. Creating interactive visualizations...")
    try:
        # Interactive clock comparison
        fig4 = visualizer.plot_aging_clock_comparison(
            aging_predictions=aging_df,
            chronological_age=chronological_age,
            save_path='./demo_clock_comparison_interactive.html',
            interactive=True
        )

        # Interactive pathway heatmap
        fig5 = visualizer.plot_rehmgb1_pathway_heatmap(
            pathway_scores=pathway_scores,
            save_path='./demo_pathway_heatmap_interactive.html',
            interactive=True
        )

        print("   ✅ Interactive plots created")
        print("   🌐 Open .html files in browser for interactive features")

    except Exception as e:
        print(f"   ⚠️  Interactive plots failed: {e}")
        print("   💡 Install plotly for interactive visualizations: pip install plotly")

    print("\\n🎉 Aging biomarker visualization demo completed!")
    print("\\n📁 Generated files:")
    print("   - demo_clock_comparison.png")
    print("   - demo_pathway_heatmap.png")
    print("   - demo_clinical_dashboard.png")
    print("   - demo_clock_comparison_interactive.html (if plotly available)")
    print("   - demo_pathway_heatmap_interactive.html (if plotly available)")

    print("\\n💡 Visualization Applications:")
    print("   - Research presentations and publications")
    print("   - Clinical report generation")
    print("   - Patient communication tools")
    print("   - Longitudinal aging monitoring")
    print("   - ReHMGB1 pathway analysis")

if __name__ == "__main__":
    main()
'''

    return demo_script
    return demo_script
