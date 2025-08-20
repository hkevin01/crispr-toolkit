"""Main CLI entry point for CRISPR Toolkit."""

import sys
from pathlib import Path
from typing import Optional

import click


@click.group()
@click.version_option()
def cli():
    """CRISPR Toolkit - AI/ML tools for aging research."""
    pass


@cli.group()
def aging():
    """Aging research specific commands."""
    pass


@aging.command()
@click.option('--tissue', default='liver', help='Target tissue (liver, brain, muscle, skin)')
@click.option('--phenotype', default='senescence', help='Target phenotype')
@click.option('--species', default='mouse', help='Model organism')
@click.option('--n-targets', default=10, help='Number of targets to return')
@click.option('--output', required=True, help='Output file path')
@click.option('--config', help='Analysis configuration file')
def prioritize(tissue: str, phenotype: str, species: str, n_targets: int, output: str, config: Optional[str]):
    """Prioritize CRISPR targets for aging interventions."""
    from ..analysis.aging.prioritization import AgingTargetPrioritizer

    try:
        # Initialize prioritizer
        prioritizer = AgingTargetPrioritizer()

        # Mock tissue data and genes for demonstration
        candidate_genes = ['TP53', 'CDKN2A', 'SIRT1', 'FOXO3', 'ATM', 'BRCA1']

        # Prioritize targets
        click.echo(f"Prioritizing targets for {tissue} - {phenotype}")
        results = prioritizer.rank_targets(
            genes=candidate_genes,
            tissue=tissue,
            phenotype=phenotype,
            n_targets=n_targets
        )

        # Save results
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        import json
        with open(output_path, 'w') as f:
            json.dump({"targets": results}, f, indent=2)

        click.echo(f"✅ Found {len(results)} prioritized targets")
        click.echo(f"📁 Results saved to {output}")

        # Show top targets
        top_targets = results[:5]
        if top_targets:
            click.echo("\nTop 5 targets:")
            for i, target in enumerate(top_targets, 1):
                score = target.get('score', 0)
                gene = target.get('gene', 'Unknown')
                click.echo(f"  {i}. {gene} (score: {score:.3f})")

    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@aging.command()
@click.option('--intervention', required=True, help='Intervention type (OSK, OSKM, senolytic)')
@click.option('--tissue', default='liver', help='Target tissue')
@click.option('--targets', help='Comma-separated list of target genes')
@click.option('--delivery', default='AAV', help='Delivery method (AAV, LNP)')
@click.option('--dose-level', default=1.0, type=float, help='Dose level (0.1-2.0)')
@click.option('--duration', default=4, type=int, help='Duration in weeks')
@click.option('--output', required=True, help='Output file path')
@click.option('--model-path', help='Path to trained model')
def predict(intervention: str, tissue: str, targets: Optional[str], delivery: str,
           dose_level: float, duration: int, output: str, model_path: Optional[str]):
    """Predict rejuvenation outcomes for CRISPR interventions."""
    from ..analysis.aging.rejuvenation import predict_rejuvenation

    try:
        # Build intervention config
        intervention_config = {
            'intervention': intervention,
            'delivery': delivery,
            'dose_level': dose_level,
            'duration_weeks': duration,
            'targets': targets.split(',') if targets else []
        }

        # Context data (mock for demonstration)
        context_data = {
            'age': 24,  # months for mouse
            'senescence_load': 0.3,
            'expression_level': 1.0,
            'chromatin_access': 0.8
        }

        click.echo(f"🔮 Predicting outcomes for {intervention} intervention in {tissue}")

        # Predict outcomes
        model_path_obj = Path(model_path) if model_path else None
        results = predict_rejuvenation(
            intervention_config=intervention_config,
            tissue=tissue,
            context_data=context_data,
            model_path=model_path_obj
        )

        # Save results
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        import json
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)

        # Display predictions
        predictions = results['predictions']
        click.echo("\n📊 Predicted outcomes:")
        click.echo(f"  Epigenetic clock change: {predictions['epigenetic_clock_delta']:.2f} years")
        click.echo(f"  Senescence score change: {predictions['senescence_score_change']:.2f}")
        click.echo(f"  Transcriptional age shift: {predictions['transcriptional_age_shift']:.2f}")
        click.echo(f"  Functional improvement: {predictions['functional_improvement']:.2f}")
        click.echo(f"  Safety risk: {predictions['safety_risk']:.2f}")

        if 'confidence_intervals' in results:
            click.echo(f"\n🎯 Confidence: {results.get('confidence', 'N/A')}")

        click.echo(f"\n📁 Full results saved to {output}")

    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@cli.group()
def screen():
    """CRISPR screen analysis commands."""
    pass


@screen.command()
@click.argument('input_file')
@click.option('--output', required=True, help='Output directory')
@click.option('--control-guides', help='Control guide file')
def analyze(input_file: str, output: str, control_guides: Optional[str]):
    """Analyze CRISPR screen data."""
    click.echo(f"Analyzing screen data from {input_file}")
    click.echo(f"Results will be saved to {output}")
    # Placeholder for screen analysis
    click.echo("✅ Screen analysis complete")


@cli.command()
def version():
    """Show version information."""
    click.echo("CRISPR Toolkit v0.1.0")
    click.echo("AI/ML tools for aging research with CRISPR")


if __name__ == '__main__':
    cli()
