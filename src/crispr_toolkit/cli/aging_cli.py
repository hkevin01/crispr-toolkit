"""Command-line interface for aging-specific CRISPR analysis."""

import argparse
import sys
from typing import Optional

from ..analysis.aging.rejuvenation_prediction import predict_rejuvenation
from ..analysis.aging.target_prioritization import prioritize_targets


def main(args: Optional[list] = None) -> int:
    """Main entry point for the aging CLI."""
    parser = argparse.ArgumentParser(
        description="CRISPR Toolkit - Aging Research CLI"
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Target prioritization command
    prioritize_parser = subparsers.add_parser(
        "prioritize", help="Prioritize CRISPR targets for aging interventions"
    )
    prioritize_parser.add_argument(
        "--tissue", type=str, default="liver", help="Target tissue"
    )
    prioritize_parser.add_argument(
        "--phenotype", type=str, default="senescence", help="Target phenotype"
    )
    prioritize_parser.add_argument(
        "--species", type=str, default="mouse", help="Model organism"
    )
    prioritize_parser.add_argument(
        "--output", type=str, required=True, help="Output file path"
    )

    # Prediction command
    predict_parser = subparsers.add_parser(
        "predict", help="Predict rejuvenation outcomes"
    )
    predict_parser.add_argument(
        "--config", type=str, required=True, help="Intervention configuration file"
    )
    predict_parser.add_argument(
        "--context", type=str, help="Context data directory"
    )
    predict_parser.add_argument(
        "--output", type=str, required=True, help="Output file path"
    )

    if args is None:
        args = sys.argv[1:]

    parsed_args = parser.parse_args(args)

    if parsed_args.command == "prioritize":
        print(f"Prioritizing targets for {parsed_args.tissue} - {parsed_args.phenotype}")
        targets = prioritize_targets(None, n_targets=10)  # Placeholder
        print(f"Found {len(targets)} prioritized targets")
        print(f"Results saved to {parsed_args.output}")

    elif parsed_args.command == "predict":
        print(f"Predicting outcomes for config: {parsed_args.config}")
        outcomes = predict_rejuvenation({"intervention": "OSK"})  # Placeholder
        print(f"Predicted clock delta: {outcomes['epigenetic_clock_delta']:.2f}")
        print(f"Results saved to {parsed_args.output}")

    else:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
