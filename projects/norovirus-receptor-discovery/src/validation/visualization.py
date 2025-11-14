"""
Visualization Tools

Create plots and visualizations for pipeline results.
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Optional
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10


def plot_score_distribution(
    ranked_file: Path,
    output_file: Path,
    metric: str = "overall_score"
):
    """
    Plot distribution of scores

    Args:
        ranked_file: Path to ranked_candidates.json
        output_file: Where to save plot
        metric: Score metric to plot
    """
    with open(ranked_file) as f:
        candidates = json.load(f)

    scores = [c[metric] for c in candidates]

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(scores, bins=30, edgecolor='black', alpha=0.7)
    ax.axvline(np.median(scores), color='red', linestyle='--',
               label=f'Median: {np.median(scores):.3f}')
    ax.axvline(np.mean(scores), color='blue', linestyle='--',
               label=f'Mean: {np.mean(scores):.3f}')

    ax.set_xlabel(metric.replace('_', ' ').title())
    ax.set_ylabel('Count')
    ax.set_title(f'Distribution of {metric.replace("_", " ").title()}')
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved plot to {output_file}")


def plot_score_comparison(
    ranked_file: Path,
    output_file: Path,
    top_n: int = 50
):
    """
    Compare different scoring components for top candidates

    Args:
        ranked_file: Path to ranked_candidates.json
        output_file: Where to save plot
        top_n: Number of top candidates to show
    """
    with open(ranked_file) as f:
        candidates = json.load(f)

    top_candidates = candidates[:top_n]

    # Extract scores
    genes = [c['gene_name'] for c in top_candidates]
    ipTM = [c['ipTM_score'] for c in top_candidates]
    structural = [c['structural_confidence'] for c in top_candidates]
    biological = [c['biological_relevance'] for c in top_candidates]
    overall = [c['overall_score'] for c in top_candidates]

    # Create plot
    x = np.arange(len(genes))
    width = 0.2

    fig, ax = plt.subplots(figsize=(14, 8))

    ax.bar(x - 1.5*width, ipTM, width, label='ipTM', alpha=0.8)
    ax.bar(x - 0.5*width, structural, width, label='Structural', alpha=0.8)
    ax.bar(x + 0.5*width, biological, width, label='Biological', alpha=0.8)
    ax.bar(x + 1.5*width, overall, width, label='Overall', alpha=0.8)

    ax.set_xlabel('Candidate')
    ax.set_ylabel('Score')
    ax.set_title(f'Score Comparison for Top {top_n} Candidates')
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=90, ha='right')
    ax.legend()
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved plot to {output_file}")


def plot_confidence_tiers(
    ranked_file: Path,
    output_file: Path
):
    """
    Plot distribution of confidence tiers

    Args:
        ranked_file: Path to ranked_candidates.json
        output_file: Where to save plot
    """
    with open(ranked_file) as f:
        candidates = json.load(f)

    # Count tiers
    tiers = {'high': 0, 'medium': 0, 'low': 0}
    for c in candidates:
        tier = c.get('confidence_tier', 'low')
        tiers[tier] = tiers.get(tier, 0) + 1

    # Create pie chart
    fig, ax = plt.subplots(figsize=(8, 8))

    colors = {'high': '#2ecc71', 'medium': '#f39c12', 'low': '#e74c3c'}
    labels = [f"{tier.title()}\n({count})" for tier, count in tiers.items()]
    values = list(tiers.values())
    colors_list = [colors[tier] for tier in tiers.keys()]

    ax.pie(values, labels=labels, colors=colors_list, autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 12, 'weight': 'bold'})
    ax.set_title('Confidence Tier Distribution', fontsize=14, weight='bold')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved plot to {output_file}")


def plot_ipTM_vs_expression(
    ranked_file: Path,
    library_file: Path,
    output_file: Path
):
    """
    Scatter plot of ipTM vs intestinal expression

    Args:
        ranked_file: Path to ranked_candidates.json
        library_file: Path to candidate_library.json
        output_file: Where to save plot
    """
    with open(ranked_file) as f:
        candidates = json.load(f)

    with open(library_file) as f:
        library = json.load(f)

    # Create mapping
    expression_map = {
        c['uniprot_id']: c.get('intestinal_expression', 0)
        for c in library
    }

    # Extract data
    ipTM_scores = []
    expressions = []
    gene_names = []
    colors = []

    color_map = {'high': 'green', 'medium': 'orange', 'low': 'red'}

    for c in candidates:
        protein_id = c['protein_id']
        if protein_id in expression_map:
            ipTM_scores.append(c['ipTM_score'])
            expressions.append(expression_map[protein_id])
            gene_names.append(c['gene_name'])
            colors.append(color_map[c['confidence_tier']])

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(10, 8))

    scatter = ax.scatter(expressions, ipTM_scores, c=colors, alpha=0.6, s=100)

    # Annotate top candidates
    top_candidates = sorted(
        zip(expressions, ipTM_scores, gene_names),
        key=lambda x: x[1],
        reverse=True
    )[:10]

    for expr, iptm, gene in top_candidates:
        ax.annotate(gene, (expr, iptm), fontsize=8,
                   xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel('Intestinal Expression (TPM)', fontsize=12)
    ax.set_ylabel('ipTM Score', fontsize=12)
    ax.set_title('Binding Confidence vs Expression', fontsize=14, weight='bold')
    ax.grid(alpha=0.3)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='green', label='High Confidence'),
        Patch(facecolor='orange', label='Medium Confidence'),
        Patch(facecolor='red', label='Low Confidence')
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    logger.info(f"Saved plot to {output_file}")


def generate_all_plots(
    ranked_file: Path,
    library_file: Path,
    output_dir: Path
):
    """
    Generate all standard visualization plots

    Args:
        ranked_file: Path to ranked_candidates.json
        library_file: Path to candidate_library.json
        output_dir: Directory to save plots
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Generating visualization plots...")

    # Score distribution
    plot_score_distribution(
        ranked_file,
        output_dir / "score_distribution.png",
        "overall_score"
    )

    # Score comparison
    plot_score_comparison(
        ranked_file,
        output_dir / "score_comparison.png",
        top_n=30
    )

    # Confidence tiers
    plot_confidence_tiers(
        ranked_file,
        output_dir / "confidence_tiers.png"
    )

    # ipTM vs expression (if expression data available)
    try:
        plot_ipTM_vs_expression(
            ranked_file,
            library_file,
            output_dir / "iptm_vs_expression.png"
        )
    except Exception as e:
        logger.warning(f"Could not generate ipTM vs expression plot: {e}")

    logger.info(f"✓ All plots saved to {output_dir}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate visualization plots for pipeline results"
    )
    parser.add_argument(
        "--ranked-candidates",
        type=Path,
        required=True,
        help="Path to ranked_candidates.json"
    )
    parser.add_argument(
        "--candidate-library",
        type=Path,
        required=True,
        help="Path to candidate_library.json"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default="results/plots",
        help="Output directory for plots"
    )

    args = parser.parse_args()

    generate_all_plots(
        args.ranked_candidates,
        args.candidate_library,
        args.output_dir
    )

    print(f"\n✓ Plots generated in {args.output_dir}")


if __name__ == "__main__":
    main()
