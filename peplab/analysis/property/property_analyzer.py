from typing import Dict, Optional
import pandas as pd


class PropertyAnalysis:
    """Statistical analysis of peptide library properties."""

    @staticmethod
    def calculate_statistics(df: pd.DataFrame) -> Dict[str, Dict[str, float]]:
        """Calculates statistical measures for numerical properties.

        Args:
            df: DataFrame containing analysis results.

        Returns:
            Dictionary of property statistics.
        """
        numerical_columns = df.select_dtypes(include=['float64', 'int64']).columns
        stats = {}

        for col in numerical_columns:
            stats[col] = {
                'mean': df[col].mean(),
                'std': df[col].std(),
                'min': df[col].min(),
                'max': df[col].max(),
                'median': df[col].median()
            }

        return stats

    @staticmethod
    def generate_summary(
        df: pd.DataFrame,
        output_path: Optional[str] = None
    ) -> str:
        """Generates a text summary of the library analysis.

        Args:
            df: DataFrame containing analysis results.
            output_path: Optional path to save summary text.

        Returns:
            Summary text.
        """
        stats = PropertyAnalysis.calculate_statistics(df)

        summary_lines = [
            "Peptide Library Analysis Summary",
            "==============================",
            f"Total Peptides: {len(df)}",
            f"Linear Peptides: {sum(df['type'] == 'linear')}",
            f"Cyclic Peptides: {sum(df['type'] == 'cyclic')}",
            "\nProperty Statistics:",
            "-------------------"
        ]

        for prop, prop_stats in stats.items():
            summary_lines.extend([
                f"\n{prop}:",
                f"  Mean: {prop_stats['mean']:.2f}",
                f"  Std Dev: {prop_stats['std']:.2f}",
                f"  Min: {prop_stats['min']:.2f}",
                f"  Max: {prop_stats['max']:.2f}",
                f"  Median: {prop_stats['median']:.2f}"
            ])

        summary = "\n".join(summary_lines)

        if output_path:
            with open(output_path, 'w') as f:
                f.write(summary)

        return summary
