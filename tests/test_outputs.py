#!/usr/bin/env python3
"""
test_outputs.py

Verify that the CRAM coverage pipeline generated all expected outputs.

Checks:
1. Output files exist.
2. TSV files are readable.
3. Required columns exist.
4. Coverage values are valid.
5. PNG exists and is non-empty.
"""

import os
import sys
import pandas as pd

OUTPUT_DIR = os.getenv("OUTPUT_DIR", "/workspace/output")

REPORT = os.path.join(OUTPUT_DIR, "final_qc_report.tsv")
FILTERED = os.path.join(OUTPUT_DIR, "final_qc_report_filtered.tsv")
PLOT = os.path.join(OUTPUT_DIR, "wgs_chromosome_coverage.png")


def check_file_exists(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Missing output file: {filename}")

    if os.path.getsize(filename) == 0:
        raise RuntimeError(f"Empty file: {filename}")

    print(f"[PASS] {os.path.basename(filename)}")


def check_report(report):
    df = pd.read_csv(report, sep="\t")

    required = [
        "rname",
        "chrom_length",
        "numreads",
        "covbases",
        "meandepth",
        "computed_coverage_fraction_pct"
    ]

    missing = [c for c in required if c not in df.columns]

    if missing:
        raise RuntimeError(f"Missing columns: {missing}")

    if len(df) == 0:
        raise RuntimeError("Report contains zero rows.")

    if (df["chrom_length"] <= 0).any():
        raise RuntimeError("Invalid chromosome lengths detected.")

    if (df["covbases"] < 0).any():
        raise RuntimeError("Negative covbases detected.")

    if (df["computed_coverage_fraction_pct"] < 0).any():
        raise RuntimeError("Negative coverage detected.")

    if (df["computed_coverage_fraction_pct"] > 100).any():
        raise RuntimeError("Coverage exceeds 100%.")

    print(f"[PASS] {len(df)} chromosomes/contigs validated.")
    return df


def check_filtered(filtered):
    df = pd.read_csv(filtered, sep="\t")

    if len(df) == 0:
        raise RuntimeError("Filtered report is empty.")

    if "rname" not in df.columns:
        raise RuntimeError("Filtered report missing rname column.")

    print(f"[PASS] Filtered report contains {len(df)} chromosomes.")


def check_png(plot):
    size = os.path.getsize(plot)

    if size < 1000:
        raise RuntimeError("PNG appears corrupted or too small.")

    print(f"[PASS] PNG ({size:,} bytes)")


# ==============================================================================
# PYTEST ENTRY POINTS (Added to ensure Pytest discovers and runs your suite)
# ==============================================================================

def test_pipeline_files_exist():
    """Validates physical presence of pipeline assets."""
    check_file_exists(REPORT)
    check_file_exists(FILTERED)
    check_file_exists(PLOT)


def test_pipeline_raw_report():
    """Runs data matrices checks on raw metrics."""
    check_report(REPORT)


def test_pipeline_filtered_report():
    """Runs filtering assertions against final data tables."""
    check_filtered(FILTERED)


def test_pipeline_visualization_plot():
    """Ensures high-resolution visual layout compiled safely."""
    check_png(PLOT)


# ==============================================================================
# STANDALONE LOCAL RUNNER
# ==============================================================================
def main():
    print("=" * 60)
    print("Validating pipeline outputs")
    print("=" * 60)

    # Trigger components sequentially
    test_pipeline_files_exist()
    report = check_report(REPORT)
    test_pipeline_filtered_report()
    test_pipeline_visualization_plot()

    mapped = report[report["numreads"] > 0]

    print()
    print("Summary")
    print("-------")
    print(f"Mapped chromosomes : {len(mapped)}")
    print(f"Total chromosomes  : {len(report)}")

    if len(mapped) > 0:
        print("Top covered chromosome:")
        top = mapped.sort_values(
            "computed_coverage_fraction_pct",
            ascending=False
        ).iloc[0]

        print(
            f"  {top.rname}: "
            f"{top.computed_coverage_fraction_pct:.3f}% "
            f"({top.numreads:,} reads)"
        )

    print()
    print("[SUCCESS] All output validation tests passed.")


if __name__ == "__main__":
    main()

