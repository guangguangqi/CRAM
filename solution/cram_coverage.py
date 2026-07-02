import argparse
import os
import subprocess
import sys
from io import StringIO
import numpy as np
import pandas as pd
import matplotlib
# Force headless file-writing mode so the script runs smoothly inside Docker
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import seaborn as sns


def check_samtools_version() -> float:
    """
    Audits the system's samtools installation.
    Prints the software signature and extracts a clean float version number 
    for dynamic parameter switching.
    """
    try:
        out = subprocess.check_output(
            ["samtools", "--version"], 
            text=True, 
            stderr=subprocess.DEVNULL
        )
        lines = out.splitlines()
        if not lines:
            raise ValueError("Empty output received from samtools --version.")
            
        version_line = lines[0].strip()
        print(f"[INFO] {version_line}", file=sys.stderr)
        
        tokens = version_line.split()
        for token in tokens:
            if any(char.isdigit() for char in token):
                clean_token = token.split("-")[0]
                version_parts = clean_token.split(".")[:2]
                if len(version_parts) >= 2:
                    return float(".".join(version_parts))
                    
        return 1.12

    except subprocess.CalledProcessError as e:
        print(f"[Pipeline Error] samtools execution binary found but failed to return version status: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("[Pipeline Error] samtools binary not found in the local PATH environment.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[WARNING] Version audit encountered an unexpected anomaly: {e}. Defaulting to legacy mode.", file=sys.stderr)
        return 1.12


def execute_samtools_coverage(cram_path: str, ref_path: str) -> pd.DataFrame:
    """
    Executes 'samtools coverage' via subprocess.
    Applies polymorphic arrays based on version compliance to handle CRAM streaming.
    """
    version_num = check_samtools_version()

    # Apply adaptive parameter polymorphism (FIXED: Removed the accidental override below this block)
    if version_num >= 1.22:
        print(f"[INFO] Configuring modern command array for samtools v{version_num}.", file=sys.stderr)
        cmd = ["samtools", "coverage", "--reference", ref_path, cram_path]
        env_config = os.environ.copy()
    else:
        print(f"[INFO] Configuring environment fallback array for legacy samtools v{version_num}.", file=sys.stderr)
        cmd = ["samtools", "coverage", cram_path]
        env_config = os.environ.copy()
        env_config["REF_PATH"] = f"{ref_path}:%s"

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env_config,
            text=True
        )
        stdout, stderr = proc.communicate()

        if stderr and proc.returncode == 0:
            print(f"[SAMTOOLS WARNINGS]\n{stderr}", file=sys.stderr)

        if proc.returncode != 0:
            raise RuntimeError(f"samtools binary execution failed (exit code {proc.returncode}):\n{stderr}")

        df = pd.read_csv(StringIO(stdout), sep="\t")
        df.columns = [c.lower().replace("#", "").strip() for c in df.columns]
        
        required_cols = {"rname", "covbases"}
        if not required_cols.issubset(df.columns):
            raise KeyError(f"Expected metrics missing from samtools output: {list(df.columns)}")
            
        return df

    except Exception as e:
        print(f"[Pipeline Error] Process execution block failed: {e}", file=sys.stderr)
        sys.exit(1)


def load_assembly_lengths(ref_fai_path: str) -> dict:
    """Parses reference .fai index to load chromosome structural boundaries."""
    chrom_lengths = {}
    try:
        with open(ref_fai_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    chrom_lengths[parts[0]] = int(parts[1])
    except Exception as e:
        print(f"[Pipeline Error] Failed to read FASTA index (.fai): {e}", file=sys.stderr)
        sys.exit(1)
    return chrom_lengths


def generate_standardized_report(df: pd.DataFrame, chrom_lengths: dict) -> pd.DataFrame:
    """
    Applies reference assembly structural boundaries to cross-validate calculations.
    Includes verification checks against text-truncated samtools precision limits.
    """
    df = df.copy()
    
    # Verify contig naming schemas (Ensembl vs UCSC build patterns)
    sample_chroms = set(df[df["numreads"] > 0]["rname"].unique()) if "numreads" in df.columns else set(df["rname"].unique())
    index_chroms = set(chrom_lengths.keys())
    intersection = sample_chroms.intersection(index_chroms)
    
    if len(sample_chroms) > 0 and len(intersection) / len(sample_chroms) < 0.70:
        print("[WARNING] Low contig overlap ratio. Check for assembly build drift.", file=sys.stderr)

    # Biological QC Sanity Check: Identify empty or completely unmapped alignment streams
    read_col = [col for col in ["numreads", "num_reads"] if col in df.columns]
    if read_col and df[read_col[0]].sum() == 0:
        print("[CRITICAL BIOLOGICAL WARNING] Total mapped reads across all chromosomes is ZERO.", file=sys.stderr)

    # Align metrics to strict physical boundaries from reference index file
    df["chrom_length"] = df["rname"].map(chrom_lengths)
    if df["chrom_length"].isna().any():
        df = df.dropna(subset=["chrom_length"])
    df["chrom_length"] = df["chrom_length"].astype(int)
    
    # Safely preserve both tracking fields instead of silent semantic overwrites
    if "coverage" in df.columns:
        df = df.rename(columns={"coverage": "samtools_raw_coverage_pct"})
    
    # Calculate independent coverage fraction
    df["computed_coverage_fraction_pct"] = (df["covbases"] / df["chrom_length"]) * 100
    
    # Cross-validate calculations vs samtools tracking array
    if "samtools_raw_coverage_pct" in df.columns:
        coverage_difference = np.abs(df["computed_coverage_fraction_pct"] - df["samtools_raw_coverage_pct"])
        max_drift = coverage_difference.max()
        
        # Relaxed to 1e-4 tolerance to bypass text string truncation variances safely
        if not (coverage_difference < 1e-4).all():
            print(f"[WARNING] Calculation variance detected vs samtools. Max drift: {max_drift}", file=sys.stderr)
        else:
            print(f"[INFO] Cross-validation successful. Matches samtools within safe rounding bounds. Max drift: {max_drift}", file=sys.stderr)
            
    return df  # FIXED: Explicit return statement added



def filter_tsv_report(input_file: str, output_file: str):
    """
    Reads a raw samtools/WGS coverage TSV file and filters out 
    all unplaced, decoy, and random contigs using explicit regex matching.
    """
    if not os.path.exists(input_file):
        print(f"[ERROR] Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)
        
    print(f"[INFO] Reading raw report: {input_file}", file=sys.stderr)
    df = pd.read_csv(input_file, sep="\t")
    
    # Capture structural size before filtering
    total_raw_rows = len(df)
    
    # Regex Pattern: Matches exactly chr1-chr22, chrX, chrY, and chrM
    # Prevents items like 'chr1_random' or 'chrUn' from passing through
    core_chrom_pattern = r"^chr([1-9]|1[0-9]|2[0-2]|[XYM])$"
    
    # Apply the string match filter
    filtered_df = df[df["rname"].astype(str).str.match(core_chrom_pattern, na=False)].copy()
    
    if filtered_df.empty:
        print("[CRITICAL ERROR] Filtering resulted in an empty dataset!", file=sys.stderr)
        print("Verify your file uses the 'chr' naming prefix (e.g., 'chr1' vs '1').", file=sys.stderr)
        sys.exit(1)
        
    # Optional: Ensure the rows follow standard biological order (1 to 22, X, Y, M)
    # This keeps your final output table clean and beautiful
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    filtered_df["rname"] = pd.Categorical(filtered_df["rname"], categories=chrom_order, ordered=True)
    filtered_df = filtered_df.sort_values("rname").dropna(subset=["rname"])
    
    # Save the polished output table
    filtered_df.to_csv(output_file, sep="\t", index=False)
    
    # Print pipeline execution metrics
    removed_count = total_raw_rows - len(filtered_df)
    print(f"\n[SUCCESS] Filtering operation complete.", file=sys.stderr)
    print(f" -> Total raw sequences processed: {total_raw_rows}", file=sys.stderr)
    print(f" -> Decoy/Random contigs removed:  {removed_count}", file=sys.stderr)
    print(f" -> Core primary chromosomes kept: {len(filtered_df)}", file=sys.stderr)
    print(f" -> Clean report destination:      {output_file}\n", file=sys.stderr)
    
    # Print quick biological stats of the core genome
    avg_depth = filtered_df[filtered_df["rname"] != "chrM"]["meandepth"].mean()
    avg_cov = filtered_df[filtered_df["rname"] != "chrM"]["samtools_raw_coverage_pct"].mean()
    print(f"[QC METRICS] Core Autosome Average Depth:   {avg_depth:.2f}x", file=sys.stderr)
    print(f"[QC METRICS] Core Autosome Average Breadth: {avg_cov:.2f}%", file=sys.stderr)



def generate_publication_plot(data_path: str, output_image_path: str):
    """
    Reads the filtered TSV report and generates a publication-quality 
    horizontal bar plot displaying chromosome breadth of coverage.
    """
    if not os.path.exists(data_path):
        print(f"[ERROR] Input file '{data_path}' not found. Please run the filter script first.", file=sys.stderr)
        sys.exit(1)

    # 1. Load dataset and exclude mitochondrial genome (chrM is typically ~100% and 
    # analyzed separately from the core nuclear genome for clarity)
    df = pd.read_csv(data_path, sep="\t")
    df = df[df["rname"] != "chrM"].copy()

    # 2. Configure publication-grade styling and high-resolution parameters
    # Apply a clean, professional sans-serif typography layout
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    plt.rcParams['text.color'] = '#333333'
    plt.rcParams['axes.labelcolor'] = '#333333'
    plt.rcParams['xtick.color'] = '#333333'
    plt.rcParams['ytick.color'] = '#333333'

    # Initialize canvas dimensions fitting standard single/double column journal space
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    
    # 3. Render core horizontal bar array using a muted academic deep blue palette
    sns.barplot(
        x="samtools_raw_coverage_pct",
        y="rname",
        data=df,
        ax=ax,
        color="#2b5c8f",     # Classic journal deep blue
        edgecolor="#1a3a5f", # Defined borders to elevate visual structure
        linewidth=0.6
    )

    # 4. Polish axis margins and desaturate outer structural frames
    ax.spines['top'].set_visible(False)   
    ax.spines['right'].set_visible(False) 
    ax.spines['left'].set_linewidth(0.8)  
    ax.spines['bottom'].set_linewidth(0.8)

    # 5. Optimize alignment grids and limits
    ax.set_xlim(0, 105) # Provide safe padding layout space on the right for text metrics
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    ax.xaxis.grid(True, linestyle='--', alpha=0.4, color='#cccccc', zorder=0) 
    ax.set_axisbelow(True) # Force horizontal gridlines to sit underneath the bars

    # Apply bold typography attributes to axis titles
    ax.set_xlabel("Breadth of Coverage (%)", fontsize=11, fontweight='bold', labelpad=10)
    ax.set_ylabel("Chromosome", fontsize=11, fontweight='bold', labelpad=10)
    ax.tick_params(axis='both', labelsize=9)

    # 6. Dynamic Text Annotations: Stamp raw numerical values to the right of each element
    for p in ax.patches:
        width = p.get_width()
        if width > 0: 
            ax.text(
                width + 1.0,                    # Horizontal coordinate offset
                p.get_y() + p.get_height()/2,   # Centered vertically inside the bar frame
                f"{width:.2f}%",                # Precision restricted to two decimals
                ha='left', 
                va='center', 
                fontsize=8, 
                color='#444444',
                fontweight='normal'
            )

    # 7. Compress margins uniformly and export to drive
    plt.tight_layout()
    plt.savefig(output_image_path, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"[SUCCESS] Publication-quality visualization rendered.", file=sys.stderr)
    print(f" -> Image destination: {output_image_path}", file=sys.stderr)


def main():
    # 1. Retrieve standardized sandboxed environment variables
    default_data_dir = os.getenv("DATA_DIR", "/workspace/data")
    default_out_dir = os.getenv("OUTPUT_DIR", "/workspace/output")

    # 2. Parse command-line arguments with explicit fallback paths
    parser = argparse.ArgumentParser(description="Dockerized WGS CRAM Coverage & Visualization Gate")
    parser.add_argument("--cram", default=os.path.join(default_data_dir, "COLO829T_TEST.cram"))
    parser.add_argument("--ref", default=os.path.join(default_data_dir, "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"))
    parser.add_argument("--fai", default=os.path.join(default_data_dir, "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai"))
    
    # Primary output targets routed directly to the mounted directory
    parser.add_argument("--out_raw", default=os.path.join(default_out_dir, "final_qc_report.tsv"))
    parser.add_argument("--out_filtered", default=os.path.join(default_out_dir, "final_qc_report_filtered.tsv"))
    parser.add_argument("--out_png", default=os.path.join(default_out_dir, "wgs_chromosome_coverage.png"))
    args = parser.parse_args()

    # Ensure output directory exists before writing data
    os.makedirs(default_out_dir, exist_ok=True)

    # =========================================================================
    # PHASE 1: Execution & Independent Cross-Validation
    # =========================================================================
    print("[PIPELINE STEP 1/3] Loading genomic indexes and computing raw coverage...", file=sys.stderr)
    chrom_lengths = load_assembly_lengths(args.fai)
    raw_coverage_df = execute_samtools_coverage(args.cram, args.ref)
    final_report_df = generate_standardized_report(raw_coverage_df, chrom_lengths)
    
    # Save the initial raw baseline data matrix
    cols_to_write = ["rname", "chrom_length"] + [col for col in final_report_df.columns if col not in ["rname", "chrom_length"]]
    final_report_df[cols_to_write].to_csv(args.out_raw, sep="\t", index=False)
    print(f" -> Raw report saved to: {args.out_raw}", file=sys.stderr)

    # =========================================================================
    # PHASE 2: Core Human Chromosome Regex Filtering
    # =========================================================================
    print("\n[PIPELINE STEP 2/3] Filtering out decoy/random contig noise...", file=sys.stderr)
    filter_tsv_report(args.out_raw, args.out_filtered)

    # =========================================================================
    # PHASE 3: Publication-Quality Rendering
    # =========================================================================
    print("\n[PIPELINE STEP 3/3] Generating high-resolution publication plot...", file=sys.stderr)
    generate_publication_plot(args.out_filtered, args.out_png)
    
    print("\n[SUCCESS] Whole pipeline execution loop complete!", file=sys.stderr)



if __name__ == "__main__":
    main()

