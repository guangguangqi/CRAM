# Production-Grade Whole Genome Sequencing (WGS) Coverage & QC Pipeline

A containerized, memory-efficient bioinformatics pipeline designed to process Whole Genome Sequencing (WGS) data in CRAM format. The ecosystem performs rapid cross-validation metrics checking against standard reference assemblies, filters out structural decoy/random contig noise, executes data auditing to prevent pipeline false-negatives, and renders publication-ready high-resolution visualizations.

## Core Engineering Features
- **Memory-Optimized Architecture**: Employs streaming UNIX memory pipes (`subprocess.PIPE`) to stream metrics dynamically from `samtools` directly into `Pandas` frames, completely bypassing intermediate file disk I/O and allowing a 20GB dataset to run safely on machines with only 8GB of RAM.
- **Adaptive Parameter Polymorphism**: Dynamically audits system versions (supporting legacy fallback routines back to `samtools` v1.12) to ensure structural stability across varying high-performance computing environments.
- **Defensive Error Handling**: Incorporates explicit genomic coordinate validation and boundary intersection checks (with an adaptive `1e-4` float-precision tolerance matrix) to intercept assembly build drifts and corrupt or empty alignment streams.

---

## Getting Started: Quick Start Execution

Follow this exact 4-step sequence to build the decoupled analysis environment, download the tracking reference data, run the execution loop, and run the validation test suite.

### Step 1: Build the Docker Environment
Compile the isolated container environment layer utilizing the reproducible Micromamba package framework.
```bash
chmod +x ./dockin
./dockin
```

### Step 2: Download the Datasets and Reference Indices
Execute the automated dataset capture utility. This pulls down the heavy dataset components from cloud storage directly to your host disk and dynamically fires up a micro-container to execute reference genome indexing (`samtools faidx`) instantly.
```bash
chmod +x ./download_data.sh
./download_data.sh
```

### Step 3: Execute the Core Analytics Pipeline
Launch the production analysis orchestration framework. This mounts your local data workspace into the container, analyzes depth matrices, isolates core human chromosomes (`chr1-22, X, Y`), and saves the polished deliverables directly to your host machine's output partition.
```bash
chmod +x ./dockrun
./dockrun
```

### Step 4: Run the Automated Pytest Verification Suite
Verify the structural integrity of your output metrics using the decoupled automated validation harness. This runs a series of strict Pytest data assertions and exports a telemetry score ledger directly to your local file logs.
```bash
chmod +x ./docktest
./docktest
```

---

## File Workspace Layout

Following successful completion of the lifecycle execution loop, your workspace folders will structurally align as follows on your local host system:

```text
my_project_root/
├── .dockerignore                 # Prevents 20GB data caching on build
├── .gitignore                    # Ensures data/ and output/ bypass Git tracking
├── dockin                        # Image build script
├── download_data.sh              # Cloud-to-host data provisioning script
├── dockrun                       # Master pipeline execution runner
├── docktest                      # Automated testing validation wrapper
├── environment/
│   └── Dockerfile                # Multi-channel dependency definition layer
├── solution/
│   └── cram_coverage.py          # Unified calculation, filtering, and plotting code
├── tests/
│   ├── test.sh                   # Pytest execution bash harness
│   └── test_outputs.py           # Explicit dataset data asset assertions
├── data/                         # Local dataset sandbox (Volume mapped at runtime)
│   ├── COLO829T_TEST.cram        
│   ├── COLO829T_TEST.cram.crai   
│   ├── GCA_000001405...fa        
│   └── GCA_000001405...fa.fai    # Automatically generated indexing file
└── output/                       # Mapped analysis deliverables directory
    ├── final_qc_report.tsv       # Baseline metrics matrix across all sequences
    ├── final_qc_report_filtered.tsv # Cleansed core human genome dataset
    └── wgs_chromosome_coverage.png  # High-resolution (300 DPI) publication-grade plot
```

---

## Project Dependencies
The underlying system relies on the following structural binaries and packages isolated entirely within the container universe:
- `samtools=1.22`
- `python=3.11`
- `pandas=2.2.2`
- `numpy=1.26.4`
- `seaborn=0.13.2`
- `matplotlib=3.8.4`
- `pytest=8.2.2`

