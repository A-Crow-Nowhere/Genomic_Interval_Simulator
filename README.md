# Genomic_Interval_Siumlator
Simulate CNV-like intervals at randomly generated loci across an arbitrary length genomic segment. The program generates a specified number of samples ('cells') to mimic technical and biological varation that is typical of bionformatic approaches to determine interval start and stop sites. 

## CNV simulator usage

This tool generates **canonical CNV regions** (truth set) and **per‑cell observed CNVs** with controllable positional/size variation. It supports **clonal** (>80% prevalence) and **subclonal** (rare) events, **run‑level variance**, and **high‑variance (HV) hotspot regions** that contain overlapping canonical CNVs.

All outputs are **BED‑style** and IGV‑ready.

---

## Installation (Option A: Bash wrapper)

No global installs required. The wrapper creates a local virtual environment and installs Python dependencies for you.

```bash
cd~
mkdir -p tools/cnvsim
cd tools/cnvsim
wget https://raw.githubusercontent.com/A-Crow-Nowhere/MalariAPI/main/scripts/cnv_sim/cnv_sim_cli_bundle.zip
unzip cnv_sim_cli_bundle.zip


See available flags (first run creates .venv and installs deps)
./run_cnv_sim.sh --help
```

The wrapper:
- creates/uses a local `.venv/` alongside the scripts,
- installs from `requirements.txt` (`numpy`, `pandas`),
- runs `cnv_sim.py` with your flags.

**Note:** You need a `chrom.sizes` file containing at least the chromosome you plan to simulate (an example `chr1.sizes` is included).

---

## Quick Start

Minimal run on hg38 `chr1`:

```bash
./run_cnv_sim.sh \
  --chrom-sizes chr1.sizes --chrom chr1 \
  --outdir cnv_run1
```

A more realistic run (100 cells, mixed clonal/subclonal, medium variance, 2 HV regions):

```bash
./run_cnv_sim.sh \
  --chrom-sizes chr1.sizes --chrom chr1 \
  --n-cells 100 \
  --n-clonal 15 \
  --n-subclonal 25 \
  --variance 50 \
  --hv-n-regions 2 \
  --hv-min-events 8 \
  --hv-radius 75000 \
  --seed 42 \
  --outdir cnv_run2
```

---

## Inputs

- **`--chrom-sizes <file>`** (required): Two‑column file with `chrom<TAB>length`. Must include `--chrom`.
  - Example (hg38 chr1):  
    ```text
    chr1    248956422
    ```

---

## Outputs (BED‑style, IGV‑ready)

Files are written to `--outdir` (default: `cnv_sim_run/`).

- **`canonical_events.bed`** — truth list  
  Columns: `chrom  start  end  name(event_id)  score(0)  strand(.)  true_freq  type  hv_regions`  
  - `true_freq`: target prevalence (e.g., `0.87` → ~87/100 cells).  
  - `type`: `clonal` or `subclonal`.  
  - `hv_regions`: comma‑separated HV IDs (e.g., `HV1,HV2`) or `.`.

- **`per_cell_cnvs.bed`** — all observed CNVs across all cells  
  Columns: `chrom  start  end  name  score(0)  strand(.)  true_freq  type  cell_id`  
  - `name` encodes: `event_id|truth:<freq>|var:<setting>`.

- **`per_cell/*.bed`** — one file per cell  
  Same columns as `per_cell_cnvs.bed`.

- **`event_observed_freqs.bed`** — counts & frequencies at truth coordinates  
  Columns: `chrom  start  end  event_id  score(0)  strand(.)  true_freq  type  observed_count  observed_freq`.

- **`high_variance_regions.bed`** — HV hotspots and their events  
  Columns: `chrom  start  end  name(HV_id)  score(0)  strand(.)  event_ids`.

**IGV tip:** Load `canonical_events.bed` as the truth track, add one or more `per_cell/*.bed` to visualize observed vs truth, and `high_variance_regions.bed` to highlight hotspots.

---

## CLI Flags and Empirical Effects

Run `./run_cnv_sim.sh --help` for current defaults. Below is an explanation of each flag and how it changes the simulation.

### Genome & Cohort Size
- `--chrom-sizes <file>` *(required)*: Chrom sizes table. Must include `--chrom`.
- `--chrom <name>` (default `chr1`): Chromosome to simulate.
- `--n-cells <int>` (default `100`): Number of cells.  
  **Effect:** Increases per‑cell rows; observed frequencies get tighter with larger `n` (binomial variance ↓).

### Canonical Truth Set
- `--n-clonal <int>` (default `15`): Number of **clonal** canonical events.  
  **Effect:** More high‑prevalence events → most cells will contain most of them.
- `--n-subclonal <int>` (default `25`): Number of **subclonal** canonical events.  
  **Effect:** More rare events → more sporadic per‑cell appearance.
- `--min-size <bp>` / `--max-size <bp>` (defaults: `50,000` / `5,000,000`): CNV sizes (log‑uniform).  
  **Effect:** Larger events experience larger absolute bp shifts under the same variance; mixed sizes add realism.

### Subclonal Rarity Model
- `--p-geom <float>` (default `0.45`): Geometric‑like decay for subclonal presence. A cell‑count `k` is drawn with P(k) ∝ p·(1−p)^(k−1); then `true_freq = k / n_cells`.  
  **Empirically:**  
  - Higher `p` (e.g., 0.6–0.8) → mostly single‑cell subclonals.  
  - Lower `p` (e.g., 0.2–0.3) → subclonals more often appear in multiple cells.
- `--max-subclonal-cells <int>` (default `10`): Cap on `k` for subclonal events.  
  **Effect:** Limits how “common” subclonals can become.

### Run‑Level Variance (Single Knob)
- `--variance <1–100>` (default `50`): Controls **positional shift** and **size distortion** for *all* observed events.  
  Internally maps linearly so that at `100` you see up to roughly:  
  - **Positional shift** ≈ 30% of the event length  
  - **Size change** (truncate/extend) ≈ 70% of the event length  
  **Interpretation:**  
  - `1` → **no variance** (observed = truth).  
  - `10–25` → small shifts (often sub‑kb for 50–200 kb events), small size changes.  
  - `40–60` → moderate shifts (hundreds bp–few kb), 10–20% size variation.  
  - `70–100` → aggressive; multi‑kb shifts and 30–70% size changes; frequent partial overlaps.

### High‑Variance (HV) Regions (Hotspots)
- `--hv-n-regions <int>` (default `2`): Number of HV regions to create.
- `--hv-radius <bp>` (default `75,000`): Each region spans `[center−radius, center+radius]`.
- `--hv-min-events <int>` (default `8`): **Guarantees** ≥ one HV region contains **at least this many** canonical events.  
  **Empirical:** Produces dense clusters of overlapping CNVs; great for testing callers in complex loci.  
  - See `high_variance_regions.bed` for HV spans and participating `event_ids`.  
  - See `canonical_events.bed` for each event’s `hv_regions` membership.

### Reproducibility & Output
- `--seed <int>` (default `1337`): Controls random draws (event positions, frequencies, HV placement, jitter).  
- `--outdir <dir>` (default `cnv_sim_run`): Output folder.

---

## How It Works (Overview)

1. **Canonical events** are created:  
   - Clonal events get `true_freq ~ U(0.80, 0.99)`  
   - Subclonal events draw a cell‑count `k` from geometric‑like decay → `true_freq = k / n_cells`  
   - Positions and sizes are random (log‑uniform), overlaps allowed

2. **HV regions** are placed:  
   - At least one HV region is guaranteed to contain ≥ `--hv-min-events` canonical events  
   - Other HV regions get a smaller number of overlapping events probabilistically

3. **Per‑cell observations**:  
   - For each canonical event, we choose `round(true_freq × n_cells)` cells (without replacement)
   - For each chosen cell, we **jitter** the truth interval by run‑level `--variance`:
     - **Positional shift** moves the whole interval
     - **Size change** truncates/extends around the center
   - Larger events shift more in absolute bp terms under the same variance

4. **Outputs** are written to `--outdir` in BED format.

---

## IGV Usage

- Load `canonical_events.bed` as the truth track.  
- Load one or more `per_cell/*.bed` to compare observed vs truth.  
- Add `high_variance_regions.bed` to highlight hotspots.  
- Tip: Color by `type`, enable labels, and sort by `start` for easier review.

---

## Reproducible Recipes

**Low‑noise benchmark**  
```bash
./run_cnv_sim.sh \
  --chrom-sizes chr1.sizes --chrom chr1 \
  --n-cells 200 --n-clonal 20 --n-subclonal 30 \
  --variance 15 --hv-n-regions 1 --hv-min-events 6 \
  --seed 7 --outdir cnv_low_noise
```

**Realistic single‑cell test**  
```bash
./run_cnv_sim.sh \
  --chrom-sizes chr1.sizes --chrom chr1 \
  --n-cells 200 --n-clonal 15 --n-subclonal 35 \
  --variance 50 --hv-n-regions 2 --hv-min-events 10 \
  --seed 42 --outdir cnv_medium_noise
```

**Stress test (high complexity)**  
```bash
./run_cnv_sim.sh \
  --chrom-sizes chr1.sizes --chrom chr1 \
  --n-cells 300 --n-clonal 20 --n-subclonal 60 \
  --variance 90 --hv-n-regions 3 --hv-min-events 15 --hv-radius 100000 \
  --seed 123 --outdir cnv_high_noise
```

---

## Troubleshooting

- **Chromosome not found**: Ensure `--chrom` exists in `--chrom-sizes`.  
- **Deps missing**: Re‑run the wrapper; it creates `.venv/` and installs `numpy`, `pandas`.  
- **Variance feels off**: Adjust `--variance` (1=no variance; 70–100=aggressive).  
- **Subclonals too common/rare**: Tweak `--p-geom` and `--max-subclonal-cells`, or adjust `--n-subclonal`.

---

## Citation / Acknowledgement

If you use this simulator in a paper or benchmark, please cite your repository release and include:  
“CNV Region Simulator (single‑cell CNV regions with run‑level variance and hotspot modeling).”
