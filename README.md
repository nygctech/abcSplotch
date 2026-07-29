# abcSplotch
Approximate Bayesian inference of spatial cellular gene expression conditioned on sample covariates. 

## Features:
- Uses generalized linear model to denoise gene expression, disentangling effects of cell type, sample covariate, and local spatial autcorrelation
- Quantifies cohort-level differential expression across sample covariates and/or cell types 
- Efficient inference via integrated nested Laplace approximation (INLA) scales to tens of millions of spots/cells
- Compatible with Visium HD, Visium v1/v2, and original STv1
- Enables deconvolution of multi-cellular data (Visium v1/v2, STv1) via per-spot cellular composition estimates

![abcSplotch model](abcSplotch.png)

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/nygctech/abcSplotch.git
cd abcSplotch
```

---

### 2. Python Environment Installation

The preprocessing and downstream analysis tools are distributed as a Python package and are intended to run inside a Conda/Mamba environment.

We strongly recommend using `mamba` (or `micromamba`) with the `conda-forge` channel.

Create the environment:

```bash
mamba env create -f environment.yml
```

Activate the environment:

```bash
conda activate abcsplotch
```

---

### 3. Install the Python Package

Install the repository in editable mode:

```bash
pip install -e . --no-deps
```

This will:

- install the `abcSplotch` Python package
- create command-line executables from the package entry points
- keep the installation linked to the local repository for development

The `--no-deps` flag prevents `pip` from attempting to reinstall scientific dependencies already managed by Conda/Mamba.

### 4. Verify installation

```bash
splotch_prepare_count_files --help
splotch_generate_input_files --help
```

or

```bash
python -c "import splotch"
```

---

## R-INLA / Apptainer Backend

The Bayesian inference backend relies on a containerized R-INLA environment distributed via Docker/Apptainer.

We recommend using the provided Apptainer image directly on HPC systems.

Many HPC environments prohibit Docker execution directly on compute nodes. The recommended workflow is therefore:

1. Build Docker image locally or on a build node (must have x86 architecture)
2. Convert to Apptainer `.sif`
3. Transfer `.sif` to the cluster
4. Execute with `apptainer exec`

---

### 1. Build the Docker image

From the `containers/` directory:

```bash
cd containers

docker build -t abcsplotch-rinla:latest .
```

---

### 2. Build the Apptainer image

```bash
apptainer build \
    abcsplotch-rinla_latest.sif \
    docker-daemon://abcsplotch-rinla:latest
```

This produces a portable, immutable Apptainer image suitable for HPC execution. The `.sif` image is self-contained and reproducible across systems supporting Apptainer/Singularity.

---

### Optional: Create a Writable Sandbox

For debugging or interactive development:

```bash
apptainer build \
    --sandbox abcsplotch-rinla.sandbox \
    abcsplotch-rinla_latest.sif
```

Note that sandbox containers are mutable and therefore less reproducible than `.sif` images.

---

## Data preprocessing

## Model execution

abcSplotch will automatically select the model mode depending on the structure of the input data file:
- Compositional mode: if `N_celltypes` and per-spot compositional data (`E`) are provided
- Spatial mode: if `car` flag is set and adjacency information (`W_sparse`)
- Regional mode: if `region_list` (mapping cells to spatial regions) is provided

Additionally, the number of hierarchical levels of the model will be read from `N_levels`. 

Informed priors for top-level beta terms can be set for each cell type using `beta_prior_mean` and `beta_prior_std` (default: Normal(0,2)).

#### Usage:

```bash
apptainer exec \
    <apptainer_image> \              
    Rscript abcSplotch/inla/abcsplotch.R \
    <input_path> \
    <output_dir> \
    [--draw-samples] \ 
    [--regional-precision prec]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `apptainer_image` | Path to the .sif or .sandbox Apptainer image |
| `input_path` | Path to the Rdump-formatted input for a gene of interest, e.g., `data_[GID].R`. |
| `output_dir` | Path to top-level directory in which to store output. To keep individual directories from overcrowding, `results_[GID].R` will be stored in `output_dir/[GID//100]` subdirectory |
| `--draw-samples` | Specify that the model should draw and store samples from the joint posterior, in addition to default mode of only computing and storing marginals for each parameter. This is expensive, and should only be done when analysis requires covariance between beta terms |
| `--regional-precision=prec` | Precision of the BYM2 component used in regional model. Should be between log(100) and log(400), with lower values yielding higher impact (default: log(200)) |

#### Examples

Run storing marginals only:

```bash
apptainer exec \
    abcsplotch-rinla.sif \              # or abcsplotch-rinla.sandbox
    Rscript abcSplotch/inla/abcsplotch.R \
    data_1.R \
    output_dir \
```

Run drawing and storing joint samples from posterior:

```bash
apptainer exec \
    abcsplotch-rinla.sif \              # or abcsplotch-rinla.sandbox
    Rscript abcSplotch/inla/abcsplotch.R \
    data_1.R \
    output_dir \
    --draw-samples
```

Run with reduced importance of regional spatial component:

```bash
apptainer exec \
    abcsplotch-rinla.sif \              # or abcsplotch-rinla.sandbox
    Rscript abcSplotch/inla/abcsplotch.R \
    data_1.R \
    output_dir \
    --regional-precision=6
```

## Downstream analysis & visualization

### Differential expression testing: `compute_contrasts.r`

Computes posterior contrasts between pairs of regression coefficients using either:

- the Gaussian approximation derived from posterior marginal summaries, or
- stored joint posterior samples.

The script supports both the effective hierarchical coefficients (`beta_l1`, `beta_l2`, `beta_l3`) and the raw hierarchy components (`beta_l1`, `delta_l2`, `delta_l3`).

#### Usage

```bash
Rscript compute_contrasts.r \
    <results_rds> \
    <output_prefix> \
    <delta_or_NA> \
    <use_samples_TRUE_FALSE> \
    <beta1:beta2> [beta3:beta4 ...] \
    [--beta-set=effective|components]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `results_rds` | Path to the saved abcSplotch results object (`.rds`). |
| `output_prefix` | Prefix for the output `.csv` and `.rds` files. |
| `delta_or_NA` | Minimum effect size threshold. Use `NA` to disable effect-size thresholding. |
| `use_samples_TRUE_FALSE` | `TRUE` to compute contrasts from stored joint posterior samples; `FALSE` to use the Gaussian approximation from marginal summaries. |
| `beta1:beta2` | Contrast specification (`beta1 - beta2`). Multiple contrasts may be supplied. |
| `--beta-set` | *(Optional)* Which coefficients to use. Defaults to `effective`. |

#### Beta sets

| Option | Uses |
|--------|------|
| `effective` *(default)* | `res$beta_effective` (`beta_l1`, `beta_l2`, `beta_l3`) |
| `components` | `res$beta_components` (`beta_l1`, `delta_l2`, `delta_l3`) |

#### Examples

Level-2 effective contrast from Gaussian marginals:

```bash
Rscript compute_contrasts.r \
    results.rds \
    crc_contrasts \
    2 \
    FALSE \
    beta_l2_c1_m6:beta_l2_c1_m2
```

Level-3 effective contrast using joint posterior samples:

```bash
Rscript compute_contrasts.r \
    results.rds \
    crc_contrasts \
    2 \
    TRUE \
    beta_l3_c8_m6:beta_l3_c8_m2
```

Contrast between raw level-2 hierarchy components:

```bash
Rscript compute_contrasts.r \
    results.rds \
    delta_contrasts \
    NA \
    TRUE \
    delta_l2_c8_m6:delta_l2_c8_m2 \
    --beta-set=components
```

#### Outputs

The script writes

- `<output_prefix>.csv`
- `<output_prefix>.rds`

containing posterior means, standard deviations, credible intervals, posterior probabilities, and (optionally) effect-size threshold statistics for each requested contrast.

---

### Posterior plotting: `plot_beta_kde.R`

Plots posterior distributions for one or more regression coefficients.

Two plotting modes are available:

- **marginal** — plots the INLA posterior marginal densities.
- **joint** — plots kernel density estimates (KDEs) of stored joint posterior samples.

The script can visualize either the effective hierarchical coefficients or the raw hierarchy components.

#### Usage

```bash
Rscript plot_beta_kde.R \
    <results_rds> \
    <mode> \
    <beta_names_comma_sep> \
    [labels_comma_sep] \
    [output_png] \
    [--beta-set=effective|components]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `results_rds` | Path to the saved abcSplotch results object. |
| `mode` | Either `marginal` or `joint`. |
| `beta_names_comma_sep` | Comma-separated coefficient names. |
| `labels_comma_sep` | *(Optional)* Legend labels. Defaults to coefficient names. |
| `output_png` | *(Optional)* Output PNG filename. If omitted, the plot is shown interactively. |
| `--beta-set` | *(Optional)* Which coefficient collection to use. Defaults to `effective`. |

#### Beta sets

| Option | Uses |
|--------|------|
| `effective` *(default)* | `res$beta_effective` (`beta_l1`, `beta_l2`, `beta_l3`) |
| `components` | `res$beta_components` (`beta_l1`, `delta_l2`, `delta_l3`) |

#### Examples

Plot posterior marginals for two effective level-2 coefficients:

```bash
Rscript plot_beta_kde.R \
    results.rds \
    marginal \
    beta_l2_c1_m1,beta_l2_c2_m1 \
    Condition1,Condition2
```

Plot KDEs from joint posterior samples:

```bash
Rscript plot_beta_kde.R \
    results.rds \
    joint \
    beta_l3_c5_m1,beta_l3_c8_m1 \
    Young,Old \
    beta_l3.png
```

Plot raw hierarchy components:

```bash
Rscript plot_beta_kde.R \
    results.rds \
    joint \
    delta_l2_c1_m1,delta_l2_c2_m1 \
    Delta1,Delta2 \
    delta_l2.png \
    --beta-set=components
```

#### Notes

- `marginal` mode visualizes the one-dimensional INLA posterior marginals.
- `joint` mode visualizes kernel density estimates computed from the stored joint posterior samples.
- Multiple coefficients are overlaid on a common set of axes for comparison.

---
