# abcSplotch R-INLA container workflow

The Python preprocessing and downstream tools are installed in the conda environment.
The R-INLA model should be run from an Apptainer image built from `containers/Dockerfile`.

## Build the Docker image

```bash
docker build -t abcsplotch-rinla:latest -f containers/Dockerfile .
```

## Convert Docker image to Apptainer SIF

```bash
apptainer build containers/abcsplotch-rinla.sif docker-daemon://abcsplotch-rinla:latest
```

## Optional: create a writable sandbox for development

```bash
apptainer build --sandbox containers/abcsplotch-rinla.sandbox docker-daemon://abcsplotch-rinla:latest
apptainer shell --writable containers/abcsplotch-rinla.sandbox
```

For reproducible releases, prefer a read-only `.sif` built directly from the Docker image/tag used for that release.
