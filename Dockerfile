# Use a Rocker base image with R 4.3.2
FROM rocker/r-ver:4.3.2

# Set environment for non-interactive installs
ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=/work/Rlibs

# Create the personal library inside the container
RUN mkdir -p /work/Rlibs

# Install system libraries required by fmesher and INLA
RUN apt-get update && apt-get install -y \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install R packages in the personal library
RUN R -e "install.packages('Rcpp', repos='https://cloud.r-project.org', lib=Sys.getenv('R_LIBS_USER'))"
RUN R -e "install.packages('fmesher', repos='https://inla.r-inla-download.org/R/stable', type='source', lib=Sys.getenv('R_LIBS_USER'))"
RUN R -e "install.packages('INLA', repos='https://inla.r-inla-download.org/R/stable', type='source', lib=Sys.getenv('R_LIBS_USER'))"

# Verify that packages load
RUN R -e "library(fmesher, lib.loc=Sys.getenv('R_LIBS_USER')); library(INLA, lib.loc=Sys.getenv('R_LIBS_USER')); print(inla.version())"

# Set working directory
WORKDIR /work

# Default command
CMD ["R"]

# Set PATH for R packages
ENV PATH=$R_LIBS_USER:$PATH