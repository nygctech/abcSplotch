FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Core system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gnupg \
    ca-certificates \
    wget \
    curl \
    build-essential \
    gfortran \
    cmake \
    pkg-config \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    locales \
    && rm -rf /var/lib/apt/lists/*

# Locale
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# ---- Add CRAN apt repo (THIS IS THE KEY FIX) ----
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | gpg --dearmor -o /usr/share/keyrings/cran.gpg

RUN echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/" \
    > /etc/apt/sources.list.d/cran.list

# Install R (will pull latest 4.5.x)
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# R library location
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN mkdir -p $R_LIBS_USER

# Install rstan
RUN R --vanilla -e " \
  install.packages(c('Rcpp', 'RcppEigen'), repos='https://cloud.r-project.org'); \
  install.packages('rstan', repos='https://cloud.r-project.org', type='source'); \
"
RUN R --vanilla -e "library(rstan); rstan::stan_version()"

# Install fmesher FIRST (INLA depends on it) -- must be installed from testing repo for R 4.5.1
RUN R -e "install.packages( \
  'fmesher', \
  repos = c( \
    INLA = 'https://inla.r-inla-download.org/R/testing', \
    CRAN = 'https://cloud.r-project.org' \
  ), \
  type = 'source' \
)"

# Install INLA
RUN R --vanilla <<'EOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(
  "INLA",
  repos = c(
    INLA = "https://inla.r-inla-download.org/R/stable",
    CRAN = "https://cloud.r-project.org"
  ),
  type = "source",
  INSTALL_opts = "--no-test-load"
)
EOF
RUN R -e "library(INLA); print(inla.version())"

# Sanity check
RUN R -e "library(fmesher); library(INLA); inla.version()"

WORKDIR /work
CMD ["R"]