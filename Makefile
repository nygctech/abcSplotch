SHELL := /bin/bash
CONDA_ENV ?= abcsplotch
DOCKER_IMAGE ?= abcsplotch-rinla:latest
SIF ?= containers/abcsplotch-rinla.sif
SANDBOX ?= containers/abcsplotch-rinla.sandbox

.PHONY: env install dev docker sif sandbox clean-build

env:
	mamba env create -f environment.yml || conda env create -f environment.yml

install:
	pip install .

dev:
	pip install -e .

docker:
	docker build -t $(DOCKER_IMAGE) -f containers/Dockerfile .

sif: docker
	apptainer build $(SIF) docker-daemon://$(DOCKER_IMAGE)

sandbox: docker
	apptainer build --sandbox $(SANDBOX) docker-daemon://$(DOCKER_IMAGE)

clean-build:
	rm -rf build dist *.egg-info
