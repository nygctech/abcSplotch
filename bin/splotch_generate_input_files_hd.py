from __future__ import absolute_import, division, print_function

import sys
import os
import pickle
import logging
import argparse
import pathlib

import numpy
import pandas as pd
import scanpy as sc

import anndata as ad
import numpy as np
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.sparse import block_diag, diags, csr_matrix, coo_matrix

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

import splotch

from splotch.utils import get_variable_mappings, watershed_tissue_sections, \
  print_summary, read_aar_matrix, read_array, read_array_metadata, \
  detect_tissue_sections, get_spot_adjacency_matrix, get_counts, \
  generate_W_sparse, generate_column_labels, generate_dictionary, \
  get_tissue_section_spots, filter_arrays, n_elements_per_level, to_rdump, \
  read_cellcomp_matrix

from splotch.utils_visium import unique_annots_loupe, read_annot_matrix_loupe, \
  detect_tissue_sections_hex, get_tissue_section_spots_hex, \
  get_spot_adjacency_matrix_hex, visium_find_position_file

from splotch.utils_sc import celltype_beta_priors

from splotch.utils_visiumHD import get_region_adjacency_matrix, get_cell_adjacency_matrix, generate_dirctionary

# Initialize containers for tissue section data and CAR structures.
aar_matrix_list = []
mouse_mapping_list = []
levels_list = []
tissue_mapping_list = []
counts_list = []
size_factors_list = []
coordinates_list = []
N_spots_list = []
W_list = []
W_n_list = []
region_list = []
files_list = []
cellcomp_matrix_list = []
car = True
region = False

metadata_file = '/gpfs/commons/groups/innovation/sarah/004/abcSplotch/splotch/metadata.tsv'
# Read sample metadata with explicit NA handling.
metadata = pd.read_csv(metadata_file,header=0,sep='\t',keep_default_na=False,na_values=pd_nan_values)

n_levels = 1
if n_levels == 1:
    levels = list(metadata['Level_1'].unique())

# Build variable mappings for the R model based on hierarchical sample levels.
level_mappings,last_level_identifiers,conditions_to_variables = get_variable_mappings(
      '',metadata,levels,n_levels)    

# Determine the number of elements at each hierarchy depth.
N_levels = [0]*n_levels
n_elements_per_level(levels,0,N_levels)

adata = ad.read('/gpfs/commons/groups/innovation/sarah/004/abcSplotch/splotch/CRC14NT_FINAL.unified.h5ad')

# Global gene and cell type annotations used for all tissue sections.
genes = adata.var_names.to_numpy(dtype=str)
celltype_names = np.loadtxt("/gpfs/commons/groups/innovation/sarah/004/abcSplotch/splotch/cell_type_categories.txt", dtype=str)

N_genes = len(genes)
N_covariates = len(celltype_names)

median_sequencing_depth = 452.0

count_files =  metadata['Anndata_File'].to_numpy().astype(str)

offset = 1
for file in count_files[:3]:
    # Load the current tissue section and extract coordinate and count data.
    adata = ad.read(file)
    array_genes = adata.var_names.to_numpy(dtype=str)
    array_coordinates_float = adata.obs[['array_row', 'array_col']].to_numpy()
    array_coordinates_str = np.array([
    f"{int(r)}_{int(c)}"
    for r, c in array_coordinates_float
])
    array_counts = adata.X.toarray()
    array_counts_per_spot = array_counts.sum(axis = 1)
    array_metadata = metadata[metadata['Anndata_File'] == file]
    array_levels = [array_metadata['Level_%d'%(idx+1)].values[0] for idx in range(0,n_levels)]

    tissue_mapping = last_level_identifiers[str(array_levels)]

    # Convert region-celltype annotations to a one-hot annotation matrix.
    array_aar_matrix = adata.obs['region_celltype'].to_numpy(dtype=str) 
    array_aar_matrix = (array_aar_matrix[:, None] == celltype_names[None, :]).astype(int)

    minimum_sequencing_depth = 100.0
    good_spots = array_counts_per_spot >= minimum_sequencing_depth  

    # Filter out low-depth spots from the current tissue section.
    array_coordinates_str,array_coordinates_float,array_counts,array_counts_per_spot, array_aar_matrix = \
        filter_arrays(good_spots,coordinates_str=array_coordinates_str,
                      coordinates_float=array_coordinates_float,
                      counts=array_counts,counts_per_spot=array_counts_per_spot, aar_matrix=array_aar_matrix.T)


    if car:
        if region: 
            if region_list:
                offset = np.unique(region_list[-1])[-1] + 1
            region_ids, tissue_section_W, connected_spots = get_region_adjacency_matrix(coords = array_coordinates_float, offset = offset)   

            array_coordinates_str,array_coordinates_float,array_counts,array_counts_per_spot, array_aar_matrix = \
            filter_arrays(connected_spots,coordinates_str=array_coordinates_str,
                        coordinates_float=array_coordinates_float,
                        counts=array_counts,counts_per_spot=array_counts_per_spot, aar_matrix=array_aar_matrix)
        else:
            tissue_section_W = get_cell_adjacency_matrix(coords = array_coordinates_float)
            connected_spots = np.array(tissue_section_W.sum(0) > 0).squeeze()
            array_coordinates_str,array_coordinates_float,array_counts,array_counts_per_spot, array_aar_matrix = \
            filter_arrays(connected_spots,coordinates_str=array_coordinates_str,
                        coordinates_float=array_coordinates_float,
                        counts=array_counts,counts_per_spot=array_counts_per_spot, aar_matrix=array_aar_matrix)
        

    array_size_factors = array_counts_per_spot / median_sequencing_depth

    counts_list.append(array_counts)
    # link the current tissue section with the filename
    files_list.append(file)
    # number of spots on the current tissue section
    N_spots_list.append(len(array_coordinates_str))
    # coordinates of the spots on the current tissue section
    coordinates_list.append(array_coordinates_str)
    # size factors of the spots on the current tissue section
    size_factors_list.append(array_size_factors)
    # annotations of the spots on the current tissue section
    aar_matrix_list.append(array_aar_matrix)
    # levels of the current tissue section
    levels_list.append(array_levels)
    tissue_mapping_list.append(tissue_mapping)

    if car:
        # adjacency matrix of the spots on the current tissue section
        if region: 
            region_list.append(region_ids)
        W_list.append(tissue_section_W)
          # number of adjacent spot pairs
        W_n_list.append(int(tissue_section_W.sum()/2))
    

N_tissues = len(counts_list)

# Build the full data dictionary prior to R dump generation.
data = generate_dictionary(N_spots_list,N_tissues,N_covariates,
                             N_levels,
                             coordinates_list,size_factors_list,
                             aar_matrix_list,level_mappings,
                             tissue_mapping_list,
                             W_list,W_n_list,car,True,False,
                             True,cellcomp_matrix_list, region_list)

output_directory = '/gpfs/commons/groups/innovation/sarah/004/crc_genes_car_region_small'

for gene_idx in range(N_genes):

    data['counts']  = get_counts(gene_idx,N_tissues,counts_list)

    if (gene_idx+1) % 1000 == 0:
      logging.info('%d/%d'%(gene_idx+1,len(genes)))

    # create the directory if it does not exist
    bucket = (gene_idx+1) // 100
    dirpath = os.path.normpath(f'{output_directory}/{bucket}/')
    if not os.path.exists(dirpath):
      os.makedirs(dirpath)

    # create the input data file containing data for the current gene
    # using the R dump format
    to_rdump(data, os.path.normpath(f'{output_directory}/{bucket}/data_{gene_idx+1}.R'))

filenames_coordinates = generate_column_labels(files_list,coordinates_list)

# Save metadata and mapping information for downstream analysis.
info_dict = {'genes':genes,'filenames_and_coordinates':filenames_coordinates,'annotation_mapping':celltype_names,'beta_mapping':conditions_to_variables,
                'n_levels':n_levels,'metadata':metadata,'scaling_factor':median_sequencing_depth,'car':True,'zi':True,'nb':False}

logging.info('Finished')
output_directory = '/gpfs/commons/groups/innovation/sarah/004/crc_genes_car_region_small'
pickle.dump(info_dict,open(os.path.normpath('/gpfs/commons/groups/innovation/sarah/004/crc_genes_car_region_small/information.p'),'wb'))