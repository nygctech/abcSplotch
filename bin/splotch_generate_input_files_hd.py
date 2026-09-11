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

from splotch.utils_visiumHD import get_region_adjacency_matrix, get_cell_adjacency_matrix, generate_dictionary

pd_nan_values = ['','#N/A','#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NA', 'NULL', 'NaN', 'n/a', 'nan', 'null']


def generate_splotch_inputs(count_files,metadata_file,
                            median_sequencing_depth,
                            n_levels,minimum_sequencing_depth,
                            output_directory,car,region=None,
                            annot_col='annotation',
                            x_col='x', y_col='y'
                           ):
    aar_matrix_list = []
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


    # Read sample metadata with explicit NA handling.
    metadata = pd.read_csv(metadata_file,header=0,sep='\t',keep_default_na=False,na_values=pd_nan_values)

    # three levels
    if n_levels == 3:
        levels = {key:{} for key in list(metadata['Level 1'].unique())}
        for level_1 in levels:
            levels[level_1] = {key: {} for key in list(metadata[metadata['Level 1'] == level_1]['Level 2'].unique())}
            for level_2 in levels[level_1]:
                levels[level_1][level_2] = list(metadata[(metadata['Level 1'] == level_1) & (metadata['Level 2'] == level_2)]['Level 3'].unique())

    # two levels
    elif n_levels == 2:
        levels = {key:{} for key in list(metadata['Level 1'].unique())}
        for level_1 in levels:
            levels[level_1] = list(metadata[metadata['Level 1'] == level_1]['Level 2'].unique())

    # one level
    elif n_levels == 1:
        levels = list(metadata['Level 1'].unique())

    # Build variable mappings for the R model based on hierarchical sample levels.
    level_mappings,last_level_identifiers,conditions_to_variables = get_variable_mappings(
        '',metadata,levels,n_levels)    

    # Determine the number of elements at each hierarchy depth.
    N_levels = [0]*n_levels
    n_elements_per_level(levels,0,N_levels)

    # Just pick the first
    adata = ad.read(count_files[0])

    # Find unique annotations and shared genes across all arrays
    celltype_names = np.unique(adata.obs[annot_col].values)
    shared_genes = adata.var_names
    if len(count_files) > 1:
        for cfile in count_files[1:]:
            adata = sc.read_h5ad(cfile)
            for col in [annot_col, x_col, y_col]:
                if col not in adata.obs.columns:
                    logging.critical("%s missing from counts file %s" % (col, cfile))
                    sys.exit(1)
            celltype_names = np.union1d(celltype_names, np.unique(adata.obs[annot_col].values))
            shared_genes = np.intersect1d(shared_genes, adata.var_names)

    N_genes = len(shared_genes)
    N_covariates = len(celltype_names)
    logging.info("Detected %d unique annotations and %d shared genes across all tissue sections" % (N_covariates, N_genes))        

    offset = 1
    for file in count_files:
        # Check if file has a metadata entry
        if file not in metadata['Anndata_File'].values:
            logging.warning("Count file %s will not be considered as it has no metadata entry" % cfile)
            continue
        if np.sum(metadata['Anndata_File']==cfile) > 1:
            logging.warning("Count file %s will not be considered as it maps to MULTIPLE metadata entries!" % cfile)
            continue
        
        # Load the current tissue section and extract coordinate and count data.
        adata = ad.read(file)
        adata = adata[:, shared_genes]
        
        array_coordinates_float = adata.obs[[y_col, x_col]].to_numpy()
        array_coordinates_str = np.array([
            f"{int(r)}_{int(c)}" for r, c in array_coordinates_float
        ])
        array_counts = adata.X.toarray()
        array_counts_per_spot = array_counts.sum(axis = 1)
        array_metadata = metadata[metadata['Anndata_File'] == file]
        array_levels = [array_metadata['Level_%d'%(idx+1)].values[0] for idx in range(0,n_levels)]

        tissue_mapping = last_level_identifiers[str(array_levels)]

        # Convert region-celltype annotations to a one-hot annotation matrix.
        array_aar_matrix = adata.obs[annot_col].to_numpy(dtype=str) 
        array_aar_matrix = (array_aar_matrix[:, None] == celltype_names[None, :]).astype(int)

        good_spots = array_counts_per_spot >= minimum_sequencing_depth  

        # Filter out low-depth spots from the current tissue section.
        array_coordinates_str,array_coordinates_float,array_counts,array_counts_per_spot, array_aar_matrix = \
            filter_arrays(good_spots,coordinates_str=array_coordinates_str,
                        coordinates_float=array_coordinates_float,
                        counts=array_counts,counts_per_spot=array_counts_per_spot, aar_matrix=array_aar_matrix.T)


        if car:
            if region is not None: 
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
            if region is not None: 
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
    pickle.dump(info_dict,open(os.path.normpath('%s/information.p'%(output_directory)),'wb'))

if __name__ == '__main__':

  def levels_type(x):
    x = int(x)
    if x not in [1,2,3]:
      raise argparse.ArgumentTypeError("Number of levels should be either 1, 2 or 3")
    return x

  parser = argparse.ArgumentParser(
    description='A script for generating input data for Splotch')
  parser.add_argument('-c','--count_files',action='store',
                      dest='count_files',type=str,nargs='+',
                      required=True,help='list of count files in h5ad-serialized AnnData format')
  parser.add_argument('-m','--metadata_file',action='store',
                      dest='metadata_file',type=str,required=True,
                      help='metadata filename')
  parser.add_argument('-s','--scaling_factor',action='store',
                      dest='scaling_factor',type=float,required=True,
                      help='scaling_factor (e.g. median sequencing depth over spots)')
  parser.add_argument('-l','--n_levels',action='store',
                      dest='n_levels',type=levels_type,required=False,default=3,
                      help='number of levels in the linear model (default is 3)')
  parser.add_argument('-d','--min_depth',action='store',
                      dest='minimum_sequencing_depth',type=float,
                      required=False,default=100.0,
                      help='minimum number of UMIs per spot (default is 100)')
  parser.add_argument('-n','--no-car',action='store_false',dest='car',required=False,
                      help='disable the conditional autoregressive prior')
  parser.add_argument('-o','--output_directory',action='store',
                      dest='output_directory',type=str,
                      required=False,default='abcsplotch_inputs',
                      help='output_directory (default is "abcsplotch_inputs")')
  parser.add_argument('-r', '--region', type=float, default=None, dest='region',required=False,
                      help='Model spatial effects across regions of specified width')
  parser.add_argument('-x', '--x_col', type=str, default='x', dest=x_col, required=False,
                      help='Column in adata.obs in which x-coordinate of spatial measurements is stored')
  parser.add_argument('-y', '--y_col', type=str, default='y', dest=y_col, required=False,
                      help='Column in adata.obs in which y-coordinate of spatial measurements is stored')
  parser.add_argument('-a', '--annot_col', type=str, default='annotation', dest=annot_col, required=False,
                      help='Column in adata.obs in which annotation of spatial measurements is stored')

  parser.set_defaults(car=True)
  options = parser.parse_args()

  # check that the supplied read count files exist
  for filename in options.count_files:
    if not os.path.isfile(filename):
      logging.critical('Count file %s does not exist!'%(filename))
      sys.exit(1)
  # check that the supplied metadata file exists
  if not os.path.isfile(options.metadata_file):
    logging.critical('Metadata file %s does not exist!'%(options.metadata_file))
    sys.exit(1)


  # generate the data dictionaries 
  generate_splotch_inputs(options.count_files,options.metadata_file,
                          options.scaling_factor,options.n_levels,
                          options.minimum_sequencing_depth,options.output_directory,
                          options.car, region=options.region, annot_col=options.annot_col,
                          x_col=options.x_col, y_col=options.y_col)
                         
  sys.exit(0)