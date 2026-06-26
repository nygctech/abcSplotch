import anndata as ad
import numpy as np
import os
import pandas as pd
import logging
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.sparse import block_diag, diags, csr_matrix, coo_matrix
import pickle

# Common pandas values interpreted as missing / NA in metadata parsing.
pd_nan_values = ['','#N/A','#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NA', 'NULL', 'NaN', 'n/a', 'nan', 'null']

def get_region_adjacency_matrix(coords=None, offset = 0):
    bin_size = 50

    grid_x = np.floor(coords[:,0] / bin_size).astype(int)
    grid_y = np.floor(coords[:,1] / bin_size).astype(int)

    regions = np.vstack([grid_x, grid_y]).T

    unique_regions, region_idx = np.unique(
        regions, axis=0, return_inverse=True
    )

    n_regions = len(unique_regions)

    region_lookup = {
        tuple(coord): rid
        for rid, coord in enumerate(unique_regions)
    }

    neighbor_offsets = [
        (-1,-1), (-1,0), (-1,1),
        (0,-1),         (0,1),
        (1,-1), (1,0), (1,1)
    ]

    region_neighbors = {}

    for rid, (gx, gy) in enumerate(unique_regions):
        neighbors = []
        for dx, dy in neighbor_offsets:
            nbr_coord = (gx+dx, gy+dy)
            if nbr_coord in region_lookup:
                neighbors.append(region_lookup[nbr_coord])
        region_neighbors[rid] = neighbors

    A = np.zeros((n_regions, n_regions), dtype=int)

    for r, nbrs in region_neighbors.items():
        A[r, nbrs] = 1

    A = np.maximum(A, A.T)  # important for symmetry

    degree = A.sum(axis=1)
    keep = degree > 0


    keep_cells = keep[region_idx]
    kept_regions = np.where(keep)[0]
    A = A[np.ix_(keep, keep)]

    old_to_new = -np.ones(n_regions, dtype=int)
    old_to_new[kept_regions] = np.arange(len(kept_regions))

    region_idx_kept = old_to_new[region_idx[keep_cells]]
    region_idx_kept = region_idx_kept + offset

    return region_idx_kept, A, keep_cells


def get_cell_adjacency_matrix(adata=None, coords=None, k=4):
    if coords is None or coords.size == 0:
        coords = adata.obs[['array_row', 'array_col']].to_numpy()

    n = coords.shape[0]

    if n < 2:
        return sparse.csr_matrix((n, n))

    k_eff = min(k, n - 1)

    tree = cKDTree(coords)

    distances, indices = tree.query(coords, k=k_eff + 1)
    indices = indices[:, 1:]

    rows_list = np.repeat(np.arange(n), k_eff)
    cols_list = indices.ravel()

    data = np.ones(len(cols_list), dtype=np.int8)

    adj = sparse.csr_matrix(
        (data, (rows_list, cols_list)),
        shape=(n, n)
    )

    # symmetrize
    adj = adj.maximum(adj.T)

    return adj




def eigs_square_chunking(coords_str, chunksize=10000, k=4, overlap=0):
    """
    Compute eigenvalues of normalized adjacency in spatial chunks (robust).
    """

    # --- parse coordinates ---
    coords_int = np.array([tuple(map(int, c.split('_'))) for c in coords_str])
    xdim = coords_int[:, 0].max()
    ydim = coords_int[:, 1].max()

    eigs_all = []

    wc = int(np.rint(chunksize**0.5))

    for i in range(0, xdim + 1, wc):
        for j in range(0, ydim + 1, wc):

            # --- include overlap to fix boundary issues ---
            xinds = (coords_int[:, 0] >= i - overlap) & (coords_int[:, 0] < i + wc + overlap)
            yinds = (coords_int[:, 1] >= j - overlap) & (coords_int[:, 1] < j + wc + overlap)
            mask = xinds & yinds
            coords_sub = coords_int[mask]

            n_sub = len(coords_sub)
            if n_sub == 1:
                eigs_all.append(np.array([0.0]))
                continue
            if n_sub > chunksize:
                logging.warning(
                    'Improper use of chunking eigenvalue solution -- neighboring spots should be unit-distance apart'
                )
            k_eff = min(k, n_sub - 1)
            if k_eff <= 0:
                continue

            W_sub = get_cell_adjacency_matrix(coords=coords_sub, k=k_eff)
            W_sub = W_sub.maximum(W_sub.T)
            D_sub = W_sub.sum(1).A1

            # --- remove zero-degree nodes (prevents inf) ---
            valid = D_sub > 0
            if not np.any(valid):
                continue

            W_sub = W_sub[valid][:, valid]
            D_sub = D_sub[valid]
            D_v = diags(1.0 / np.sqrt(D_sub), 0, format='csr')
            R = D_v @ W_sub @ D_v

            eigs = np.linalg.eigvalsh(R.toarray())

            eigs_all.append(eigs)

    if len(eigs_all) == 0:
        raise ValueError("No eigenvalues computed — check inputs")

    eigs_all = np.concatenate(eigs_all)
    eigs_all.sort()

    if len(eigs_all) != len(coords_str):
        logging.warning(
            f"Eigenvalue count ({len(eigs_all)}) != number of coords ({len(coords_str)}) "
            f"(likely due to isolated nodes)"
        )

    return eigs_all

def generate_dictionary(N_spots_list,N_tissues,N_covariates,
                        N_levels,coordinates_list,
                        size_factors_list,aar_matrix_list, 
                        level_mappings,tissue_mapping_list,
                        W_list,W_n_list,car,region,zi,nb,
                        compositional,cellcomp_matrix_list, region_list):

    data = {'N_spots': N_spots_list,
            'N_tissues': N_tissues,
            'N_covariates': N_covariates,
            'tissue_mapping': tissue_mapping_list,
            'N_levels': len(N_levels),
            'zi': 1*zi,
            'nb': 1*nb,
            'car': 1*car}

    for idx in range(0,len(N_levels)):
        data['N_level_%d'%(idx+1)]  = N_levels[idx]
    for idx in range(len(N_levels),3):
        data['N_level_%d'%(idx+1)]  = 0

    for idx in range(0,len(N_levels)-1):
        data['level_%d_mapping'%(idx+2)] = level_mappings[idx]

    for idx in range(len(N_levels)-1+2,3+1):
        data['level_%d_mapping'%(idx)] = []

    concatenated_size_factors = []
    for tissue_idx in range(0,N_tissues):
        concatenated_size_factors = concatenated_size_factors + \
        list(size_factors_list[tissue_idx])
    data['size_factors'] = concatenated_size_factors

    concatenated_D = []
    for tissue_idx in range(0,N_tissues):
        concatenated_D = concatenated_D + [numpy.where(tissue_section_aar_matrix)[0][0]+1 \
        for tissue_section_aar_matrix in aar_matrix_list[tissue_idx].T]
    data['D'] = concatenated_D

    if car:
        if region:
          concatenated_region_ids = []
          for tissue_idx in range(0,N_tissues):
            concatenated_region_ids = concatenated_region_ids + \
            list(region_list[tissue_idx])
          data['region_list'] = concatenated_region_ids
        data['W_n']  = [sum(W_n_list)]
        W = block_diag(W_list,format='csr')
        data['W_sparse'] = generate_W_sparse(sum(N_spots_list),data['W_n'][0],W)
        data['D_sparse'] = W.sum(1).A1.astype(int)

        eigval_list = []
        for t, (W_tissue, coords_tissue) in enumerate(zip(W_list, coordinates_list)):
            if W_tissue.shape[0] < 10000:
              W = block_diag([W_tissue],format='csr')
              D_sparse_tissue = W.sum(1).A1.astype(int)
              D_v = diags(1.0/numpy.sqrt(D_sparse_tissue),0,format='csr')
              evs = numpy.linalg.eigvalsh(D_v.dot(W).dot(D_v).toarray())
            else:
              logging.info('Tissue %d exceeds 10k connected points -- CAR prior eigenvalues calculated by chunking' % t)
              evs = eigs_square_chunking(coords_tissue, 10000)
            eigval_list.append(evs)
        data['eig_values'] = numpy.sort(numpy.concatenate(eigval_list))
        
        V, U = sparse_to_csr_indptr(data['W_sparse'], numpy.sum(data['N_spots']))
        data['U'] = U 
        data['V'] = V

    else:
        data['W_n'] = []
        data['W_sparse'] = numpy.zeros((0,0))
        data['D_sparse'] = []
        data['eig_values'] = []
        data['U'] = []
        data['V'] = []

    return data



def filter_arrays(indices,coordinates_str=None,coordinates_float=None,
                  counts=None,counts_per_spot=None,size_factors=None,
                  aar_matrix=None,W=None,
                  cellcomp_matrix=None):
  output = []

  if coordinates_str is not None:
    output.append(coordinates_str[indices])
  if coordinates_float is not None:
    output.append(coordinates_float[indices,:])
  if counts is not None:
    output.append(counts[indices,:])
  if counts_per_spot is not None:
     output.append(counts_per_spot[indices])
  if size_factors is not None:
     output.append(size_factors[indices])
  if aar_matrix is not None:
    output.append(aar_matrix[:,indices])
  if W is not None:
    output.append(W[indices,:][:,indices])
  if cellcomp_matrix is not None:
    output.append(cellcomp_matrix[:,indices])

  return output


def generate_column_labels(files_list,coordinates_list):
  filenames = [[files_list[r]]*len(coordinates_list[r]) \
    for r in range(0,len(coordinates_list))]
  filenames = [foo for bar in filenames for foo in bar]
  coordinates = numpy.hstack((coordinates_list[:]))
  filenames_coordinates =  list(zip(*[filenames,coordinates]))

  return filenames_coordinates


def get_variable_mappings(count_files,metadata,
                             levels,n_levels):

  conditions_to_variables = {'beta_level_%d'%(idx+1):[] for idx in range(0,n_levels)}

  level_mappings = [[] for _ in range(0,n_levels-1)]
  last_level_identifiers = {}

  if n_levels == 3:
    levels_1 = levels.keys()
  
    level_2_idx = 1
    level_3_idx = 1
    for level_1_idx,level_1 in enumerate(levels_1,start=1):
      logging.info('beta_level_1[%d] := %s'%(level_1_idx,level_1))
      conditions_to_variables['beta_level_1'].append('%s'%(level_1))
  
      for level_2 in levels[level_1]:
        logging.info('beta_level_2[%d] := %s %s'%(level_2_idx,level_1,level_2))
        conditions_to_variables['beta_level_2'].append('%s %s'%(level_1,level_2))
  
        level_mappings[0].append(level_1_idx)
  
        for level_3 in levels[level_1][level_2]:
          logging.info('beta_level_3[%d] := %s %s %s'%(level_3_idx,level_1,level_2,level_3))
          conditions_to_variables['beta_level_3'].append('%s %s %s'%(level_1,level_2,level_3))
          last_level_identifiers[str([level_1,level_2,level_3])] = level_3_idx
          level_mappings[1].append(level_2_idx)
  
          level_3_idx += 1
  
        level_2_idx += 1

  elif n_levels == 2:
    levels_1 = levels.keys()
  
    level_2_idx = 1
    for level_1_idx,level_1 in enumerate(levels_1,start=1):
      logging.info('beta_level_1[%d] := %s'%(level_1_idx,level_1))
      conditions_to_variables['beta_level_1'].append('%s'%(level_1))
  
      for level_2 in levels[level_1]:
        logging.info('beta_level_2[%d] := %s %s'%(level_2_idx,level_1,level_2))
        conditions_to_variables['beta_level_2'].append('%s %s'%(level_1,level_2))
  
        level_mappings[0].append(level_1_idx)
  
        last_level_identifiers[str([level_1,level_2])] = level_2_idx
  
        level_2_idx += 1

  elif n_levels == 1:
    levels_1 = levels
  
    for level_1_idx,level_1 in enumerate(levels_1,start=1):
      logging.info('beta_level_1[%d] := %s'%(level_1_idx,level_1))
      conditions_to_variables['beta_level_1'].append('%s'%(level_1))
  
      last_level_identifiers[str([level_1])] = level_1_idx
  
  return level_mappings,last_level_identifiers,conditions_to_variables


def n_elements_per_level(node,n,tmp=None):
  if tmp is None:
    tmp = [0]*10
  tmp[n] += len(node)
  if not isinstance(node,dict):
      return tmp
  for key, item in node.items():
    if isinstance(item,dict):
      n_elements_per_level(item,n+1,tmp)
    else:
      tmp[n+1] += len(item)
  return tmp




def sparse_to_csr_indptr(W_sparse, N_spots):
  """Convert an edge list into CSR indices and indptr arrays."""
  r_inds = numpy.concatenate((W_sparse[:,0]-1, W_sparse[:,1]-1))  
  c_inds = numpy.concatenate((W_sparse[:,1]-1, W_sparse[:,0]-1))

  data = numpy.ones_like(r_inds, dtype=np.int64)
  cmat = csr_matrix((data, (r_inds, c_inds)), shape=(N_spots, N_spots))
  return cmat.indices+1, cmat.indptr+1


def to_rdump(data,filename):
    with open(filename,'w') as f:
        for key in data:
            tmp = numpy.asarray(data[key])
            if len(tmp.shape) == 0:
                f.write('%s <- %s\n'%(key,str(tmp)))
            elif len(tmp.shape) == 1:
                f.write('%s <-\n c(%s)\n'%(key,
                    ','.join(map(str,tmp))))
            else:
                f.write('%s <-\n structure(c(%s), .Dim=c(%s))\n'%(key,
                    ','.join(map(str,tmp.flatten('F'))),
                    ','.join(map(str,tmp.shape))))


def get_counts(gene_idx,N_tissues,counts_list):
  concatenated_counts = []
  for tissue_idx in range(0,N_tissues):
    concatenated_counts = concatenated_counts + list(counts_list[tissue_idx][:,gene_idx])

  return concatenated_counts


def generate_W_sparse(N,W_n,W):
  """Generate an undirected edge list from a sparse adjacency matrix."""

  spr, spc = csr_matrix(W).nonzero()
  W_sparse = numpy.vstack((spr+1, spc+1)).T
  W_sparse = numpy.array([x for x in W_sparse if x[0] < x[1]])

  if not W_sparse.shape[0] == W_n:
    logging.critical('Adjacency matrix W does not contain expected number of adjacent pairs W_n!')
    import sys
    sys.exit(1)

  return W_sparse

def get_region_adjacency_matrix(coords=None, offset = 0):
    """Bin spatial coordinates into regions and build region adjacency.

    This function groups spots into spatial bins, identifies neighboring regional
    bins, removes isolated regions, and remaps region indices for downstream use.
    """
    bin_size = 50

    # Assign each spot to a grid cell based on its x/y coordinates.
    grid_x = np.floor(coords[:,0] / bin_size).astype(int)
    grid_y = np.floor(coords[:,1] / bin_size).astype(int)

    regions = np.vstack([grid_x, grid_y]).T

    unique_regions, region_idx = np.unique(
        regions, axis=0, return_inverse=True
    )

    n_regions = len(unique_regions)

    region_lookup = {
        tuple(coord): rid
        for rid, coord in enumerate(unique_regions)
    }

    # Define the 8-connected neighborhood offsets for adjacent bins.
    neighbor_offsets = [
        (-1,-1), (-1,0), (-1,1),
        (0,-1),         (0,1),
        (1,-1), (1,0), (1,1)
    ]

    region_neighbors = {}

    for rid, (gx, gy) in enumerate(unique_regions):
        neighbors = []
        for dx, dy in neighbor_offsets:
            nbr_coord = (gx+dx, gy+dy)
            if nbr_coord in region_lookup:
                neighbors.append(region_lookup[nbr_coord])
        region_neighbors[rid] = neighbors

    # Create the region adjancency matrix from neighbor relationships.
    A = np.zeros((n_regions, n_regions), dtype=int)

    for r, nbrs in region_neighbors.items():
        A[r, nbrs] = 1

    A = np.maximum(A, A.T)  # important for symmetry

    degree = A.sum(axis=1)
    keep = degree > 0

    # Keep only spots that belong to non-isolated regions.
    keep_cells = keep[region_idx]

    # Filter adjacency to retained regions only.
    kept_regions = np.where(keep)[0]
    A = A[np.ix_(keep, keep)]

    # Remap remaining regions to a contiguous index set.
    old_to_new = -np.ones(n_regions, dtype=int)
    old_to_new[kept_regions] = np.arange(len(kept_regions))

    region_idx_kept = old_to_new[region_idx[keep_cells]]

    # OPTIONAL: apply offset ONLY at the very end.
    region_idx_kept = region_idx_kept + offset

    return region_idx_kept, A, keep_cells

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