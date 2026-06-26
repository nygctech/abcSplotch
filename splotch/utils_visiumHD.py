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