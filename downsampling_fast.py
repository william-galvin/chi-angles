from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer
import h5py
from tqdm import tqdm


so_file = "/gscratch/scrubbed/wgalvin/python/utils.so"
c_functions = CDLL(so_file)

filename = "/gscratch/scrubbed/wgalvin/clean_clone/protein_holography-pytorch/protein_holography_pytorch/neighborhoods_tiny.hdf5"

max_atoms = 1000
dt = np.dtype([
    ('res_id','S6', (6)), # S5, 5 (old) ; S6, 6 (new with 2ndary structure)
    ('atom_names', 'S4', (max_atoms)),
    ('elements', 'S1', (max_atoms)),
    ('res_ids', 'S6', (max_atoms, 6)), # S5, 5 (old) ; S6, 6 (new with 2ndary structure)
    ('coords', 'f8', (max_atoms, 3)),
    ('SASAs', 'f8', (max_atoms)),
    ('charges', 'f8', (max_atoms)),
])

BACKBONE_ATOMS = np.array([b' N  ', b' CA ', b' C  ', b' O  '])
bb_len = len(BACKBONE_ATOMS)
BACKBONE_ATOMS = BACKBONE_ATOMS.tobytes()
MAX_SITES = 10_000
MAX_CHAINS = 10

get_mask = c_functions.get_mask
get_mask.restype = None
get_mask.argtypes = [
    c_int, 
    c_int, 
    c_char_p,
    ndpointer(c_long, flags="C_CONTIGUOUS"),
    c_void_p, 
    ndpointer(c_bool, flags="C_CONTIGUOUS"),
    ndpointer(c_bool, flags="C_CONTIGUOUS"), 
    ndpointer(c_double, flags="C_CONTIGUOUS"),
    ndpointer(c_bool, flags="C_CONTIGUOUS"), 
    c_float, 
    c_void_p, 
    c_int,
    c_char,
    c_long
]


def get_mask_c(n, p):
    # unpadding only necessary because converting sites -> int
    unpad_mask = n['atom_names'] != b''
    res_ids = n['res_ids'][unpad_mask]
    chains = res_ids[:, 2].tobytes()
    sites = np.array(res_ids[:, 3], dtype=int)
    SIZE = len(res_ids)

    names = n['atom_names'].tobytes()

    rand = np.random.rand((max_atoms))

    seen = np.zeros(MAX_CHAINS * MAX_SITES, dtype=bool)
    removed = np.zeros(MAX_CHAINS * MAX_SITES, dtype=bool)

    mask = np.zeros(max_atoms, dtype=bool)

    central_chain = n['res_id'][2].tobytes()
    central_site = int(n['res_id'][3])

    get_mask(
        SIZE, MAX_SITES, chains, sites, names, 
        seen, removed, rand, mask, p, BACKBONE_ATOMS, bb_len,
        central_chain, central_site
    )

    return mask


def pad(arrays, max_atoms):
    """
    Returns LIST of ndarrays padded to max_atoms
    """
    return [pad_arr(arr, max_atoms) for arr in arrays]

def pad_arr(arr, padded_length):
    # get dtype of input array
    dt = arr.dtype

    # shape of sub arrays and first dimension (to be padded)
    shape = arr.shape[1:]
    orig_length = arr.shape[0]

    # check that the padding is large enough to accomdate the data
    if padded_length < orig_length:
        print('Error: Padded length of {}'.format(padded_length),
              'is smaller than original length of array {}'.format(orig_length))

    # create padded array
    padded_shape = (padded_length,*shape)
    mat_arr = np.zeros(padded_shape, dtype=dt)

    # add data to padded array
    mat_arr[:orig_length] = np.array(arr)

    return mat_arr

# @profile
def downsample(neighborhood, p, max_atoms=1000, remove_central=True):
    """
    Takes a neighborhood, removes p proportion
    of sidechains.
    
    Leaves backbone atoms in place. 
    
    Returns a copy of the neighborhood modifies
    
    ```
    # USAGE: 
    for neighborhood in data:
        neighborhood = downsample(neighborhood, .5)
        ...
    ```
    """
    mask = get_mask_c(neighborhood, p)

    info = [neighborhood['res_id'], 
    neighborhood['atom_names'][mask],
    neighborhood['elements'][mask],
    neighborhood['res_ids'][mask],
    neighborhood['coords'][mask],
    neighborhood['SASAs'][mask],
    neighborhood['charges'][mask]]
    
    info[1:] = pad(info[1:], max_atoms)
    
    x = np.zeros(shape=(1), dtype=dt)
    x[0] = (*info, )
    return x[0]



if __name__ == '__main__':
    with h5py.File(filename, "r") as f:
        data = np.unique(np.array(f['data'], dtype=dt), axis=0)

    p = .5 # proportion of sidechains to remove

    # Benchmark
    for n in tqdm(data):
        downsample(n, p)

    # Sanity check
    BACKBONE = np.array([b' N  ', b' CA ', b' C  ', b' O  '])
    def count_backbone_atoms(neighborhood):
        count = 0
        for a in neighborhood['atom_names']:
            if a in BACKBONE: count += 1
        return count

    def count_unique_sidechains(neighborhood):
        seen = set()
        for name, res_id in zip(neighborhood['atom_names'], neighborhood['res_ids']):
            if name not in BACKBONE: seen.add(tuple(res_id))
        return len(seen)

    def assert_central_removed(neighborhood):
        central = neighborhood['res_id']
        for name, res_id in zip(neighborhood['atom_names'], neighborhood['res_ids']):
            if name not in BACKBONE: 
                if(tuple(res_id) == tuple(central)): 
                    print(res_id)

    def report_confidence_interval(Ps, p, alpha=.95):
        import scipy.stats as st
        i = st.t.interval(alpha=alpha, df=len(Ps)-1,
              loc=np.mean(Ps),
              scale=st.sem(Ps))
        print(f'{100 * alpha}% confident that proportion of removed sidechains is in {i} (True p={p})')

    Ps = []
    for i in range(100, 200):
        n = data[i]
        n_1 = downsample(n, p)
        assert(count_backbone_atoms(n) == count_backbone_atoms(n_1))
        assert_central_removed(n_1)
        Ps.append(1 - (count_unique_sidechains(n_1) / (count_unique_sidechains(n) - 1))) # -1 because central is always removed

    print("\nSANITY CHECKS PASSED\n\t-All backbone atoms remain\n\t-central chain removed\n\t-", end="")
    report_confidence_interval(Ps, p)
