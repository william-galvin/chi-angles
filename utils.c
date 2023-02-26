#include <stdio.h>


int in_backbone(int atom, int* backbone_atoms, int len) 
{
    for (int i = 0; i < len; i++) 
    {
        if (atom == backbone_atoms[i]) 
        {
            return 1;
        }
    }
    return 0;
}

/**
Params:
    - SIZE <- number of atoms (NOT padded length)
    - MAX_SITES <- assumed maximum number of amino acids in a protein
    - chains <- bytes of chains (b'A', b'B', etc)
    - sites <- int[] of site numbers
    - names <- Bytes of atom names (note: happen to be 32 bytes => int*)
    - seen, removed <- zeroed bool arrays of size MAX_CHAINS * MAX_SITES
    - random <- float[] (max atom) array on [0, 1]
    - mask <- zeroed bool array of size (max_atoms)
    - p <- proportion of sidechains to keep
    - backbone_atoms <- bytes of backbone atoms
    - bb_len <- number of backbone atoms
    - central_chain, central_site <- site and chain of central amino acid

Returns:
    - None
    - Modifies mask

Notes:
    - res_ids are unique within protein by 1) chain and 2) site. 
      By assuming there are no more than 10 chains and 10,000 sites, 
      we can created a collision-free "hash" map for `seen` and `removed`
      in 100,000 bytes each (0.1 MB)
*/
void get_mask(int SIZE, int MAX_SITES, 
              char* chains, long* sites, int* names, 
              char* seen, char* removed, double* random, char* mask, 
              float p, int* backbone_atoms, int bb_len, 
              char central_chain, long central_site)
{
    int central_index = (central_chain - 'A') * MAX_SITES + central_site;
    for (int i = 0; i < SIZE; i++) 
    {
        int id_index = (chains[i * 6] - 'A') * MAX_SITES + sites[i];
        if (! (seen[id_index]))
        {
            seen[id_index] = 1;
            removed[id_index] = random[i] < p || id_index == central_index;
        }
        mask[i] = in_backbone(names[i], backbone_atoms, bb_len) | (!removed[id_index]);
    }
}
