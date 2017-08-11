#!/usr/bin/env python
"""
Merges multiple mappings into a single file.

Hazen 08/17
"""
import pickle

def mergeMaps(mapping_files):
    """
    mapping_files is a list of mapping file names. 

    The result is a dictionary containing all the mappings. 
    Channels are converted to the order in the mapping_files 
    list. It is assumed that channel 0 is the same in all
    of the files.
    """
    assert(len(mapping_files) > 0)
    
    mappings = {}
    for i, filename in enumerate(mapping_files):

        # Load channel to channel mapping.
        with open(filename, 'rb') as fp:
            mp_transform = pickle.load(fp)

        # Add to overall mappings.
        mappings["0_" + str(i+1) + "_x"] = mp_transform["0_1_x"]
        mappings["0_" + str(i+1) + "_y"] = mp_transform["0_1_y"]
        mappings[str(i+1) + "_0_x"] = mp_transform["1_0_x"]
        mappings[str(i+1) + "_0_y"] = mp_transform["1_0_y"]

    return mappings


if (__name__ == "__main__"):
    import argparse

    parser = argparse.ArgumentParser(description = 'Merge mappings.')

    parser.add_argument('--results', dest='results', type=str, required=True,
                        help = "The name of the file to save the merged mappings in.")
    parser.add_argument('--maps', dest='maps', type=str, required=True, nargs = '*',
                        help = "The name of the mapping files to merge.")

    args = parser.parse_args()

    merged_map = mergeMaps(args.maps)

    with open(args.results, 'wb') as fp:
        pickle.dump(merged_map, fp)
        
    if False:
        for key in sorted(merged_map):
            print(key)

    
