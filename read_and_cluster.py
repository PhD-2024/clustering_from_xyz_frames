import argparse
import matplotlib.pyplot as plt


def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {"true", "1", "yes"}:
        return True
    elif value.lower() in {"false", "0", "no"}:
        return False
    else:
        raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")

def read_xyz_file(filename):
    """
    Reads an XYZ file and returns the number of atoms, the atom types, and their coordinates.

    :param filename: Path to the XYZ file.
    :return: A tuple containing the number of atoms, a list of atom types, and a list of coordinates.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    num_atoms = int(lines[0].strip())

    atom_dict={}
    number_in_index_file=0

    for line in lines[2:num_atoms + 2]:
        parts = line.split()
        atom_types=parts[0]
        coordinates=[float(coord) for coord in parts[1:]]
        atom_dict[number_in_index_file] = (atom_types, coordinates)
        number_in_index_file += 1
    return atom_dict, num_atoms

def analyze_numbers(cluster_dict):
    """
    Analyzes the number of atoms in each cluster and returns a dictionary with the counts.

    :param cluster_dict: Dictionary where keys are cluster indices and values are lists of atom indices.
    :return: Dictionary with cluster sizes.
    """
    cluster_sizes = {}
    for key, atoms in cluster_dict.items():
        cluster_sizes[key] = len(atoms)
    return cluster_sizes

def get_distances(coord1, coord2):
    """currently without pbc, later replace this by pbc if necessary"""

    """
    Calculates the Euclidean distance between two points in 3D space.

    :param coord1: First coordinate as a list [x, y, z].
    :param coord2: Second coordinate as a list [x, y, z].
    :return: The Euclidean distance between the two coordinates.
    """
    return ((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (coord1[2] - coord2[2]) ** 2) ** 0.5

def connected_atoms(atom_dict, num_atoms, cutoff=1.65):
    """calculates all distances between atom_dict[i] and atom_dict[j] for i <j 
    Then compares those with the cutoff. If they are below the cutoff the indices
    [i,j] get added to the connected_atoms list

    """
    connected_atoms = []

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = get_distances(atom_dict[i][1], atom_dict[j][1])
            if distance < cutoff:
                connected_atoms.append([i, j])

    return connected_atoms

def from_connected_atoms_to_cluster(connected_atoms):
    """
    Converts a list of connected atoms into a dictionary of clusters.

    This function takes a list of connected atoms and organizes them into clusters.
    Each cluster is represented as a key-value pair in a dictionary, where the key
    is the cluster index (starting from 0) and the value is a list of atoms belonging
    to that cluster. The function ensures that all connected atoms are grouped together
    into the same cluster.

    Args:
        connected_atoms (list of lists): A list where each sublist contains indices
            of atoms that are directly connected to each other.
      

    Returns:
        dict: A dictionary where each key is a cluster index (int) and the value is a
        list of atom indices (list of int) that belong to that cluster.

    Example:
        >>> connected_atoms = [[0, 1], [1, 2], [3, 4]]
        >>> num_atoms = 5
        >>> from_connected_atoms_to_cluster(connected_atoms, num_atoms)
        {0: [0, 1, 2], 1: [3, 4]}

    Notes:
        - The function uses a helper function `merge_entries` to iteratively merge
          overlapping sublists in `connected_atoms` until no further merging is possible.
        - The merging process ensures that all atoms that are directly or indirectly
          connected are grouped into the same cluster.
    """

    cluster_dict = {}
    def merge_entries(list_of_lists):
        """merges lists with identical entries, stats with a list of lists"""

        for i in list_of_lists:
            for j in list_of_lists:
                if i !=j:
                    #check whether the two lists have common entries
                    if set(i) & set(j):
                        # merge the two lists
                        merged_entry= list(set(i) | set(j))
                        # remove the old entries
                        index_of_i = list_of_lists.index(i)
                        index_of_j = list_of_lists.index(j)
                        list_of_lists[index_of_i] = merged_entry
                        list_of_lists.pop(index_of_j)
                        return list_of_lists, True
        return list_of_lists, False
    not_yet_finished= True
    while not_yet_finished:
        connected_atoms, not_yet_finished = merge_entries(connected_atoms)
    for i in range(len(connected_atoms)):
        cluster_dict[i] = connected_atoms[i]
    return cluster_dict

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Read an XYZ file, calculate connected atoms, and cluster them.")
    parser.add_argument(
        "--filename",
        type=str,
        #default="tests/testfile_h.xyz",
        default="tests/Dimer-2metal_OptionA.xyz",
        help="Path to the XYZ file"
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.65,
        help="Cutoff distance for connecting atoms (default: 1.65)"
    )
    parser.add_argument(
        "--outname",
        type=str,
        default="tests/outputfile.txt",
        help="Name of the output file (default: outputfile.txt)"
    )
    parser.add_argument("--pass_to_Multiwfn", type=str_to_bool, default=True, help="If True, pass the clusters to Multiwfn (default: False)")
    parser.add_argument("--select_N_largest_clusters", type=int, default=2, help="If >0, select the N largest clusters (default: 2) for passing to Multiwfn")
    parser.add_argument("--indexing_1", type=str_to_bool, default=True, help="If True, output indices will be 1-indexed (default: True)")
    #parser.add_argument("--logname", type=str, default="spectrum.log", help="Name of the log file (default: spectrum.log)")
    #parser.add_argument("--fchk_name", type=str, default="spectrum.fchk", help="Name of the fchk file (default: spectrum.fchk)")
    parser.add_argument("--states", default=[1,2], help="List of states to analyze (default: [1, 2])")
    args = parser.parse_args()

    atoms, num_atoms = read_xyz_file(args.filename)
    connected_list = connected_atoms(atoms, num_atoms, cutoff=args.cutoff)
    cluster_dict = from_connected_atoms_to_cluster(connected_list)


    if args.indexing_1:
        print("Output will be 1-indexed.")
        # Convert indices to 1-indexed
        for key in cluster_dict:
            cluster_dict[key] = [index + 1 for index in cluster_dict[key]]
    else:  
        print("Output will be 0-indexed.")
    for key in cluster_dict:
        output_with_key=args.outname.replace(".txt", f"_{key}.txt")
        with open(output_with_key, 'w') as output_file:
           print(*cluster_dict[key], file=output_file, sep=",")



    clust_sizes=analyze_numbers(cluster_dict)
    print("Cluster sizes:", clust_sizes)
    #make histogram of cluster sizes
    plt.hist(clust_sizes.values(), bins=range(1, max(clust_sizes.values()) + 2), align='left', rwidth=0.8)
    plt.xlabel('Cluster Size')
    plt.ylabel('Frequency')
    plt.title('Histogram of Cluster Sizes')
    plt.savefig('cluster_sizes_histogram.pdf')
    #plt.show()

#select the N-largest clusters

if args.pass_to_Multiwfn and args.select_N_largest_clusters > 0:
    d={}
    counter=0
    sorted_clusters = sorted(clust_sizes.items(), key=lambda item: item[1], reverse=True)
    largest_clusters = sorted_clusters[:args.select_N_largest_clusters]


    for cluster_index, size in largest_clusters:

        output_with_key=args.outname.replace(".txt", f"_{cluster_index}.txt")
        with open(output_with_key, 'r') as output_file:
            atoms_in_cluster = output_file.read() #.strip()
            d[counter] = atoms_in_cluster.strip()
            counter += 1
           
    tmp_string=""
    for i in range(args.select_N_largest_clusters):
        tmp_string+=" '" + str(d[i]) + "'"

    print("suggested Multiwfn run:")

    print("-----------------------------------------")
    for state in args.states:
        print(f"for {args.select_N_largest_clusters} largest clusters and state {state}:")
        
        print(f"bash preliminary_bash_script.sh {args.select_N_largest_clusters} {state} {tmp_string} ")
        print("-----------------------------------------")
