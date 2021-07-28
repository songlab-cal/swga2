import src.utility
import src.parameter
import multiprocessing
import h5py
import os

# #everything should be 5' to 3' written
def get_all_positions_per_k(kmer_list, seq_fname, circular, fname_prefix=None):
    """
    Gets all the positions of k-mers for one value of k. It assumes the first kmer in the list is the
    same length as all kmers in the list. Everything is written in the 5' to 3' direction.

    Args:
        kmer_list: The list of all kmers where k is a specific and single value.
        seq_fname: The path to the fasta file of the genome.
        fname_prefix: This argument is only necessary for knowing which fasta file is currently being processed.

    Returns:
        kmer_dict: Dictionary of kmer to count in the fasta file seq_fname.

    """
    if len(kmer_list) == 0:
        return {}
    k = len(kmer_list[0])
    kmer_dict = {}
    for kmer in kmer_list:
        kmer_dict[kmer] = []

    seq = src.utility.read_fasta_file(seq_fname)
    current_text = ''.join([next(seq) for _ in range(k)])
    beginning_text = current_text[:-1]

    if fname_prefix is not None:
        print("Starting the search for " + fname_prefix + ' ' + str(k) + 'mers...')

    i = 0
    if current_text in kmer_dict:
        kmer_dict[current_text].append(i)
        i += 1

    for c in seq:
        current_text = current_text[1:] + c.upper()
        if current_text in kmer_dict:
            kmer_dict[current_text].append(i)
        i += 1

    if circular:
        for c in beginning_text:
            current_text = current_text[1:] + c.upper()
            if current_text in kmer_dict:
                kmer_dict[current_text].append(i)
            i += 1
    return kmer_dict

def write_to_h5py(kmer_dict, fname_prefix):
    """
    Writes the kmer counts to an h5py file, which allows for efficient access in terms of looking up the
    frequency of a particular k-mer. If the kmer already exists in the dataset, the entry in the h5py file
    is overwritten with the new data.

    Args:
        kmer_dict: A dictionary of all the k-mers and their respective frequencies.
        fname_prefix: The file path prefix for the output h5py file. If k=6, for example, '_6mer_positions.h5'
        will be appended to the file name.
    """
    k = len(list(kmer_dict.keys())[0])
    f = h5py.File(fname_prefix + '_' + str(k) + 'mer_positions.h5', 'r+')
    for kmer, positions in kmer_dict.items():

        if kmer not in f:
            f.create_dataset(kmer, data=positions)
        else:
            del f[kmer]
            f.create_dataset(kmer, data=positions)
    f.close()

def check_which_primers_absent_in_h5py(primer_list, fname_prefix):
    """
    This function checks if the primers in a given list are missing from an h5py file.

    Args:
        primer_list: The list of primers that need to be checked exist in the h5py file.
        fname_prefix: The prefix of the h5py file--basically the path minus '_6mer_positions.h5' where k = 6.

    Returns:
        filtered_primer_list: All the primers from the given list of primers that are missing from the h5py file.
    """
    if len(primer_list) == 0:
        return primer_list

    k = len(primer_list[0])

    if not os.path.exists(fname_prefix + '_' + str(k) + 'mer_positions.h5'):
        f = h5py.File(fname_prefix + '_' + str(k) + 'mer_positions.h5', 'a')
        f.close()
        return primer_list

    # print("Printing in check for " + fname_prefix + ' ' + str(k) + 'mer with ' + str(len(primer_list)) + " primers")

    all_present_kmers_in_genome = []
    with open(fname_prefix + '_' + str(k) + 'mer_all.txt', 'r') as txt_f:
        for line in txt_f:
            curr_kmer = line.split(" ")[0]
            all_present_kmers_in_genome.append(curr_kmer)
    all_present_kmers_in_genome = set(all_present_kmers_in_genome)

    # print('READING ' + fname_prefix + '_' + str(k) + 'mer_positions.h5')
    f = h5py.File(fname_prefix + '_' + str(k) + 'mer_positions.h5', 'r')
    keys = set(f.keys())
    f.close()

    filtered_primer_list = []
    for primer in primer_list:
        if primer not in keys:
            if primer in all_present_kmers_in_genome:
                filtered_primer_list.append(primer)
    return filtered_primer_list

def get_positions(primer_list, fname_prefixes, fname_genomes, circular, overwrite=False, no_all_primer_files=False):
    """
    Launches a multiprocessing pool to check if all primers exists in their relevant h5py file and modifies the file
    if frequencies for that k-mer is missing.

    Args:
        primer_list: List of k-mers to be checked that exist in the h5py files or to modify them if not.
        fname_prefixes: The path prefixes for the h5py files, basically the path minus '_6mer_positions.h5' where k = 6.
        fname_genomes: A list of paths to the fasta files.
        overwrite: Boolean which when set to true means overwrite the k-mer entries in the h5py file if it already exists.
    """
    tasks = []
    for i, fg_prefix in enumerate(fname_prefixes):
        for k in [6,7,8,9,10,11,12]:
            tasks.append(([primer for primer in primer_list if len(primer) == k], fg_prefix, fname_genomes[i], k, circular, overwrite))
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(append_positions_to_h5py_file, tasks)

def append_positions_to_h5py_file(task):
    """
    Check if the primer exists in the h5py file and modify all
    the file if frequencies for that k-mer is missing.

    Args:
        task: A tuple consisting of the following arguments:
            primer_list: List of k-mers to be checked that exist in the h5py files or to modify them if not.
            fname_prefix: The path prefix for the h5py files, basically the path minus '_6mer_positions.h5' where k = 6.
            fname_genome: A list of paths to the fasta files.
            k: The length of the k-mers.
            overwrite: Boolean which when set to true means overwrite the k-mer entries in the h5py file if it already exists.

    """
    primer_list, fname_prefix, fname_genome, k, circular, overwrite = task

    if len(primer_list) == 0:
        return

    new_list = set(primer_list)
    new_list.update([src.utility.reverse_complement(primer) for primer in new_list])
    one_k_list = [primer for primer in sorted(list(new_list)) if len(primer) == k]

    if not overwrite:
        filtered_list = check_which_primers_absent_in_h5py(one_k_list, fname_prefix)
    else:
        filtered_list = one_k_list

    if len(filtered_list) > 0:
        kmer_dict = get_all_positions_per_k(list(set(filtered_list)), fname_genome, circular, fname_prefix)
        write_to_h5py(kmer_dict, fname_prefix)

if __name__ == "__main__":
    k = 6
    # output_prefix = 'kmer_files/myco'

    # test_txt = src.parameter.data_dir + output_prefix + '_' + str(k) + 'mer.txt'

    kmer_dict = get_all_positions_per_k(["AACC"], seq_fname='/Users/janeyu/Documents/primer_stuff/soapswga/examples/simple_plasmid_example/pcDNA.fasta', circular=False)
    print(len(kmer_dict["AACC"]))