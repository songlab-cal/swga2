import src.utility
import src.parameter
import multiprocessing
import h5py
import os

# #everything should be 5' to 3' written
def get_all_positions_per_k(kmer_list, seq_fname=None, fname_prefix=None):
    """


    Args:

    Returns:
    """
    if len(kmer_list) == 0:
        return {}
    k = len(kmer_list[0])
    kmer_dict = {}
    for kmer in kmer_list:
        kmer_dict[kmer] = []

    seq = src.utility.read_fasta_file(seq_fname)
    current_text = ''.join([next(seq) for _ in range(k)])

    if fname_prefix is not None:
        print("Starting the search for " + fname_prefix + ' ' + str(k) + 'mers...')

    i = 0
    for c in seq:
        if current_text in kmer_dict:
            kmer_dict[current_text].append(i)
        current_text = current_text[1:] + c.upper()
        i += 1
    return kmer_dict

def write_to_h5py(kmer_dict, fname_prefix):
    k = len(list(kmer_dict.keys())[0])
    # print('APPEND ' + fname_prefix + '_' + str(k) + 'mer_positions.h5')
    f = h5py.File(fname_prefix + '_' + str(k) + 'mer_positions.h5', 'r+')
    for kmer, positions in kmer_dict.items():

        if kmer not in f:
            f.create_dataset(kmer, data=positions)
        else:
            del f[kmer]
            f.create_dataset(kmer, data=positions)
    f.close()

def check_which_primers_absent_in_h5py(primer_list, fname_prefix):
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

def get_positions(primer_list, fname_prefixes, fname_genomes, overwrite=False, no_all_primer_files=False):
    tasks = []
    for i, fg_prefix in enumerate(fname_prefixes):
        for k in [6,7,8,9,10,11,12]:
            tasks.append(([primer for primer in primer_list if len(primer) == k], fg_prefix, fname_genomes[i], k, overwrite, no_all_primer_files))
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    pool.map(append_positions_to_h5py_file, tasks)

def append_positions_to_h5py_file(task):
    primer_list, fname_prefix, fname_genome, k, overwrite, no_all_primer_files = task

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
        kmer_dict = get_all_positions_per_k(list(set(filtered_list)), fname_genome, fname_prefix)
        write_to_h5py(kmer_dict, fname_prefix)

if __name__ == "__main__":
    k = 6
    output_prefix = 'kmer_files/myco'

    test_txt = src.parameter.data_dir + output_prefix + '_' + str(k) + 'mer.txt'

    kmer_dict = get_all_positions_per_k(["CCGAATCG"], seq_fname=src.parameter.data_dir + 'genomes/chr2.fa')
    print(len(kmer_dict["CCGAATCG"]))