import src.optimize
import src.utility
import os
import h5py
import melting
import src.string_search
import multiprocessing
import numpy as np

def get_melting_tm(primer):
    """
    Gets the predicted melting temperature of a primer using the melt package (https://pypi.org/project/melt/).

    Args:
        primer: The sequence of the primer.

    Returns:
        The predicted melting temperature.
    """
    return melting.temp(primer)

def get_gini(primer, fname_prefixes):
    """
    Get the gini index of the gaps between all adjacent positions a primer may bind. This computed from positional gaps
    both the forward and reverse strand.

    Args:
        primer: The sequence of the primer.
        fname_prefixes: The list of path prefixes to the h5py files.

    Returns:
        gini: The computed gini index of the positional gaps.

    """
    k = len(primer)
    positions_diffs = []
    for i, fname_prefix in enumerate(fname_prefixes):
        if os.path.exists(fname_prefix + '_' + str(k) + 'mer_positions.h5'):
            db=h5py.File(fname_prefix + '_' + str(k) + 'mer_positions.h5','r')
            if primer in db:
                position_diffs_forward = src.optimize.get_positional_gap_lengths(db[primer])
                positions_diffs.extend(position_diffs_forward)
            if src.utility.reverse_complement(primer) in db:
                position_diffs_reverse = src.optimize.get_positional_gap_lengths(db[src.utility.reverse_complement(primer)])
                positions_diffs.extend(position_diffs_reverse)
        else:
            print("Cannot find file: " + fname_prefix)
    gini = src.utility.gini(positions_diffs)
    return gini

# This probably belongs in a different file--its not a primer attribute.
def get_gini_from_txt_for_one_k(primer_list, fname_prefix, fname_genome, seq_length, circular):
    """
    Measures the gini index of the gaps between all adjacent positions any primer in primer_list may bind to.

    Args:
        primer_list: The list of primers to consider.
        fname_prefix: The path prefix to the h5py file.
        fname_genome: The path to the fasta file.
        seq_length: The length of the genome contained in fname_genome.

    Returns:
        primer_to_ginis: A dictionary of primers to a tuple of the gini indices (the first being computed from the
        forward strand, and the second being computed from the reverse strand).
    """
    # primer_list, fname_prefix, fname_genome, seq_length, k = task
    rc_primer_list = [src.utility.reverse_complement(primer) for primer in primer_list]
    all_primer_list = list(set(primer_list + rc_primer_list))
    kmer_dict = src.string_search.get_all_positions_per_k(kmer_list=all_primer_list, seq_fname=fname_genome, circular=circular, fname_prefix=fname_prefix)
    ginis = []

    for primer in primer_list:
        position_diffs_forward = src.optimize.get_positional_gap_lengths(kmer_dict[primer], circular, seq_length=seq_length)
        gini_forward = src.utility.gini_exact(position_diffs_forward)
        position_diffs_reverse = src.optimize.get_positional_gap_lengths(kmer_dict[src.utility.reverse_complement(primer)], circular, seq_length=seq_length)
        gini_reverse = src.utility.gini_exact(position_diffs_reverse)
        ginis.append((gini_forward, gini_reverse))

    primer_to_ginis = dict(zip(primer_list, ginis))
    return primer_to_ginis

def get_gini_from_txt_for_one_k_helper(args):
    primer_list, fname_prefix, fname_genome, seq_length, circular = args
    return get_gini_from_txt_for_one_k(primer_list, fname_prefix, fname_genome, seq_length, circular)

def get_gini_from_txt(primer_list, fname_prefixes, fname_genomes, seq_lengths, circular):
    """
    This runs get_gini_from_txt_for_one_k in a multiprocessed fashion where the task is divided based on the length
    of the primers.

    Args:
        primer_list: The list of primers to consider.
        fname_prefixes: The list of path prefixes to the h5py file.
        fname_genomes: The list of paths to the fasta files.
        seq_lengths: The length of all the genomes in fname_genomes (in the same order!).

    Returns:
        The average gini_index across all gini indices computed for each primer.
    """
    tasks = []
    for i, fg_prefix in enumerate(fname_prefixes):
        for k in [6, 7, 8, 9, 10, 11, 12]:
            primer_list_a = [primer for primer in primer_list if len(primer) == k]
            if len(primer_list_a) > 0:
                # print([primer_list_a, fg_prefix, fname_genomes[i], seq_lengths[i]])
                tasks.append([primer_list_a, fg_prefix, fname_genomes[i], seq_lengths[i], circular])

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    results = pool.map(get_gini_from_txt_for_one_k_helper, tasks)

    primer_to_all_ginis = {}

    for primer_to_ginis in results:
        for primer, (gini_forward, gini_reverse) in primer_to_ginis.items():
            if primer not in primer_to_all_ginis:
                primer_to_all_ginis[primer] = []
            primer_to_all_ginis[primer].append(np.mean([gini_forward, gini_reverse]))
    gini_mean = [np.mean(primer_to_all_ginis[primer]) for primer in primer_list]
    return gini_mean

def get_rate_from_h5py(primer, fname_prefixes):
    """
    Gets the frequency of a primer in both the forward and reverse strand.

    Args:
        primer: The sequence of the primer.
        fname_prefixes: The list of path prefixes to all the relevant h5py files.

    Returns:
        count: The frequency of a primer in both the forward and reverse strand.
    """
    k = len(primer)
    count = 0
    for i, fname_prefix in enumerate(fname_prefixes):
        if os.path.exists(fname_prefix + '_' + str(k) + 'mer_positions.h5'):
            db=h5py.File(fname_prefix + '_' + str(k) + 'mer_positions.h5','r')
            if primer in db:
                count += len(db[primer])
            if src.utility.reverse_complement(primer) in db:
                count += len(db[src.utility.reverse_complement(primer)])
        else:
            print("Cannot find file: " + fname_prefix)
    return count