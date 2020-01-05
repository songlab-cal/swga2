import src.optimize
import src.utility
import os
import h5py
import melting
import src.string_search
import multiprocessing
import numpy as np

def get_melting_tm(primer):
    return melting.temp(primer)

def get_gini(primer, fname_prefixes):
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

def get_gini_from_txt_for_one_k(task):
    primer_list, fname_prefix, fname_genome, seq_length, k = task
    rc_primer_list = [src.utility.reverse_complement(primer) for primer in primer_list]
    all_primer_list = list(set(primer_list + rc_primer_list))
    kmer_dict = src.string_search.get_all_positions_per_k(kmer_list=all_primer_list, seq_fname=fname_genome,fname_prefix=fname_prefix)
    ginis = []

    for primer in primer_list:
        position_diffs_forward = src.optimize.get_positional_gap_lengths(kmer_dict[primer], seq_length=seq_length)
        gini_forward = src.utility.gini_exact(position_diffs_forward)
        position_diffs_reverse = src.optimize.get_positional_gap_lengths(kmer_dict[src.utility.reverse_complement(primer)], seq_length=seq_length)
        gini_reverse = src.utility.gini_exact(position_diffs_reverse)
        ginis.append((gini_forward, gini_reverse))

    return dict(zip(primer_list, ginis))

def get_gini_from_txt(primer_list, fname_prefixes, fname_genomes, seq_lengths):
    tasks = []
    for i, fg_prefix in enumerate(fname_prefixes):
        for k in [6, 7, 8, 9, 10, 11, 12]:
            tasks.append(([primer for primer in primer_list if len(primer) == k], fg_prefix, fname_genomes[i], seq_lengths[i], k))

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = pool.map(get_gini_from_txt_for_one_k, tasks)

    primer_to_all_ginis = {}


    for primer_to_ginis in results:
        for primer, (gini_forward, gini_reverse) in primer_to_ginis.items():
            if primer not in primer_to_all_ginis:
                primer_to_all_ginis[primer] = []
            primer_to_all_ginis[primer].append(np.mean([gini_forward, gini_reverse]))
    return [np.mean(primer_to_all_ginis[primer]) for primer in primer_list]

def get_rate_from_h5py(primer, fname_prefixes):
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