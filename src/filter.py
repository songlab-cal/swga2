import src.parameter
import src.primer_attributes
from collections import Counter
import pickle
import pandas as pd
import melting
import multiprocessing
import src.dimer

def filter_extra(primer):
    """Filters a primer based on the following five rules:
    1. Avoid using three G or C nucleotides in a row at the 3′-end of a primer.
    2. The GC content (the number of G's and C's in the primer as a percentage of the total bases) of primer should be 40-60%.
    3. The presence of G or C bases within the last five bases from the 3' end of primers (GC clamp) helps promote specific binding at the 3' end due to the stronger bonding of G and C bases. More than 3 G's or C's should be avoided in the last 5 bases at the 3' end of the primer.
    4. A repeat is a di-nucleotide occurring many times consecutively and should be avoided because they can misprime. For example: ATATATAT. A maximum number of di-nucleotide repeats acceptable in an oligo is 4 di-nucleotides.
    5. Primers with long runs of a single base should generally be avoided as they can misprime.For example, AGCGGGGGATGGGG has runs of base 'G' of value 5 and 4. A maximum number of runs accepted is 4bp.

    Args:
        primer: The primer to be evaluated, written in the 5' to 3' direction.

    Returns:
        filter_bool: True if it passes all filters; false otherwise.
    """

    primer_tm = melting.temp(primer)
    # print(primer_tm)
    if primer_tm > 45 or primer_tm < 15:
        if verbose:
            print("Melting temperature: " + str(primer_tm))
        return False

    if verbose:
        print("Checking 5 in a row.")

    #Rule 5
    if 'GGGGG' in primer:
        return False
    if 'CCCCC' in primer:
        return False
    if 'AAAAA' in primer:
        return False
    if 'TTTTT' in primer:
        return False

    if src.parameter.verbose:
        print("Checking the GC content")

    #Rule 2
    all_cnt = Counter(primer)
    # print(all_cnt)
    if 'C' not in all_cnt:
        all_cnt['C'] = 0
    if 'G' not in all_cnt:
        all_cnt['G'] = 0
    GC_content = (all_cnt['G'] + all_cnt['C'])/float(len(primer))
    if GC_content <= 0.375 or GC_content >= 0.625:
        return False

    if src.parameter.verbose:
        print("Checking last five bases")

    #Rule 3
    last_five_count = Counter(primer[-5:])
    # print(last_five_count)
    if 'C' not in last_five_count:
        last_five_count['C'] = 0
    if 'G' not in last_five_count:
        last_five_count['G'] = 0
    if last_five_count['G'] + last_five_count['C'] > 3 or last_five_count['G'] + last_five_count['C'] == 0:
        return False

    if src.parameter.verbose:
        print("Checking last three bases")

    # Rule 1
    last_three_count = Counter(primer[-3:])
    # print(last_three_count)
    if 'C' not in last_three_count:
        last_three_count['C'] = 0
    if 'G' not in last_three_count:
        last_three_count['G'] = 0
    if last_three_count['G'] + last_three_count['C'] == 3:
        return False

    if src.parameter.verbose:
        print("Checking dinucleotide repeats")

    #Rule 4
    if len(primer) >= 10:
        more_than_5 = [nucleo for nucleo, count in all_cnt.items() if count >= 5]
        # print(more_than_5)
        if len(more_than_5) > 1:
            if 'A' in more_than_5:
                if 'T' in more_than_5:
                    if 'ATATATATAT' in primer or 'TATATATATA' in primer:
                        return False
                if 'G' in more_than_5:
                    if 'AGAGAGAGAG' in primer or 'GAGAGAGAGA' in primer:
                        return False
                if 'C' in more_than_5:
                    if 'ACACACACAC' in primer or 'CACACACACA'in primer:
                        return False
            if 'T' in more_than_5:
                if 'C' in more_than_5:
                    if 'TCTCTCTCTC' in primer or 'CTCTCTCTCT' in primer:
                        return False
                if 'G' in more_than_5:
                    if 'GTGTGTGTGT' in primer or 'TGTGTGTGTG' in primer:
                        return False
            if 'C' in more_than_5:
                if 'G' in more_than_5:
                    if 'CGCGCGCGCG' in primer or 'GCGCGCGCGC' in primer:
                        return False

    if src.dimer.is_dimer(primer, primer, src.parameter.default_max_self_dimer_bp):
        # print("here")
        return False

    return True


def get_all_rates(primer_list, fg_prefixes, bg_prefixes, fg_total_length, bg_total_length):
    """
    Computes the foreground and background binding site frequencies normalized by their respective genome lengths.

    Args:
        primer_list: The list of primers to compute frequencies for.
        fg_prefixes: The list of foreground path prefixes used for creating the kmer files.
        bg_prefixes: The list of background path prefixes used for creating the kmer files.
        fg_total_length: The total number of base pairs in the foregound genome.
        bg_total_length: The total number of base pairs in the background genome.

    Returns:
        df: A pandas dataframe with the sequence, unnormalized counts, and  columns fg_bool and bg_bool which indicate if the sequence passes the respective filters.
    """

    primer_to_fg_count = get_rates_for_one_species(primer_list, fg_prefixes)
    primer_to_bg_count = get_rates_for_one_species(primer_list, bg_prefixes)

    # print(primer_to_fg_count)
    # print(primer_to_bg_count)

    results = []

    for primer in primer_list:
        fg_count = primer_to_fg_count[primer]
        fg_bool = (fg_count is None or fg_count / fg_total_length > src.parameter.min_fg_freq)
        bg_count = primer_to_bg_count[primer]
        bg_bool = (bg_count is None or bg_count / bg_total_length < src.parameter.max_bg_freq)
        results.append([primer, fg_count, bg_count, fg_bool, bg_bool])

    df = pd.DataFrame(results, columns=['primer', 'fg_count', 'bg_count', 'fg_bool', 'bg_bool'])

    return df


def get_rates_for_one_species(primer_list, fname_prefixes):
    """
    Computes the binding site frequencies for all ppsth prefixes in fname_prefixes.

    Args:
        primer_list: The list of primers to compute frequencies for.
        fg_prefixes: The list of foreground path prefixes used for creating the kmer files.

    Returns:
        all_primer_to_count: A dictonary of primer to frequency.
    """
    stratified_primer_list = {}

    for primer in primer_list:
        k = len(primer)
        if k not in stratified_primer_list:
            stratified_primer_list[k] = []
        stratified_primer_list[k].append(primer)

    tasks = []

    for fname_prefix in fname_prefixes:
        for k, primer_list_k in stratified_primer_list.items():
            tasks.append((primer_list_k, fname_prefix, k))
            # get_rate_for_one_file((primer_list_k, fname_prefix, k))

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = pool.map(_get_rate_for_one_file, tasks)

    all_primer_to_count = {}

    for primer_to_count in results:
        for primer, count in primer_to_count.items():
            if primer not in all_primer_to_count:
                all_primer_to_count[primer] = count
            else:
                all_primer_to_count[primer] += count
    return all_primer_to_count


def _get_rate_for_one_file(task):
    primer_list, fname_prefix, k = task
    primer_to_count = {}
    with open(fname_prefix + '_' + str(k) + 'mer_all.txt', 'r') as f_in:
        for line in f_in:
            primer = line.split(" ")[0]
            count = line.split(" ")[1]
            primer_to_count[primer] = int(count)

    all_counts = []
    for primer in primer_list:
        if primer in primer_to_count:
            all_counts.append(primer_to_count[primer])
        else:
            all_counts.append(0)
    return dict(zip(primer_list, all_counts))


def get_gini(fg_prefixes, fg_genomes, fg_seq_lengths, df, circular):
    """Computes the Gini index of the gap distances between binding sites.

    Args:
        fg_prefixes: List of path prefixes to the kmer files of the foreground genome.
        fg_genomes: List of paths to the foreground fasta files.
        fg_seq_lengths: List of sequence length(s) of the foreground genome(s).
        df: Pandas dataframe with column primer containing the primer sequences.

    Returns:
        df: Input dataframe with new column 'gini' for the computed Gini indices.

    """
    df['gini'] = src.primer_attributes.get_gini_from_txt(df['primer'].values, fg_prefixes, fg_genomes, fg_seq_lengths, circular)

    if len(df['gini']) == 0:
        df['gini_bool'] = []
        return df

    df['gini_bool'] = df.apply(lambda x: x['gini'] is not None and x['gini'] < src.parameter.max_gini, axis=1)

    return df[df['gini_bool']]
