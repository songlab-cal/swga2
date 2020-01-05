import src.parameter
import src.primer_attributes
from collections import Counter
import pickle
import pandas as pd
import melting
import multiprocessing
import src.dimer

def filter_extra(primer, verbose=False):
    #1 Avoid using three G or C nucleotides in a row at the 3′-end of a primer.
    #2 The GC content (the number of G's and C's in the primer as a percentage of the total bases) of primer should be 40-60%.
    #3 The presence of G or C bases within the last five bases from the 3' end of primers (GC clamp) helps promote specific binding at the 3' end due to the stronger bonding of G and C bases. More than 3 G's or C's should be avoided in the last 5 bases at the 3' end of the primer.
    #4 A repeat is a di-nucleotide occurring many times consecutively and should be avoided because they can misprime. For example: ATATATAT. A maximum number of di-nucleotide repeats acceptable in an oligo is 4 di-nucleotides.
    #5 Primers with long runs of a single base should generally be avoided as they can misprime.For example, AGCGGGGGATGGGG has runs of base 'G' of value 5 and 4. A maximum number of runs accepted is 4bp.

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

    if verbose:
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

    if verbose:
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

    if verbose:
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

    if verbose:
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

    if dimer.is_dimer(primer, primer, parameter.default_max_self_dimer_bp):
        # print("here")
        return False

    return True

#fg_rate and bg_rate should be the number of binding sites scaled by the genome length
def get_all_rates(primer_list, fg_prefixes, bg_prefixes, fg_total_length, bg_total_length, output_df_fname=src.parameter.data_dir + 'primer_candidate_list_myco_human.p'):

    primer_to_fg_count = get_rates_for_one_species(primer_list, fg_prefixes)
    primer_to_bg_count = get_rates_for_one_species(primer_list, bg_prefixes)

    results = []

    for primer in primer_list:
        fg_count = primer_to_fg_count[primer]
        fg_bool = (fg_count is None or fg_count/fg_total_length > parameter.min_fg_freq)

        bg_count = primer_to_bg_count[primer]
        bg_bool = (bg_count is None or bg_count/bg_total_length < parameter.max_bg_freq)
        results.append([primer, fg_count, bg_count, fg_bool, bg_bool])

    df = pd.DataFrame(results, columns=['primer', 'fg_count', 'bg_count', 'fg_bool', 'bg_bool'])
    if output_df_fname:
        pickle.dump(df, open(output_df_fname, 'wb'))
    return df

def get_rates_for_one_species(primer_list, fname_prefixes):
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
    results = pool.map(get_rate_for_one_file, tasks)

    all_primer_to_count = {}

    for primer_to_count in results:
        for primer, count in primer_to_count.items():
            if primer not in all_primer_to_count:
                all_primer_to_count[primer] = count
            else:
                all_primer_to_count[primer] += count
    return all_primer_to_count

def get_rate_for_one_file(task):
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

def get_gini(fg_prefixes, fg_genomes, fg_seq_lengths, df=None, input_df_fname=src.parameter.data_dir + 'primer_candidate_list_myco_human.p'):
    if df is None:
        df = pickle.load(open(input_df_fname, 'rb'))
    df['gini'] = primer_attributes.get_gini_from_txt(df['primer'].values, fg_prefixes, fg_genomes, fg_seq_lengths)
    print(df['gini'])

    if len(df['gini']) == 0:
        df['gini_bool'] = []
        return df

    df['gini_bool'] = df.apply(lambda x: x['gini'] is not None and x['gini'] < src.parameter.max_gini, axis=1)

    return df

if __name__ == "__main__":
    human_chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
                      'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM',
                      'chrX', 'chrY']
    primer_list = ['AACGATGAC', 'AAGAACGAC', 'AAGATGTCG', 'AATGACGAC']

    myco_seq_lengths = [4411532]
    human_seq_lengths = [parameter.seq_len['human_' + chr] for chr in human_chr_list]
    myco_prefixes = [parameter.data_dir + 'kmer_files/myco']
    human_prefixes = [parameter.data_dir + 'kmer_files/human_' + chr for chr in human_chr_list]
    human_genomes = [parameter.data_dir + 'genomes/' + chr + '.fa' for chr in human_chr_list]
    myco_genomes = [parameter.data_dir + 'genomes/MTBH37RV.fasta']

    primer_to_fg_count = get_rates_for_one_species(['CTAGGTAGTAG'], myco_prefixes)
    primer_to_bg_count = get_rates_for_one_species(['CTAGGTAGTAG'], human_prefixes)
    print(primer_to_fg_count)
    print(primer_to_bg_count)

    print(get_all_rates(['CTAGGTAGTAG'], myco_prefixes, human_prefixes, sum(myco_seq_lengths), sum(human_seq_lengths)))

    print(filter_extra('TTCGTACCG', verbose=True))
