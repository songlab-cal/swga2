from optparse import OptionParser
import parameter
import multiprocessing
import os
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

parser = OptionParser()
parser.add_option('-k', '--kmer_fore', help='')
parser.add_option('-l', '--kmer_back', help='')
parser.add_option('-x', '--fasta_fore', help='')
parser.add_option('-y', '--fasta_back', help='')
parser.add_option('-j', '--json-file', type=str, help='')
parser.add_option('-z', '--data_dir', default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data/project'), type=str, help='')
parser.add_option('-u', '--min_fg_freq', default=float(1/100000), type=float, help='')
parser.add_option('-v', '--max_bg_freq', default=float(1/200000), type=float, help='')
parser.add_option('g', '--max_gini', default=0.6, type=float, help='')
parser.add_option('--max_primer', default=500, type=int, help='')
parser.add_option('--min_amp_pred', default=5, type=float, help='')
parser.add_option('-t', '--min_tm', default=15, type=int, help='')
parser.add_option('-u', '--max_tm', default=45, type=int, help='')
parser.add_option('--max_dimer_bp', default=3, type=int, help='')
parser.add_option('--max_self_dimer_bp', default=4, type=int, help='')
parser.add_option('--selection_metric', default='deterministic', type=str, help='')
parser.add_option('--iterations', default=8, type=int, help='')
parser.add_option('--max_sets', default=5, type=int, help='')
parser.add_option('--fg_circular', action="store_true", help='')
parser.add_option('--bg_circular', action="store_true", help='')
parser.add_option('--drop_iterations', default=[4], type=int, help='')
parser.add_option('-c', '--cpus', default=int(multiprocessing.cpu_count()), type=int, help='')
(options, args) = parser.parse_args()
print(options)
params = parameter.write_args_to_json(options)
print(params)

import kmer
import rf_preprocessing
import utility
import pickle
import optimize
import numpy as np
import filter
import string_search
import time
import os

fg_prefixes = params['fg_prefixes']
bg_prefixes = params['bg_prefixes']

fg_genomes = params['fg_genomes']
bg_genomes = params['bg_genomes']

fg_seq_lengths = params['fg_seq_lengths']
bg_seq_lengths = params['bg_seq_lengths']

def step1():
    """
        Creates files of all k-mers of length 6 to 12 which is located in the path specificed by --kmer_fore and --kmer_b .
    """
    for i, fg_prefix in enumerate(fg_prefixes):
        kmer.run_jellyfish(fg_genomes[i], fg_prefix)

    for i, bg_prefix in enumerate(bg_prefixes):
        kmer.run_jellyfish(bg_genomes[i], bg_prefix)

    return fg_seq_lengths, bg_seq_lengths

def step2(all_primers=None):
    """Flattens a multidimensional array.

    Args:
        l: The list to flatten.

    Returns:
        flattened_list: The flattened array.
    """
    kwargs = {'fg_prefixes': fg_prefixes, 'bg_prefixes': bg_prefixes, 'fg_total_length': sum(fg_seq_lengths),
              'bg_total_length': sum(bg_seq_lengths)}
    # print(kwargs)
    if all_primers is None:
        all_primers = kmer.get_primer_list_from_kmers(fg_prefixes)

    # all_primers = all_primers[:100]

    print("Retreiving foreground and background rates.")
    rate_df = filter.get_all_rates(all_primers, **kwargs)
    filtered_rate_df = rate_df[(rate_df['fg_bool']) & (rate_df['bg_bool'])]
    filtered_rate_df = filtered_rate_df.drop(columns=['fg_bool', 'bg_bool'])
    print(rate_df)
    gini_df = filter.get_gini(fg_prefixes, fg_genomes, fg_seq_lengths, df = filtered_rate_df)
    print("Filtered " + str(len(filtered_rate_df) - len(gini_df)) + " number of primers based on gini rate.")

    gini_df['ratio'] = gini_df['bg_count']/gini_df['fg_count']
    filtered_gini_df = gini_df.sort_values(by=['ratio'], ascending=False)[:parameter.max_primer]

    pickle.dump(filtered_gini_df, open(os.path.join(parameter.data_dir, 'step2_df.p'), 'wb'))
    string_search.get_positions(filtered_gini_df['primer'], fg_prefixes, fg_genomes)
    string_search.get_positions(filtered_gini_df['primer'], bg_prefixes, bg_genomes)
    return filtered_gini_df

#RANK BY RANDOM FOREST
def step3(out_file_prefix=None):
    print(parameter.src_dir)
    step2_df = pickle.load(open(os.path.join(parameter.data_dir, 'step2_df.p'), 'rb'))

    print("Length before filtering: " + str(step2_df.shape[0]))
    primer_list = step2_df['primer']
    fg_scale = sum(fg_seq_lengths)/6200

    df_pred = rf_preprocessing.create_augmented_df(fg_prefixes, primer_list)
    df_pred = rf_preprocessing.scale_delta_Gs(df_pred, on_scale=fg_scale)
    df_pred['molarity'] = 2.5

    results = rf_preprocessing.predict_new_primers(df_pred)
    results.sort_values(by=['on.target.pred'], ascending=[False], inplace=True)

    results['on.target.pred'] = results['on.target.pred']
    step3_df = results[results['on.target.pred'] >= parameter.min_amp_pred]

    step2_df = step2_df.set_index('primer')
    step3_df = step3_df.set_index('sequence')

    joined_step3_df = step3_df.join(step2_df[['ratio', 'gini', 'fg_count', 'bg_count']], how='left').sort_values(by='gini')

    if out_file_prefix:
        pickle.dump(joined_step3_df, open(out_file_prefix + "df_pred.p", "wb"))
        joined_step3_df.to_csv(out_file_prefix+'_regression_results.csv')
    else:
        pickle.dump(joined_step3_df, open(os.path.join(parameter.data_dir, 'step3_df.p'), "wb"))

    print("Length after filtering: " + str(joined_step3_df.shape[0]))
    print(joined_step3_df)

    return joined_step3_df

# OPTIMIZE AND SEARCH
def step4(primer_list=None, scores = None, target_var='coverage', selection_metric='deterministic', initial_primer_sets=None, iterations=10, max_sets=10, in_file_prefix=None, out_file=None, fg_circular=True, bg_circular=False, banned_primers=[]):
    if primer_list is None:
        step3_df = pickle.load(open(os.path.join(parameter.data_dir, 'step3_df.p'), 'rb'))
        primer_list = step3_df.index
        scores = step3_df['on.target.pred']
    if type(primer_list) is str:
        primer_list = pickle.load(open(primer_list, 'rb'))['sequence']

        # if scores is None:


    all_results = []
    all_scores = []
    cache = {}

    drop_indices_per_level = [[4], [7], [6], [7], [6], [7], [5], [6], [6], [6], [6]]

    for i in range(1):
        ts = time.time()
        print(ts)
        print("Repeat #: " + str(i+1))

        if initial_primer_sets is not None:
            primer_list.extend(utility.flatten(initial_primer_sets))
        else:
            initial_primer_sets = optimize.random_intial_start(np.asarray(primer_list).reshape((-1, 1)), scores, max_sets=max_sets)

        print("Initial primers: ")
        print(initial_primer_sets)
        print(cache)

        primer_list = sorted(list(set(primer_list)))
        random_primers = utility.flatten(initial_primer_sets)
        combined_primer_list = list(set(primer_list+random_primers))

        print("Banned primers: " + str(banned_primers))
        results, scores, cache = optimize.bfs(combined_primer_list, fg_prefixes, bg_prefixes, fg_seq_lengths, bg_seq_lengths, initial_primer_sets=initial_primer_sets, iterations=iterations, max_sets=max_sets, target_var=target_var, selection_method=selection_metric, banned_primers=banned_primers, fg_circular=fg_circular, bg_circular=bg_circular, cache=cache, drop_indices=drop_indices_per_level[i])
        print(results)
        banned_primers = list(set(utility.flatten(results)))

        all_results.append(results)
        all_scores.append(scores)

        banned_primers=[utility.most_frequent(utility.flatten(results)), utility.most_frequent(utility.flatten(results))]

    if out_file:
        pickle.dump(results, open(out_file, 'wb'))

    all_results = utility.flatten(all_results)
    all_scores = utility.flatten(all_scores)

    results = []
    scores = []

    for i, result in enumerate(all_results):
        if result not in result:
            results.append(result)
            scores.append(all_scores[i])

    idx = np.argsort(scores)[::-1]
    print(list(np.array(results)[idx.astype(int)]))
    print(list(np.array(scores)[idx.astype(int)]))

    return np.array(results)[idx.astype(int)], np.array(scores)[idx.astype(int)]

def step5(primer_sets):

    primer_scores = optimize.stand_alone_score_primer_sets(primer_sets)

    argsort_indices = np.argsort(primer_scores)[::-1]

    selection_indices = [0]

    avoid_indices = []

    for index in argsort_indices:
        can_add = True
        for selected_index in selection_indices:
            num_intersect = len(utility.intersection(primer_sets[selected_index], primer_sets[index]))
            if num_intersect > 4:
                can_add = False
        if can_add and index not in avoid_indices:
            selection_indices.append(index)

    selected_sets = [primer_sets[selected_index] for selected_index in selection_indices]
    selected_scores = [primer_scores[selected_index] for selected_index in selection_indices]

    for i, selected_set in enumerate(selected_sets):
        print(str(selected_set) + ", score=%0.5f" % selected_scores[i])

if __name__ == "__main__":
    print("pipeline.py")

    # global params


    # step1()
    step2()
    # step3()
    # step4()
    # print(parameter.min_fg_freq)

