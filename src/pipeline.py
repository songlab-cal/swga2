from optparse import OptionParser
import src.parameter
import src.utility
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
parser.add_option('-p', '--min_fg_freq', default=float(1/100000), type=float, help='')
parser.add_option('-q', '--max_bg_freq', default=float(1/200000), type=float, help='')
parser.add_option('-g', '--max_gini', default=0.6, type=float, help='')
parser.add_option('-e', '--max_primer', default=500, type=int, help='')
parser.add_option('-a', '--min_amp_pred', default=5, type=float, help='')
parser.add_option('-m', '--min_tm', default=15, type=int, help='')
parser.add_option('-n', '--max_tm', default=45, type=int, help='')
parser.add_option('-t', '--max_dimer_bp', default=3, type=int, help='')
parser.add_option('-u', '--max_self_dimer_bp', default=4, type=int, help='')
parser.add_option('-s', '--selection_metric', default='deterministic', type=str, help='')
parser.add_option('-i', '--iterations', default=8, type=int, help='')
parser.add_option('-#', '--top_set_count', default=10, type=int, help='')
parser.add_option('-r', '--retries', default=5, type=int, help='')
parser.add_option('-w', '--max_sets', default=5, type=int, help='')
parser.add_option('-f', '--fg_circular', default=True, help='')
parser.add_option('-b', '--bg_circular', default=False, help='')
parser.add_option('-d', '--drop_iterations', default=[5], type=int, help='')
parser.add_option('-v', '--verbose', default=False, help='')
parser.add_option('-c', '--cpus', default=int(multiprocessing.cpu_count()), type=int, help='')
(options, args) = parser.parse_args()
params = src.parameter.write_args_to_json(options)

print('======================================================================================================')
print("Parameters:")
for param_key, param_val in params.items():
    print(param_key + ": " + str(param_val))
print('======================================================================================================')

import src.kmer
import src.rf_preprocessing
import src.utility
import pickle
import src.optimize
import numpy as np
import src.filter
import src.string_search
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
    for prefix in fg_prefixes + bg_prefixes:
        if not os.path.exists(os.path.dirname(prefix)):
            os.makedirs(os.path.dirname(prefix))

    print("Running jellyfish for foreground...")
    for i, fg_prefix in enumerate(fg_prefixes):
        src.kmer.run_jellyfish(fg_genomes[i], fg_prefix)

    print("Running jellyfish for background...")
    for i, bg_prefix in enumerate(bg_prefixes):
        src.kmer.run_jellyfish(bg_genomes[i], bg_prefix)

    print('Done running jellyfish')

    return fg_seq_lengths, bg_seq_lengths

def step2(all_primers=None):
    """Filters all candidate primers according to primer design principles (http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html)
    and minimum foreground frequency, maximum background frequency, and maximum Gini index of distances between binding sites.

    Args:
        all_primers: The list of candidate primers to consider.

    Returns:
        filtered_gini_df: .
    """
    kwargs = {'fg_prefixes': fg_prefixes, 'bg_prefixes': bg_prefixes, 'fg_total_length': sum(fg_seq_lengths),
              'bg_total_length': sum(bg_seq_lengths)}
    if all_primers is None:
        all_primers = src.kmer.get_primer_list_from_kmers(fg_prefixes)

    print("Computing foreground and background rates...")
    rate_df = src.filter.get_all_rates(all_primers, **kwargs)
    filtered_rate_df = rate_df[(rate_df['fg_bool']) & (rate_df['bg_bool'])]
    filtered_rate_df = filtered_rate_df.drop(columns=['fg_bool', 'bg_bool'])
    print("Filtered " + str(len(rate_df) - len(filtered_rate_df)) + " number of primers based on foreground/background rate.")

    if src.parameter.verbose:
        print(rate_df)

    print("Computing Gini index...")
    gini_df = src.filter.get_gini(fg_prefixes, fg_genomes, fg_seq_lengths, df = filtered_rate_df)
    print("Filtered " + str(len(filtered_rate_df) - len(gini_df)) + " number of primers based on gini rate.")
    gini_df['ratio'] = gini_df['bg_count']/gini_df['fg_count']
    filtered_gini_df = gini_df.sort_values(by=['ratio'], ascending=False)[:src.parameter.max_primer]

    pickle.dump(filtered_gini_df, open(os.path.join(src.parameter.data_dir, 'step2_df.p'), 'wb'))
    src.string_search.get_positions(filtered_gini_df['primer'], fg_prefixes, fg_genomes)
    src.string_search.get_positions(filtered_gini_df['primer'], bg_prefixes, bg_genomes)
    return filtered_gini_df

#RANK BY RANDOM FOREST
def step3(out_file_prefix=None):
    step2_df = pickle.load(open(os.path.join(src.parameter.data_dir, 'step2_df.p'), 'rb'))

    primer_list = step2_df['primer']
    fg_scale = sum(fg_seq_lengths)/6200

    print("Computing random forest features...")
    df_pred = src.rf_preprocessing.create_augmented_df(fg_prefixes, primer_list)
    df_pred = src.rf_preprocessing.scale_delta_Gs(df_pred, on_scale=fg_scale)
    df_pred['molarity'] = 2.5

    print("Computing amplification efficacy...")
    results = src.rf_preprocessing.predict_new_primers(df_pred)
    results.sort_values(by=['on.target.pred'], ascending=[False], inplace=True)

    results['on.target.pred'] = results['on.target.pred']
    step3_df = results[results['on.target.pred'] >= src.parameter.min_amp_pred]

    step2_df = step2_df.set_index('primer')
    step3_df = step3_df.set_index('sequence')

    joined_step3_df = step3_df.join(step2_df[['ratio', 'gini', 'fg_count', 'bg_count']], how='left').sort_values(by='gini')

    if out_file_prefix:
        pickle.dump(joined_step3_df, open(out_file_prefix + "df_pred.p", "wb"))
        joined_step3_df.to_csv(out_file_prefix+'_regression_results.csv')
    else:
        pickle.dump(joined_step3_df, open(os.path.join(src.parameter.data_dir, 'step3_df.p'), "wb"))

    print("Filtered  " + str(step2_df.shape[0] - joined_step3_df.shape[0]) + " number of primers based on efficacy.")

    if src.parameter.verbose:
        print(joined_step3_df)

    return joined_step3_df

# src.optimize AND SEARCH
def step4(primer_list=None, scores = None, target_var='coverage', selection_metric='deterministic', initial_primer_sets=None, iterations=10, max_sets=10, in_file_prefix=None, out_file=None, fg_circular=True, bg_circular=False, banned_primers=[]):
    if primer_list is None:
        step3_df = pickle.load(open(os.path.join(src.parameter.data_dir, 'step3_df.p'), 'rb'))
        primer_list = step3_df.index
        scores = step3_df['on.target.pred']

    if type(primer_list) is str:
        primer_list = pickle.load(open(primer_list, 'rb'))['sequence']

    all_results = []
    all_scores = []
    cache = {}

    for i in range(src.parameter.iterations):
        print("Repeat #: " + str(i+1))

        if initial_primer_sets is not None:
            primer_list.extend(src.utility.flatten(initial_primer_sets))
        else:
            initial_primer_sets = src.optimize.random_initial_start(np.asarray(primer_list).reshape((-1, 1)), scores, max_sets=max_sets)

        print("Initial primers: ")
        print(initial_primer_sets)

        primer_list = sorted(list(set(primer_list)))
        random_primers = src.utility.flatten(initial_primer_sets)
        combined_primer_list = list(set(primer_list+random_primers))

        print("Banned primers: " + str(banned_primers))
        results, scores, cache = src.optimize.bfs(combined_primer_list, fg_prefixes, bg_prefixes, fg_seq_lengths, bg_seq_lengths, initial_primer_sets=initial_primer_sets, iterations=iterations, max_sets=max_sets, target_var=target_var, selection_method=selection_metric, banned_primers=banned_primers, fg_circular=fg_circular, bg_circular=bg_circular, cache=cache, drop_indices=src.parameter.drop_iterations)

        all_results.append(results)
        all_scores.append(scores)

        banned_primers=[src.utility.most_frequent(src.utility.flatten(results)), src.utility.most_frequent(src.utility.flatten(results))]

    if out_file:
        pickle.dump(results, open(out_file, 'wb'))

    all_results = src.utility.flatten(all_results)
    all_scores = src.utility.flatten(all_scores)

    results = []
    scores = []

    for i, result in enumerate(all_results):
        if result not in result:
            results.append(result)
            scores.append(all_scores[i])

    # scores = src.utility.softmax(scores)

    idx = np.argsort(scores)[::-1]
    final_sets = list(np.array(results)[idx.astype(int)])
    final_scores = list(np.array(scores)[idx.astype(int)])

    # final_sets = final_sets[:src.parameter.top_set_count]
    # final_scores = final_scores[:src.parameter.top_set_count]

    min_score = min(final_scores)

    if min_score < 0:
        final_scores = [score + abs(min_score) for score in final_scores]
    elif min_score > 100:
        final_scores = [score - abs(min_score) for score in final_scores]

    print("FINAL TOP:" + str(src.parameter.top_set_count))
    for i, final_set in enumerate(final_sets[:src.parameter.top_set_count]):
        print('[' + ', '.join(map(str, final_set)) + '], ' + str(round(final_scores[i], 2)))

def step5(primer_sets):

    primer_scores = src.optimize.stand_alone_score_primer_sets(primer_sets)

    argsort_indices = np.argsort(primer_scores)[::-1]

    selection_indices = [0]

    avoid_indices = []

    for index in argsort_indices:
        can_add = True
        for selected_index in selection_indices:
            num_intersect = len(src.utility.intersection(primer_sets[selected_index], primer_sets[index]))
            if num_intersect > 4:
                can_add = False
        if can_add and index not in avoid_indices:
            selection_indices.append(index)

    selected_sets = [primer_sets[selected_index] for selected_index in selection_indices]
    selected_scores = [primer_scores[selected_index] for selected_index in selection_indices]

    for i, selected_set in enumerate(selected_sets):
        print(str(selected_set) + ", score=%0.5f" % selected_scores[i])

if __name__ == "__main__":
    # print("pipeline.py")

    # global params


    # step1()
    # step2()
    # step3()
    step4()
    # print(parameter.min_fg_freq)

    # src.optimize.stand_alone_score_primer_sets([['AAAGTATGG', 'ATTACAGCA', 'GATAAGAGA', 'GTGATGAAA', 'TCAACTATA']], fg_prefixes, bg_prefixes, fg_seq_lengths, bg_seq_lengths)

