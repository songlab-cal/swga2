import pandas as pd
import src.parameter
import src.utility
import scipy.sparse
import numpy as np
import pickle
import src.dimer
from functools import partial
import h5py
import scipy.stats
from sortedcontainers import SortedSet
import os

ridge_regression_model = 'evaluation_models/on_mean_entropy_off_mean_kurtosis.p'

def search_initialization(fg_fname_prefixes, bg_fname_prefixes, fg_seq_lengths, bg_seq_lengths, initial_primer_sets=None, max_sets=10, fg_circular=True, bg_circular=False):
    if initial_primer_sets is None or len(initial_primer_sets) == 0:
        top_sets = [[] for x in range(max_sets)]
        top_scores = [-np.inf for x in range(max_sets)]
        top_fg_fname_to_positions = [
            initialize_fname_to_positions_dict(fg_fname_prefixes, fg_seq_lengths) for x in
            range(max_sets)]
        top_bg_fname_to_positions = [
            initialize_fname_to_positions_dict(bg_fname_prefixes, bg_seq_lengths) for x in
            range(max_sets)]
    else:
        print("Starting foreground search initialization...")
        fg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=fg_fname_prefixes,
                                       seq_lengths_fname=fg_seq_lengths, circular=fg_circular)
        top_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, initial_primer_sets, src.parameter.cpus)
        print("Done with foreground initialization.")

        print("Starting background search initialization...")
        bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_fname_prefixes,
                                       seq_lengths_fname=bg_seq_lengths, circular=bg_circular)
        top_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, initial_primer_sets, src.parameter.cpus)
        print("Done with background search initialization.")

        partial_initialized_f = partial(evaluate_wrapper, fg_seq_lengths=fg_seq_lengths, bg_seq_lengths=bg_seq_lengths)
        top_scores = src.utility.create_pool(partial_initialized_f, list(zip(top_fg_fname_to_positions, top_bg_fname_to_positions)), src.parameter.cpus)

        top_sets = initial_primer_sets

    return top_sets, top_scores, top_fg_fname_to_positions, top_bg_fname_to_positions

def random_initial_start(primers, scores, max_sets):
    selection = make_selection(primers, scores, max_sets, 'normalized')
    next_top_sets = [list(primers[i]) for i in selection]
    return next_top_sets

def get_compressed_string(primer_set):
    compressed_primer_string = ",".join(sorted(primer_set))
    return compressed_primer_string

def drop_out_layer(primer_sets, cache, fg_prefixes=None, bg_prefixes=None, fg_seq_lengths=None, bg_seq_lengths=None, fg_circular=True, bg_circular=False, max_sets=10):
    returning_sets = []
    for primer_set in primer_sets:
        tasks = []
        for i, primer in enumerate(primer_set):
            primer_set_sub = list(primer_set)
            del primer_set_sub[i]
            compressed_primer_string = get_compressed_string(primer_set_sub)
            if compressed_primer_string not in cache:
                tasks.append(primer_set_sub)
                cache[compressed_primer_string] = 0

        fg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=fg_prefixes,
                                       seq_lengths_fname=fg_seq_lengths, circular=fg_circular)
        curr_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, tasks, cpus=src.parameter.cpus)

        bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_prefixes,
                                       seq_lengths_fname=bg_seq_lengths, circular=bg_circular)
        curr_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, tasks, cpus=src.parameter.cpus)

        partial_initialized_f = partial(evaluate_wrapper, fg_seq_lengths=fg_seq_lengths,
                                        bg_seq_lengths=bg_seq_lengths)
        curr_scores = src.utility.create_pool(partial_initialized_f, list(zip(curr_fg_fname_to_positions, curr_bg_fname_to_positions)), cpus=src.parameter.cpus)

        for i, primer_set_sub in enumerate(tasks):
            compressed_primer_string = get_compressed_string(primer_set_sub)
            cache[compressed_primer_string] = curr_scores[i]

    top_sets = [[] for i in range(max_sets)]
    top_scores = [float("-inf") for i in range(max_sets)]

    for primer_set in primer_sets:
        for i, primer in enumerate(primer_set):
            primer_set_sub = list(primer_set)
            del primer_set_sub[i]
            compressed_primer_string = get_compressed_string(primer_set_sub)
            score = cache[compressed_primer_string]
            if score > min(top_scores):
                j = top_scores.index(min(top_scores))
                top_sets[j] = primer_set_sub
                top_scores[j] = score

    fg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=fg_prefixes,
                                   seq_lengths_fname=fg_seq_lengths, circular=fg_circular)
    top_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, top_sets, cpus=src.parameter.cpus)

    bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_prefixes,
                                   seq_lengths_fname=bg_seq_lengths, circular=bg_circular)
    top_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, top_sets, cpus=src.parameter.cpus)

    return top_sets, top_scores, top_fg_fname_to_positions, top_bg_fname_to_positions, cache

#Returns cache of compressed primer string to score, i,
def bfs_one_top_set(**kwargs):
    primer_list = kwargs['primer_list']
    banned_primers = kwargs['banned_primers']
    dimer_mat = kwargs['dimer_mat']
    top_set = kwargs['top_set']
    primer_to_index_dict = kwargs['primer_to_index_dict']
    compressed_string_to_score = kwargs['compressed_string_to_score']
    top_fg_fname_to_positions_sub = kwargs['top_fg_fname_to_positions_sub']
    top_bg_fname_to_positions_sub = kwargs['top_bg_fname_to_positions_sub']
    fg_fname_prefixes = kwargs['fg_fname_prefixes']
    bg_fname_prefixes = kwargs['bg_fname_prefixes']
    fg_seq_lengths = kwargs['fg_seq_lengths']
    bg_seq_lengths = kwargs['bg_seq_lengths']
    max_sets = kwargs['max_sets']
    normalize_metric = kwargs['normalize_metric']


    tasks_sub = []

    # Make all the sets combining the i-th top set and a primer.
    for primer in primer_list:
        # if primer not in banned_primers and dimer.has_no_dimer_risk(top_set + [primer], max_dimer_bp=parameter.default_max_dimer_bp):
        if primer not in banned_primers and src.dimer.compatible(dimer_mat, top_set, primer, primer_to_index_dict):
            # compressed_primer_string = get_compressed_string(top_set + [primer])
            # if compressed_primer_string not in compressed_string_to_score and primer not in top_set:
            if primer not in top_set:
                tasks_sub.append(primer)

    # print("TASKS: " + str(tasks_sub))

    max_sets_sub = 2

    # Run all the tasks for the i-th top set
    if len(tasks_sub) > 0:
        next_top_sets_sub, next_top_scores_sub, next_top_fg_fname_to_positions_sub, next_top_bg_fname_to_positions_sub, cache = parallel_set_search_one_iteration(tasks_sub, top_set,
                                                                                   top_fg_fname_to_positions_sub,
                                                                                   top_bg_fname_to_positions_sub,
                                                                                   fg_fname_prefixes, bg_fname_prefixes,
                                                                                   fg_seq_lengths, bg_seq_lengths,
                                                                                   max_sets_sub, normalize_metric=normalize_metric, cache=compressed_string_to_score)
    else:
        next_top_sets_sub = [[] for x in range(max_sets)]
        next_top_scores_sub = [0 for x in range(max_sets)]
        next_top_fg_fname_to_positions_sub = [
            initialize_fname_to_positions_dict(fg_fname_prefixes, fg_seq_lengths) for x in
            range(max_sets)]
        next_top_bg_fname_to_positions_sub = [
            initialize_fname_to_positions_dict(bg_fname_prefixes, bg_seq_lengths) for x in
            range(max_sets)]
    return next_top_sets_sub, next_top_scores_sub, next_top_fg_fname_to_positions_sub, next_top_bg_fname_to_positions_sub, compressed_string_to_score

# Searches for primer sets in a greedy fashion for a specified number of iterations, keeping track of ten sets at a time. The method of selection for each iteration can also be specified as "normalized" (probability of selection is the score divided by the sum of all scores), "softmax", or "deterministic" (no randomization).
def bfs(primer_list, fg_fname_prefixes, bg_fname_prefixes, fg_seq_lengths, bg_seq_lengths, scores=None, initial_primer_sets=None, iterations=10, max_sets=10, target_var='coverage', selection_method='deterministic', banned_primers=[], fg_circular=True, bg_circular=False, cache={}, drop_indices=[4]):
    print("Method: " + target_var + ', ' + selection_method)
    primer_list = sorted(list(primer_list))
    primer_to_index_dict = dict(zip(primer_list, range(len(primer_list))))
    dimer_mat = src.dimer.heterodimer_matrix(primer_list, max_dimer_bp=src.parameter.max_dimer_bp)

    top_sets, top_scores, top_fg_fname_to_positions, top_bg_fname_to_positions = search_initialization(fg_fname_prefixes, bg_fname_prefixes,
                                                                                 fg_seq_lengths, bg_seq_lengths,
                                                                                 initial_primer_sets=initial_primer_sets,
                                                                                 max_sets=max_sets,
                                                                                 fg_circular=fg_circular,
                                                                                 bg_circular=bg_circular)

    finished_sets = []  # Keeps track of the best sets out of all iterations.
    finished_scores = []

    if src.parameter.verbose:
        print("The banned primers: " + str(banned_primers))

    # choose the one with the maximum score
    for iter in range(1, iterations):

        next_top_sets_all = []  #Selected sets of the iteration using the ith top set.
        next_top_scores_all = []
        next_top_fg_fname_to_positions_all = []
        next_top_bg_fname_to_positions_all = []
        prev_top_set_index_all = []

        print("Iteration #: " + str(iter))

        if src.parameter.verbose:
            print(top_scores)

        for i, top_set in enumerate(top_sets):

            # Get the positions for the ith top set
            top_fg_fname_to_positions_sub = top_fg_fname_to_positions[i]
            top_bg_fname_to_positions_sub = top_bg_fname_to_positions[i]

            next_top_sets_sub, next_top_scores_sub, next_top_fg_fname_to_positions_sub, next_top_bg_fname_to_positions_sub, cache = \
                bfs_one_top_set(primer_list=primer_list,
                                banned_primers=banned_primers,
                                dimer_mat=dimer_mat,
                                top_set=top_set,
                                primer_to_index_dict=primer_to_index_dict,
                                top_fg_fname_to_positions_sub=top_fg_fname_to_positions_sub,
                                top_bg_fname_to_positions_sub=top_bg_fname_to_positions_sub,
                                fg_fname_prefixes=fg_fname_prefixes,
                                bg_fname_prefixes=bg_fname_prefixes,
                                fg_seq_lengths=fg_seq_lengths,
                                bg_seq_lengths=bg_seq_lengths,
                                max_sets=max_sets,
                                normalize_metric='deterministic',
                                compressed_string_to_score=cache)

            next_top_sets_all.extend(next_top_sets_sub)
            next_top_scores_all.extend(next_top_scores_sub)
            next_top_fg_fname_to_positions_all.extend(next_top_fg_fname_to_positions_sub)
            next_top_bg_fname_to_positions_all.extend(next_top_bg_fname_to_positions_sub)
            prev_top_set_index_all.extend([i for x in range(len(next_top_sets_sub))])

        if max(next_top_scores_all) < min(top_scores) + 0.01:
            finished_sets.extend(top_sets)
            finished_scores.extend(top_scores)

            print("Finished sets:")
            print(', '.join(map(str, finished_sets)))
            print(finished_scores)
            print()

            return finished_sets, finished_scores, cache

        selection = make_selection(next_top_sets_all, next_top_scores_all, max_sets, selection_method)

        print("Selections for iteration " + str(iter))
        for j in selection:
            # print(prev_top_set_index_all[j])
            old_set = top_sets[prev_top_set_index_all[j]]
            new_set = next_top_sets_all[j]
            new_primer = set(new_set) - set(old_set)
            print("From top set number " + str(prev_top_set_index_all[j]) + ", added primer " + ', '.join(map(str, new_primer)) + " to [" + ', '.join(map(str, old_set)) + ']')

        threshold = max(next_top_scores_all)
        for i, top_set in enumerate(top_sets):
            #CHANGE THIS BEFORE SHIPPING
            if top_scores[i] > threshold or len(top_set) == 7:
                finished_scores.append(top_scores[i])
                finished_sets.append(top_set)

        if iter in drop_indices:
            top_sets, top_scores, top_fg_fname_to_positions, top_bg_fname_to_positions, cache = drop_out_layer([next_top_sets_all[i] for i in selection], cache, fg_prefixes=fg_fname_prefixes, bg_prefixes=bg_fname_prefixes, fg_seq_lengths=fg_seq_lengths, bg_seq_lengths=bg_seq_lengths, fg_circular=fg_circular, bg_circular=bg_circular, max_sets=max_sets)
            banned_primers = []

        else:
            top_sets = [next_top_sets_all[i] for i in selection]
            top_scores = [next_top_scores_all[i] for i in selection]
            top_fg_fname_to_positions = [next_top_fg_fname_to_positions_all[i] for i in selection]
            top_bg_fname_to_positions = [next_top_bg_fname_to_positions_all[i] for i in selection]

    finished_sets.extend(top_sets)
    finished_scores.extend(top_scores)

    # if sum(finished_scores) < 0:
    # finished_scores = [src.utility.sigmoid(x) for x in finished_scores]
    # print(finished_scores)

    print("Finished sets:")
    for finished_set in finished_sets:
        print('[' + ', '.join(map(str, finished_set))+']')

    return finished_sets, finished_scores, cache

#For each genome file, create a tuple of two sorted sets--one for forward strand and the other for the reverse strand.
def initialize_fname_to_positions_dict(fname_prefixes, seq_lengths):
    curr_fname_to_positions = []

    for prefix, seq_length in zip(fname_prefixes, seq_lengths):
        curr_fname_to_positions.append((SortedSet([]),SortedSet([])))

    return curr_fname_to_positions

#Parallelizes evaluating adding all valid primers to a top set from the previous iteration.
def parallel_set_search_one_iteration(tasks, top_set, top_fg_fname_to_positions, top_bg_fname_to_positions, fg_fname_prefixes, bg_fname_prefixes, fg_seq_lengths, bg_seq_lengths, max_sets, normalize_metric='softmax', cache={}):
    seq_initialized_f = partial(evaluate_pairs_helper, curr_fg_fname_to_positions=top_fg_fname_to_positions, curr_bg_fname_to_positions = top_bg_fname_to_positions, fg_fname_prefixes=fg_fname_prefixes, bg_fname_prefixes=bg_fname_prefixes, fg_seq_lengths=fg_seq_lengths, bg_seq_lengths=bg_seq_lengths)
    results = src.utility.create_pool(seq_initialized_f, tasks, src.parameter.cpus)

    selection=make_selection([[task] + top_set for task in tasks], [score for score, temp_fg, temp_bg in results], min(max_sets, len(tasks)), normalize_metric)

    for i, result in enumerate(results):
        primer_set = sorted([tasks[i]] + top_set)
        compressed_string = get_compressed_string(primer_set)
        cache[compressed_string] = results[i][0]

    return [sorted([tasks[i]] + top_set) for i in selection], [results[i][0] for i in selection], [results[i][1] for i in selection], [results[i][2] for i in selection], cache

def make_selection(all_sets, all_scores, max_sets, normalize_metric):
    sets = []
    scores = []
    all_compressed_string = set()

    for i, primer_set in enumerate(all_sets):
        compressed_string = get_compressed_string(primer_set)
        if compressed_string not in all_compressed_string:
            sets.append(primer_set)
            scores.append(all_scores[i])
            all_compressed_string.add(compressed_string)

    num_tasks = len(scores)
    if normalize_metric == 'deterministic':
        argsort_indices = np.argsort(scores)
        selection = argsort_indices[-max_sets:]
    elif normalize_metric == 'softmax':
        probabilities = src.utility.softmax(scores)
        nonzero_count = np.count_nonzero(probabilities)
        if nonzero_count != 0:
            selection = np.random.choice(range(num_tasks), size=min(nonzero_count,max_sets), replace=False, p=probabilities)
        else:
            selection=[]
    elif normalize_metric == 'normalized':
        probabilities =  scores/np.sum(scores)
        nonzero_count = np.count_nonzero(probabilities)
        if nonzero_count != 0:
            selection = np.random.choice(range(num_tasks), size=min(nonzero_count,max_sets), replace=False, p=probabilities)
        else:
            selection=[]
    selection = list(selection)
    if len(selection) < max_sets and num_tasks> max_sets - len(selection):
        extra_samples = max_sets - len(selection)
        argsort_indices = np.argsort([result for result in scores])
        diff = list(set(argsort_indices).difference(set(selection)))
        if extra_samples <= len(diff):
            selection.extend(diff[-extra_samples:])
        elif len(diff) > 0:
            selection.extend(diff)

    return selection

def get_features_simplified(fname_to_positions, seq_lengths, target=True, circular=False):

    matches = 0
    gap_sizes = []
    gap_ginis = []

    within_mean_gaps = []
    agnostic_mean_gaps = []

    for i, positions in enumerate(fname_to_positions):
        positions_forward = positions[0]
        positions_reverse = positions[1]

        matches += len(positions_forward)
        matches += len(positions_reverse)

        data_forward = get_position_features(positions_forward, seq_lengths[i])
        data_reverse = get_position_features(positions_reverse, seq_lengths[i])

        gap_sizes.append((data_forward['mean'] + data_reverse['mean'])/seq_lengths[i]/2)
        gap_ginis.append((data_forward['gini'] + data_reverse['gini'])/2)

        alternating_gaps_data = get_positional_gap_lengths_alternating(positions_forward, positions_reverse, seq_length=seq_lengths[i])
        within_mean_gaps.append(alternating_gaps_data['within_mean']/seq_lengths[i])

        agnostic_gaps_data = get_positional_gap_lengths_agnostic(positions_forward, positions_reverse, seq_length=seq_lengths[i])
        agnostic_mean_gaps.append(agnostic_gaps_data['mean']/seq_lengths[i])

    if target:
        prefix = 'on_'
    else:
        prefix = 'off_'

    total_seq_length = sum(seq_lengths)

    data = {}

    data['sum'] = matches/total_seq_length

    if len(gap_sizes) > 0:
        data['mean_gap'] = np.mean(gap_sizes)
        data['gap_gini'] = np.mean(gap_ginis)
    else:
        data['mean_gap'] = np.sum(seq_lengths)
        data['gap_gini'] = 0

    if len(agnostic_mean_gaps) > 0:
        data['agnostic_mean_gap'] = np.mean(agnostic_mean_gaps)
    else:
        data['agnostic_mean_gap'] = 0

    if len(within_mean_gaps) > 0:
        data['within_mean_gap'] = np.mean(within_mean_gaps)
    else:
        data['within_mean_gap'] = 0

    position_features = pd.DataFrame([data.values()],columns=[prefix + column_name for column_name in data.keys()])

    return position_features

def get_features(fname_to_positions, seq_lengths, target=True):

    matches = 0
    gap_sizes = []
    gap_ginis = []

    within_mean_gaps = []
    agnostic_mean_gaps = []
    agnostic_mean_ginis = []
    agnostic_mean_entropies = []

    for i, positions in enumerate(fname_to_positions):
        positions_forward = positions[0]
        positions_reverse = positions[1]

        matches += len(positions_forward)
        matches += len(positions_reverse)

        data_forward = get_position_features(positions_forward, seq_lengths[i])
        data_reverse = get_position_features(positions_reverse, seq_lengths[i])

        gap_sizes.append((data_forward['mean'] + data_reverse['mean'])/seq_lengths[i]/2)
        gap_ginis.append((data_forward['gini'] + data_reverse['gini'])/2)

        alternating_gaps_data = get_positional_gap_lengths_alternating(positions_forward, positions_reverse, seq_length=seq_lengths[i])
        within_mean_gaps.append(alternating_gaps_data['within_mean']/seq_lengths[i])

        agnostic_gaps_data = get_positional_gap_lengths_agnostic(positions_forward, positions_reverse, seq_length=seq_lengths[i], include_gini=True)

        agnostic_mean_gaps.append(agnostic_gaps_data['mean']/seq_lengths[i])
        agnostic_mean_ginis.append(agnostic_gaps_data['gini'])
        agnostic_mean_entropies.append(agnostic_gaps_data['entropy']/seq_lengths[i])

    if target:
        prefix = 'on_'
    else:
        prefix = 'off_'

    total_seq_length = sum(seq_lengths)

    data = {}

    data['sum'] = matches/total_seq_length

    if len(gap_sizes) > 0:
        data['mean_gap'] = np.mean(gap_sizes)
        data['gap_gini'] = np.mean(gap_ginis)
    else:
        data['mean_gap'] = np.sum(seq_lengths)
        data['gap_gini'] = 0

    if len(agnostic_mean_gaps) > 0:
        data['agnostic_mean_gap'] = np.mean(agnostic_mean_gaps)
    else:
        data['agnostic_mean_gap'] = 0

    if len(agnostic_mean_ginis) > 0:
        data['agnostic_mean_gini'] = np.mean(agnostic_mean_ginis)
    else:
        data['agnostic_mean_gini'] = 1

    if len(agnostic_mean_entropies) > 0:
        data['agnostic_mean_entropy'] = np.mean(agnostic_mean_entropies)
    else:
        data['agnostic_mean_entropy'] = 1

    if len(within_mean_gaps) > 0:
        data['within_mean_gap'] = np.mean(within_mean_gaps)
    else:
        data['within_mean_gap'] = 0

    position_features = pd.DataFrame([data.values()],columns=[prefix + column_name for column_name in data.keys()])

    return position_features

def get_positional_gap_lengths_alternating(positions_forward, positions_reverse, seq_length=None, threshold=100000, extra=False):

    within_differences = []

    i = 0
    j = 0

    while i < len(positions_forward) and j < len(positions_reverse):
        while j < len(positions_reverse) and positions_reverse[j] < positions_forward[i]:
            j += 1

        sub_j = j

        while sub_j < len(positions_reverse) and positions_reverse[sub_j] - positions_forward[i] <= threshold and positions_reverse[sub_j] >= positions_forward[i]:
            within_differences.append(positions_reverse[sub_j] - positions_forward[i])
            sub_j += 1

        i += 1

    data = {}

    if len(within_differences) == 0:
        data['within_mean'] = 0
    else:
        data['within_mean'] = np.sum(within_differences)  # mean/seq_length

    return data

def get_positional_gap_lengths_agnostic(positions_forward, positions_reverse, seq_length=None, include_gini=False):
    positions = SortedSet(positions_forward)
    positions.update(positions_reverse)

    position_differences = get_positional_gap_lengths(positions, seq_length = seq_length)

    data = {}

    if len(position_differences) <= 1:
        return {'mean': seq_length, 'entropy': 1, 'gini': 1}
    else:
        data['mean'] = np.mean(position_differences)
        if include_gini:
            data['gini'] = src.utility.gini_exact(positions)
            data['entropy'] = scipy.stats.entropy(positions)
    return data


def get_positional_gap_lengths(positions, seq_length=None):

    differences = [a_i - b_i for a_i, b_i in zip(positions[1:], positions[:-1])]

    if differences == []:
        return positions

    if seq_length is None:
        print("Where is the seq length in get_positional_gap_lengths?")

    differences.append(seq_length-positions[-1])
    differences.append(positions[0])

    return differences

def get_position_features(positions, seq_length):

    position_differences = get_positional_gap_lengths(positions, seq_length=seq_length)
    data = {}

    if len(position_differences) <= 1:
        return {'mean': seq_length, 'gini': 1}
    else:
        data['mean'] = np.mean(position_differences)
        data['gini'] = src.utility.gini_exact(position_differences)

    return data

def get_exact_count_features_one_fname(primer_set, fname_prefix):
    curr_pos_results_forward = []
    curr_pos_results_reverse = []

    for primer in primer_set:
        k = len(primer)
        if os.path.exists(fname_prefix + '_' + str(k) + 'mer_positions.h5'):
            db = h5py.File(fname_prefix + '_' + str(k) + 'mer_positions.h5', 'r')
            if primer in db:
                positions = db[primer]
                curr_pos_results_forward.extend(positions)
            if src.utility.reverse_complement(primer) in db:
                positions = db[src.utility.reverse_complement(primer)]
                curr_pos_results_reverse.extend(positions)
        else:
            print(fname_prefix + '_' + str(k) + 'mer_positions.h5 does not exist.')
    return curr_pos_results_forward, curr_pos_results_reverse

def get_fname_to_positions_dict(fname_prefixes, primer_list):

    fname_to_positions = {}

    for fname_prefix in fname_prefixes:
        positions_forward, positions_reverse = get_exact_count_features_one_fname(primer_list, fname_prefix)
        fname_to_positions[fname_prefix] = (positions_forward, positions_reverse)
        # matches += len(positions)
    return fname_to_positions

def incorporate_positions_row(curr_fname_to_positions, fname_prefixes, primer):
    # print(curr_fname_to_positions)
    fname_to_positions = []
    total_matches = 0

    for curr_positions, fname_prefix in zip(curr_fname_to_positions, fname_prefixes):

        curr_positions_forward = curr_positions[0].copy()
        curr_positions_reverse = curr_positions[1].copy()

        positions_forward, positions_reverse = get_exact_count_features_one_fname([primer], fname_prefix)

        total_matches += len(positions_forward)
        total_matches += len(positions_reverse)

        curr_positions_forward.update(positions_forward)
        curr_positions_reverse.update(positions_reverse)

        fname_to_positions.append((curr_positions_forward, curr_positions_reverse))
    return fname_to_positions, total_matches

def evaluate(temp_fg_fname_to_positions, temp_bg_fname_to_positions, fg_seq_lengths=None, bg_seq_lengths=None, target_var='percent', return_X=False, fg_circular=True, bg_circular=False):
    fg_features = get_features_simplified(temp_fg_fname_to_positions, fg_seq_lengths, circular=fg_circular, target=True)

    bg_features = get_features_simplified(temp_bg_fname_to_positions, bg_seq_lengths, circular=bg_circular, target=False)

    pd.set_option('display.max_columns', 500)

    if fg_features.empty or bg_features.empty:
        return [-np.inf], None

    all_features = pd.concat([fg_features, bg_features], axis=1)
    all_features = all_features.reset_index(drop=True)

    all_features['off_sum'] = all_features['off_sum'].replace(0, 0.0000000001)
    all_features['on_within_mean_gap'] = all_features['on_within_mean_gap'].replace(0, 0.0000000001)

    all_features['ratio'] = all_features['on_sum']/all_features['off_sum']
    all_features['within_mean_gap_ratio'] = all_features['off_within_mean_gap'] / all_features['on_within_mean_gap']
    all_features['agnostic_mean_gap_ratio'] = all_features['on_agnostic_mean_gap'] / all_features['off_agnostic_mean_gap']

    # print(all_features)

    X = all_features[['ratio', 'agnostic_mean_gap_ratio', 'on_gap_gini', 'off_gap_gini', 'within_mean_gap_ratio']]
    clf = pickle.load(open(os.path.join(src.parameter.src_dir, 'ratio_agnostic_mean_gap_ratio_all_ginis_within_mean_gap_ratio.p'),
        'rb'))

    score = clf.predict(X)[0]

    if return_X:
        return score, X

    return score

def evaluate_pairs_helper(tasks, curr_fg_fname_to_positions=None, curr_bg_fname_to_positions=None, fg_fname_prefixes=None, bg_fname_prefixes=None, fg_seq_lengths=None, bg_seq_lengths=None, target_var=None):
    new_primer = tasks

    temp_fg_fname_to_positions, _ = incorporate_positions_row(curr_fg_fname_to_positions, fg_fname_prefixes, new_primer)
    temp_bg_fname_to_positions, _ = incorporate_positions_row(curr_bg_fname_to_positions, bg_fname_prefixes, new_primer)

    score = evaluate(temp_fg_fname_to_positions, temp_bg_fname_to_positions, fg_seq_lengths=fg_seq_lengths, bg_seq_lengths=bg_seq_lengths, target_var=target_var)
    return score, temp_fg_fname_to_positions, temp_bg_fname_to_positions

def build_fname_to_positions(primer_set, fname_prefixes=None, seq_lengths_fname=None, circular=True):
    temp_fname_to_positions = initialize_fname_to_positions_dict(fname_prefixes, seq_lengths_fname)
    for primer in primer_set:
        temp_fname_to_positions, total_matches = incorporate_positions_row(temp_fname_to_positions, fname_prefixes, primer)
    return temp_fname_to_positions

def evaluate_wrapper(fname_to_positions, fg_seq_lengths=None, bg_seq_lengths=None, target_var='coverage'):
    return evaluate(fname_to_positions[0], fname_to_positions[1], fg_seq_lengths=fg_seq_lengths,bg_seq_lengths=bg_seq_lengths, target_var=target_var)

def stand_alone_score_primer_sets(primer_sets, fg_prefixes, bg_prefixes, fg_seq_lengths, bg_seq_lengths, fg_circular=True, bg_circular=False):
    human_chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                      'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM',
                      'chrX', 'chrY']

    fg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=fg_prefixes,
                                   seq_lengths_fname=fg_seq_lengths, circular=fg_circular)
    initial_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, primer_sets, cpus=src.parameter.cpus)

    bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_prefixes,
                                   seq_lengths_fname=bg_seq_lengths, circular=bg_circular)
    initial_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, primer_sets, cpus=src.parameter.cpus)

    partial_initialized_f = partial(evaluate_wrapper, fg_seq_lengths=fg_seq_lengths,
                                    bg_seq_lengths=bg_seq_lengths)
    top_scores = src.utility.create_pool(partial_initialized_f,
                                     list(zip(initial_fg_fname_to_positions, initial_bg_fname_to_positions)), cpus=src.parameter.cpus)
    print(primer_sets)
    print(top_scores)
    return top_scores

