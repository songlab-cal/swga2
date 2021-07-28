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
import warnings

# Path to the saved model for evaluating the primer set.
ridge_regression_model = 'evaluation_models/on_mean_entropy_off_mean_kurtosis.p'

def search_initialization(fg_fname_prefixes, bg_fname_prefixes, fg_seq_lengths, bg_seq_lengths, initial_primer_sets=None, max_sets=10, fg_circular=True, bg_circular=False):
    """
    Initializes the optimization process by doing the following steps:
    1. Sets the top_sets to a list of empty sets, where the number of top sets to keep in every iteration is determined by max_sets.
    If the user wants to input a set of initial primer_sets to build off of, the top_sets will be initialized to this.
    2. Computes the scores for the top_sets and intializes these as top_scores.
    3. Collects the positions that have an exact match in all the genomes to the primers in the top_sets.

    Args:
        fg_fname_prefixes: The path prefixes to all the h5py files relevant to the on-target genome.
        bg_fname_prefixes: The path prefixes to all the h5py files relevant to the off-target genome.
        fg_seq_lengths: The lengths of the on-target genomes in the fasta files (in the same order as fg_fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes in the fasta files (in the same order as bg_fname_prefixes).
        initial_primer_sets: The primer sets that we want to build off of. If set to None, the program will start from scratch.
        max_sets: The number of top sets to keep in every iteration. Set to 10 by default.
        fg_circular: Boolean variable indicating whether the on-target genome is circular. Set to true by default.
        bg_circular: Boolean variable indicating whether the off-target genome is circular. Set to false by default.

    Returns:
        top_sets: Either a list of empty lists or the initial primer sets input by the user.
        top_scores: A list of scores for the top_sets.
        top_fg_fname_to_positions: List of lists of tuples of sorted sets containing positions of exact binding
        locations in the forward and reverse strand, respectively, where the first index corresponds to the top set
        and the second index corresponds to the h5py path prefix.
        top_bg_fname_to_positions: List of lists of tuples of sorted sets containing positions of exact binding
        locations in the forward and reverse strand, respectively, where the first index corresponds to the top set
        and the second index corresponds to the h5py path prefix.
    """
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
                                       seq_lengths=fg_seq_lengths, circular=fg_circular)
        top_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, initial_primer_sets, src.parameter.cpus)
        print("Done with foreground initialization.")

        print("Starting background search initialization...")
        bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_fname_prefixes,
                                       seq_lengths=bg_seq_lengths, circular=bg_circular)
        top_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, initial_primer_sets, src.parameter.cpus)
        print("Done with background search initialization.")

        partial_initialized_f = partial(evaluate_wrapper, fg_seq_lengths=fg_seq_lengths, bg_seq_lengths=bg_seq_lengths, fg_circular=fg_circular, bg_circular=bg_circular)
        top_scores = src.utility.create_pool(partial_initialized_f, list(zip(top_fg_fname_to_positions, top_bg_fname_to_positions)), src.parameter.cpus)
        top_sets = initial_primer_sets

    return top_sets, top_scores, top_fg_fname_to_positions, top_bg_fname_to_positions

def random_initial_start(primer_sets, scores, max_sets):
    """
    This function picks the top set according to a probability based on the scores.

    Args:
        primer_sets: The primer_sets from which to choose from.
        scores: The corresponding scores of the primer sets.
        max_sets: The number of sets to choose.

    Returns:
        next_top_sets: A list of length max_sets constituting the selected primer sets.
    """
    selection = make_selection(primer_sets, scores, max_sets, 'normalized')
    next_top_sets = [list(primer_sets[i]) for i in selection]
    return next_top_sets

def get_compressed_string(primer_set):
    """
     This function returns the primers joined together as a single string delimited by a comma. This is mostly to check
     for redundant sets.

     Args:
         primer_set: The primer set to be compressed

     Returns:
         compressed_primer_string: Single string joining the primers in the set by commas.
     """
    compressed_primer_string = ",".join(sorted(primer_set))
    return compressed_primer_string

def drop_out_layer(primer_sets, cache, fg_prefixes=None, bg_prefixes=None, fg_seq_lengths=None, bg_seq_lengths=None, fg_circular=True, bg_circular=False, max_sets=10):
    """
    This function computes the scores of all subsets of primer_sets one size smaller than than those in primer_sets.
    It returns the subsets with the best the highest scores. The "drop out" layer is primarily for allowing backtracking
    in case the greedy optimization prematurely added a non-optimal primer and got stuck in a local optima. 

    Args:
        primer_sets: The list of primer sets to compute subsets from.
        cache: This is a dictionary mapping the compressed string representation of the primer set to its score.
        This prevents recomputing the score for the same primer set everytime.
        fg_prefixes: The path prefixes to all the h5py files relevant to the on-target genome.
        bg_prefixes: The path prefixes to all the h5py files relevant to the off-target genome.
        fg_seq_lengths: The lengths of the on-target genomes in the fasta files (in the same order as fg_fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes in the fasta files (in the same order as bg_fname_prefixes).
        fg_circular: Boolean variable indicating whether the on-target genome is circular. Set to true by default.
        bg_circular: Boolean variable indicating whether the off-target genome is circular. Set to false by default.
        max_sets: The number of top sets to keep in every iteration. Set to 10 by default.

    Returns:
        top_sets: The top subsets selected.
        top_scores: The scores of the top subsets selected.
        top_fg_fname_to_positions: List of lists of tuples of sorted sets containing positions of exact binding
        locations in the forward and reverse strand, respectively, where the first index corresponds to the top set
        and the second index corresponds to the h5py path prefix.
        top_bg_fname_to_positions: List of lists of tuples of sorted sets containing positions of exact binding
        locations in the forward and reverse strand, respectively, where the first index corresponds to the top set
        and the second index corresponds to the h5py path prefix.
        cache: Dictionary mapping the compressed string representation of the primer set to its score.
        This prevents recomputing the score for the same primer set everytime. May or may not have been modified
        during the function call depending on if a primer set was missing from the cache or not.
    """
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
                                       seq_lengths=fg_seq_lengths, circular=fg_circular)
        curr_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, tasks, cpus=src.parameter.cpus)

        bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_prefixes,
                                       seq_lengths=bg_seq_lengths, circular=bg_circular)
        curr_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, tasks, cpus=src.parameter.cpus)

        partial_initialized_f = partial(evaluate_wrapper, fg_seq_lengths=fg_seq_lengths,
                                        bg_seq_lengths=bg_seq_lengths, fg_circular=fg_circular, bg_circular=bg_circular)
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
                                   seq_lengths=fg_seq_lengths, circular=fg_circular)
    top_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, top_sets, cpus=src.parameter.cpus)

    bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_prefixes,
                                   seq_lengths=bg_seq_lengths, circular=bg_circular)
    top_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, top_sets, cpus=src.parameter.cpus)

    return top_sets, top_scores, top_fg_fname_to_positions, top_bg_fname_to_positions, cache

def bfs_one_top_set(**kwargs):
    """
    This function runs the optimization search by building off of one primer set. It creates the new primer sets from
    adding individual candidate primers to the input primer set. We compute the scores (checking that we haven't
    computed it already and return the top two primer sets built from adding to the primer set.

    Args:
        kwargs: Dictionary of arguments containing the following keys and values.
            primer_list: The primer set from which to build new primer sets from.
            banned_primers: A list of primers we don't want to bother adding to the primer set.
            dimer_mat: A 2D binary, symmetric array where a 1 in row i and column j indicates primer i and primer j
            (from primer_list)may form a heterodimer.
            top_set: The primer set off of which to build new primer sets from.
            primer_to_index_dict: A map from primer sequence to primer ID.
            compressed_string_to_score: Dictionary mapping the compressed string representation of the primer set to its score.
            This prevents recomputing the score for the same primer set everytime.
            top_fg_fname_to_positions_sub: List of a tuple of sorted sets containing positions of exact binding
            locations in the forward and reverse strand of the on-target genome, respectively, where the index
            corresponds to the h5py path prefix.
            top_bg_fname_to_positions_sub: List of a tuple of sorted sets containing positions of exact binding
            locations in the forward and reverse strand of the off-target genome, respectively, where the index
            corresponds to the h5py path prefix.
            fg_fname_prefixes: The path prefixes to all the h5py files relevant to the on-target genome.
            bg_fname_prefixes: The path prefixes to all the h5py files relevant to the off-target genome.
            fg_seq_lengths: The lengths of the on-target genomes in the fasta files (in the same order as fg_fname_prefixes).
            bg_seq_lengths: The lengths of the off-target genomes in the fasta files (in the same order as bg_fname_prefixes).
            max_sets: The number of top sets to keep in every iteration. Set to 10 by default.
            normalize_metric:
    Returns:
        next_top_sets_sub: The top sets selected.
        next_top_scores_sub: The scores of the top sets selected.
        next_top_fg_fname_to_positions_sub: List of lists of a tuple of sorted sets containing positions of exact binding
        locations in the forward and reverse strand, respectively, of the on-target where the first index corresponds to the top set
        and the second index corresponds to the h5py path prefix.
        next_top_bg_fname_to_positions_sub: List of lists of tuples of sorted sets containing positions of exact binding
        locations in the forward and reverse strand, respectively, of the off-target where the first index corresponds to the top set
        and the second index corresponds to the h5py path prefix.
        compressed_string_to_score: Dictionary mapping the compressed string representation of the primer set to its score.
        This prevents recomputing the score for the same primer set everytime. May or may not have been modified
        during the function call depending on if a primer set was missing from the cache or not.
    """
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
    fg_circular = kwargs['fg_circular']
    bg_circular = kwargs['bg_circular']

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
                                                                                   max_sets_sub, normalize_metric=normalize_metric, cache=compressed_string_to_score, fg_circular=fg_circular, bg_circular=bg_circular)
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
def bfs(primer_list, fg_fname_prefixes, bg_fname_prefixes, fg_seq_lengths, bg_seq_lengths, initial_primer_sets=None, iterations=10, max_sets=10, selection_method='deterministic', banned_primers=[], fg_circular=True, bg_circular=False, cache={}, drop_indices=[4]):
    """
    This is the master function for launching optimization of primer sets given a list of candidate primers.

    Args:
        primer_list:
        fg_fname_prefixes: The path prefixes to all the h5py files relevant to the on-target genome.
        bg_fname_prefixes: The path prefixes to all the h5py files relevant to the off-target genome.
        fg_seq_lengths: The lengths of the on-target genomes in the fasta files (in the same order as fg_fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes in the fasta files (in the same order as bg_fname_prefixes).
        initial_primer_sets: The primer sets that we want to build off of. If set to None, the program will start from scratch.
        iterations: The number of iterations to try to run the optimization. It will prematurely halt if the scores of the primer sets are no longer increasing. By default set to 10.
        max_sets: The number of top sets to keep in every iteration. Set to 10 by default.
        selection_method: The method by which to pick the next iterations set of best primers. See the function make_selection.
        banned_primers: A list of primers which we do not want to consider. This was added because primers in the 
        initial_primer_sets were also undesiredly being considered for other future primer sets. 
        fg_circular: Boolean variable indicating whether the on-target genome is circular. Set to true by default.
        bg_circular: Boolean variable indicating whether the off-target genome is circular. Set to false by default.
        cache: Dictionary mapping the compressed string representation of the primer set to its score.
        This prevents recomputing the score for the same primer set everytime. 
        drop_indices: A list of iterations which will have a drop out layer at the end. See drop_out_layer. 
        If the indices are too early, it may be meaningless, but computational time increases as the index grows.

    Returns:
        finished_sets: The primer sets finally chosen at the end. 
        finished_scores: The scores of the chosen primer sets. 
        cache: Dictionary mapping the compressed string representation of the primer set to its score.
        This prevents recomputing the score for the same primer set everytime. May or may not have been modified
        during the function call depending on if a primer set was missing from the cache or not, but it may be useful 
        to see the optimization paths.
    """
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
                                compressed_string_to_score=cache, 
                                fg_circular=fg_circular, bg_circular=bg_circular)

            next_top_sets_all.extend(next_top_sets_sub)
            next_top_scores_all.extend(next_top_scores_sub)
            next_top_fg_fname_to_positions_all.extend(next_top_fg_fname_to_positions_sub)
            next_top_bg_fname_to_positions_all.extend(next_top_bg_fname_to_positions_sub)
            prev_top_set_index_all.extend([i for x in range(len(next_top_sets_sub))])

        if max(next_top_scores_all) < min(top_scores) + 0.01:
            finished_sets.extend(top_sets)
            finished_scores.extend(top_scores)

            print("Finished sets for this restart:")
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

            if len(new_primer) == 1:
                print("From top set number " + str(prev_top_set_index_all[j]) + ", added primer " + ', '.join(map(str, new_primer)) + " to [" + ', '.join(map(str, old_set)) + ']')

        threshold = max(next_top_scores_all)
        for i, top_set in enumerate(top_sets):
            if top_scores[i] > threshold:
                finished_scores.append(top_scores[i])
                finished_sets.append(top_set)

        if iter == drop_indices:
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

    print("Finished sets for this restart:")
    for finished_set in finished_sets:
        if len(finished_set) > 0:
            print('[' + ', '.join(map(str, finished_set))+']')
    print(finished_scores)

    return finished_sets, finished_scores, cache

def initialize_fname_to_positions_dict(fname_prefixes, seq_lengths):
    """
    For each path prefix to a set of h5py files corresponding to one genome, create a tuple of two sorted sets--one for
    forward strand and the other for the reverse strand.

    Args:
        fname_prefixes: The h5py path prefixes for the genomes of interest. Each path prefix corresponds to one fasta file.
        seq_lengths: The lengths of the genomes in the corresond fasta files (in the same order as fname_prefixes).

    Returns:
        curr_fname_to_positions: List of tuples of empty two sorted sets, where the index corresponds to the h5py path prefix.
    """
    curr_fname_to_positions = []

    for prefix, seq_length in zip(fname_prefixes, seq_lengths):
        curr_fname_to_positions.append((SortedSet([]),SortedSet([])))

    return curr_fname_to_positions

#Parallelizes evaluating adding all valid primers to a top set from the previous iteration.
def parallel_set_search_one_iteration(tasks, top_set, top_fg_fname_to_positions, top_bg_fname_to_positions, fg_fname_prefixes, bg_fname_prefixes, fg_seq_lengths, bg_seq_lengths, max_sets, normalize_metric='softmax', cache={}, fg_circular=True, bg_circular=False):
    """
    This is a helper function for multiprocessing computation and evaluation of all primer sets built from one of the 
    best primer sets from the previous iteration.

    Args:
        tasks: The list of primers to consider adding to top_set.
        top_set: One of the best primer sets from the previous which we want to build upon. 
        top_fg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand, respectively, of the on-target genome to any of the primers in top set, where the index corresponds to the h5py path prefix.
        top_bg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand, respectively, of the off-target genome to any of the primers in top set, where the index corresponds to the h5py path prefix.
        fg_fname_prefixes: The path prefixes to all the h5py files relevant to the on-target genome.
        bg_fname_prefixes: The path prefixes to all the h5py files relevant to the off-target genome.
        fg_seq_lengths: The lengths of the on-target genomes in the fasta files (in the same order as fg_fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes in the fasta files (in the same order as bg_fname_prefixes).
        initial_primer_sets: The primer sets that we want to build off of. If set to None, the program will start from scratch.
        iterations: The number of iterations to try to run the optimization. It will prematurely halt if the scores of the primer sets are no longer increasing. By default set to 10.
        max_sets: The number of top sets to keep in every iteration. Set to 10 by default.
        selection_method: The method by which to pick the next iterations set of best primers. See the function make_selection.
        cache: Dictionary mapping the compressed string representation of the primer set to its score.
        This prevents recomputing the score for the same primer set everytime. 
        drop_indices: A list of iterations which will have a drop out layer at the end. See drop_out_layer. 
        If the indices are too early, it may be meaningless, but computational time increases as the index grows.
        
    Returns:
        selected_primer_sets: The selected new primer sets.
        selected_scores: The scores of the selected primer sets.
        selected_fg_fname_to_positions: List of lists of tuples of sorted sets containing positions of exact binding
        locations in the forward and reverse strand of the on-target, respectively, where the first index corresponds to the selected
        primer set and the second index corresponds to the h5py path prefix.
        selected_bg_fname_to_positions: List of lists of tuples of sorted sets containing positions of exact binding
        locations in the forward and reverse strand of the off-target , respectively, where the first index corresponds to the selected
        primer set and the second index corresponds to the h5py path prefix.
        cache: Dictionary mapping the compressed string representation of the primer set to its score.
        This prevents recomputing the score for the same primer set everytime. May or may not have been modified
        during the function call depending on if a primer set was missing from the cache or not.
    """
    seq_initialized_f = partial(evaluate_pairs_helper, curr_fg_fname_to_positions=top_fg_fname_to_positions, curr_bg_fname_to_positions = top_bg_fname_to_positions, fg_fname_prefixes=fg_fname_prefixes, bg_fname_prefixes=bg_fname_prefixes, fg_seq_lengths=fg_seq_lengths, bg_seq_lengths=bg_seq_lengths, fg_circular=fg_circular, bg_circular=bg_circular)
    results = src.utility.create_pool(seq_initialized_f, tasks, src.parameter.cpus)

    selection=make_selection([[task] + top_set for task in tasks], [score for score, temp_fg, temp_bg in results], min(max_sets, len(tasks)), normalize_metric)

    for i, result in enumerate(results):
        primer_set = sorted([tasks[i]] + top_set)
        compressed_string = get_compressed_string(primer_set)
        cache[compressed_string] = results[i][0]

    selected_primer_sets = [sorted([tasks[i]] + top_set) for i in selection]
    selected_scores = [results[i][0] for i in selection]
    selected_fg_fname_to_positions =  [results[i][1] for i in selection]
    selected_bg_fname_to_positions = [results[i][2] for i in selection]

    return selected_primer_sets, selected_scores, selected_fg_fname_to_positions, selected_bg_fname_to_positions, cache 

def make_selection(all_sets, all_scores, max_sets, normalize_metric):
    """
    Chooses max_sets number of sets from all_sets probabilistically based on the scores.

    Args:
       all_sets: A list of all the primer sets to be considered.
       all_scores: The corresponding scores of all the primer sets (same order as in all_sets)
       max_sets: The maximum number of sets to choose.
       normalize_metric: The method be which to choose the sets.
            "deterministic": strictly chooses the largest based on scores
            "softmax": Uses the softmax function to compute probabilities and picks randomly according to these probabilities.
            "normalized": Normalizes the scores by the total sum of the scores and picks randomly according to these probabilities.

    Returns:
        The indices corresponding to all_sets of the sets that were chosen.
    """
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

def get_features_simplified(fname_to_positions, circular, seq_lengths, target=True):
    """
    A function which computes various features for regression.

    Args:
        fname_to_positions: List of lists of tuples of two sorted sets (one for the forward strand and the second for
        the reverse strand) containing the positions of exact binding locations, where the primary of the index corresponds
        to which h5py path prefix it comes from.
        seq_lengths: The lengths of the genomes correspond to the h5py path prefix (in the same order as fname_to_positions.
        target: Boolean indicating whether fname_to_positions corresponds to the on-target genome. We would usually
        not do this for the off-target so default is True.

    Returns:
        data: Dictionary with 'mean_gap' (the mean gap length averaged across individual strands first, then across forward and reverse,
        and then across genomes), 'sum' (The total exact matches across all genomes scaled by the total length of all the fasta files),
        'gap_gini' (The is the gini index computed from each individual strand, then averaged across strands, and then averaged across genomes),
        'within_mean_gap' (This is the average alternating gap mean. See the function get_positional_gap_lengths_alternating),
        'agnostic_mean_gap' (This is mean gap length using all positions (regardless of which strand its on. This is then averaged across genomes).
    """
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

        data_forward = get_position_features(positions_forward, circular, seq_lengths[i])
        data_reverse = get_position_features(positions_reverse, circular, seq_lengths[i])

        gap_sizes.append((data_forward['mean'] + data_reverse['mean'])/seq_lengths[i]/2)
        gap_ginis.append((data_forward['gini'] + data_reverse['gini'])/2)

        alternating_gaps_data = get_positional_gap_lengths_alternating(positions_forward, positions_reverse, circular)
        within_mean_gaps.append(alternating_gaps_data['within_mean']/seq_lengths[i])

        agnostic_gaps_data = get_positional_gap_lengths_agnostic(positions_forward, positions_reverse, circular, seq_length=seq_lengths[i])
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

def get_features(fname_to_positions, circular, seq_lengths, target=True):
    """
    A function which computes various features for regression.

    Args:
        fname_to_positions: List of lists of tuples of two sorted sets (one for the forward strand and the second for
        the reverse strand) containing the positions of exact binding locations, where the primary of the index corresponds
        to which h5py path prefix it comes from.
        seq_lengths: The lengths of the genomes correspond to the h5py path prefix (in the same order as fname_to_positions.
        target: Boolean indicating whether fname_to_positions corresponds to the on-target genome. We would usually
        not do this for the off-target so default is True.

    Returns:
        data: Dictionary with 'mean_gap' (the mean gap length averaged across individual strands first, then across forward and reverse,
        and then across genomes), 'sum' (The total exact matches across all genomes scaled by the total length of all the fasta files),
        'gap_gini' (The is the gini index computed from each individual strand, then averaged across strands, and then averaged across genomes),
        'within_mean_gap' (This is the average alternating gap mean. See the function get_positional_gap_lengths_alternating),
        'agnostic_mean_gap' (This is mean gap length using all positions (regardless of which strand its on. This is then averaged across genomes),
        'agnostic_mean_gini' (This is mean gini index of the gap length using all positions (regardless of which strand its on, averaged across genomes),
        'agnostic_mean_entropy' (This is mean entropy of the gap length using all positions (regardless of which strand its on, averaged across genomes).
    """
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

        data_forward = get_position_features(positions_forward, circular, seq_lengths[i])
        data_reverse = get_position_features(positions_reverse, circular, seq_lengths[i])

        gap_sizes.append((data_forward['mean'] + data_reverse['mean'])/seq_lengths[i]/2)
        gap_ginis.append((data_forward['gini'] + data_reverse['gini'])/2)

        alternating_gaps_data = get_positional_gap_lengths_alternating(positions_forward, positions_reverse, circular, seq_length=seq_lengths[i])
        within_mean_gaps.append(alternating_gaps_data['within_mean']/seq_lengths[i])

        agnostic_gaps_data = get_positional_gap_lengths_agnostic(positions_forward, positions_reverse, circular, seq_length=seq_lengths[i], include_gini=True)

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

def get_positional_gap_lengths_alternating(positions_forward, positions_reverse, circular, threshold=100000):
    """
    A helper function that gets all the lengths of the gaps between adjacent positions on opposite strands which are
    within a threshold. The idea is that this should be more faithful to the actual process of exponential amplification
    which needs a binding site on both strands. This is because a primer needs to be able to amplify the the copies
    of the 5' to 3' strand or vice versa.

    Args:
        positions_forward: A list of all the positions of exact binding locations on the forward strand.
        positions_reverse: A list of all the positions of exact binding locations on the reverse strand.
        threshold: If the gap length between adjacent positions is longer than this value, we assume the phi29 enzyme
        will not elongate this far, and this gap is not added. By default this is 100000.

    Returns:
        data: Dictionary where key 'within_mean' is the mean gap length between adjacent positions on opposite strands.
    """
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

def get_positional_gap_lengths_agnostic(positions_forward, positions_reverse, circular, seq_length=None, include_gini=False):
    """
    A helper function that gets all the lengths of the gaps between adjacent positions, without regard for which strand it occurs on.

    Args:
        positions_forward: A list of all the positions of exact binding locations on the forward strand.
        positions_reverse: A list of all the positions of exact binding locations on the reverse strand.
        seq_length: The length of the genome. By default None.
        include_gini: Boolean indicating whether to compute the gini index and entropy of the gap lengths. False by default.

    Returns:
        data: Dictionary where key 'mean' is the mean gap length, 'gini' is the gini index of the gap lengths,
        and 'entropy' is the entropy of the gap lengths.
    """
    positions = SortedSet(positions_forward)
    positions.update(positions_reverse)

    position_differences = get_positional_gap_lengths(positions,circular, seq_length = seq_length)

    data = {}

    if len(position_differences) <= 1:
        return {'mean': seq_length, 'entropy': 1, 'gini': 1}
    else:
        data['mean'] = np.mean(position_differences)
        if include_gini:
            data['gini'] = src.utility.gini_exact(positions)
            data['entropy'] = scipy.stats.entropy(positions)
    return data


def get_positional_gap_lengths(positions, circular, seq_length=None):
    """
    Gets all the lengths of the gaps between adjacent positions.

    Args:
        positions: A list of all the positions of exact binding locations.
        seq_length: The length of the genome.

    Returns:
        differences: A list of all the gap lengths.
    """
    differences = [a_i - b_i for a_i, b_i in zip(positions[1:], positions[:-1])]

    if differences == []:
        return positions

    if circular:
        differences.append(seq_length - positions[-1] + positions[0])
    else:
        differences.append(seq_length-positions[-1])
        differences.append(positions[0])

    if seq_length is None:
        print("Where is the seq length in get_positional_gap_lengths?")

    return differences

def get_position_features(positions, circular, seq_length):
    """
    Gets the features derived from the positions of the exact binding locations. Namely, The mean gap length and the
    gini index of the positional gap lengths.

    Args:
        positions: A list of all the positions of exact binding locations.
        seq_length: The length of the genome.

    Returns:
        data: Dictionary where key 'mean' is the mean gap length and 'gini' is the gini index of the gap lengths.
    """
    position_differences = get_positional_gap_lengths(positions, circular, seq_length=seq_length)

    data = {}

    if len(position_differences) <= 1:
        return {'mean': seq_length, 'gini': 1}
    else:
        data['mean'] = np.mean(position_differences)
        data['gini'] = src.utility.gini_exact(position_differences)

    return data

def get_exact_count_features_one_fname(primer_set, fname_prefix):
    """
    For each primer in primer_set, this helper function gets the positions of exact binding (all bases in the primer
    and genome bind).

    Args:
        primer_set: The list of primer sequences of which to get the positions for.
        fname_prefix: The h5py path prefix for the genome, basically the path minus '_6mer_positions.h5' where k = 6.

    Returns:
        curr_pos_results_forward: A list of exact binding locations on the forward strand.
        curr_pos_results_reverse: A list of exact binding locations on the reverse strand.
    """
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

# def get_fname_to_positions_dict(fname_prefixes, primer_list):
#
#     fname_to_positions = {}
#
#     for fname_prefix in fname_prefixes:
#         positions_forward, positions_reverse = get_exact_count_features_one_fname(primer_list, fname_prefix)
#         fname_to_positions[fname_prefix] = (positions_forward, positions_reverse)
#         # matches += len(positions)
#     return fname_to_positions

def incorporate_positions_row(curr_fname_to_positions, fname_prefixes, primer):
    """
    Incorporates the set of positions that a primer may exactly bind to the genome.

    Args:
        curr_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand, respectively, of the genome to any of the primers in primer set, where the index corresponds to the h5py path prefix.
        fname_prefixes: The h5py path prefixes for the genomes of interest.
        primer: The primer of interest.

    Returns:
        fname_to_positions:  List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand, respectively, updated with the exact binding
        positions of 'primer', where the index corresponds to the h5py path prefix.
        total_matches: The total number of exact matches in the genome with the primer (on the forward and reverse strand).
    """
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

def evaluate(temp_fg_fname_to_positions, temp_bg_fname_to_positions, fg_seq_lengths=None, bg_seq_lengths=None, fg_circular=True, bg_circular=False, return_X=False):
    """
    Predicts the amplification score of the entire primer set using a pretrained ridge regression model. Current features are:
        'ratio': The ratio of total exact binding locations in the on-target vs. the off-target. Each total is scaled by their total genome lengths first.
        'agnostic_mean_gap_ratio': This ratio between that of the on-target to off-target of the mean gap lengths
        using all positions (regardless of which strand its on. This is then averaged across genomes
        'on_gap_gini': The is the gini index computed from each individual strand, then averaged across strands, and then averaged across the on-target genomes.
        'off_gap_gini': The is the gini index computed from each individual strand, then averaged across strands, and then averaged across off-target genomes.
        'within_mean_gap_ratio': The ratio of the within_mean_gaps of the OFF-target to ON-target, where the
        within_mean_gap is defined as the is the average alternating gap mean. See the function get_positional_gap_lengths_alternating.
        The ratio is not computed as ON-target to OFF-target because we want this to be 0 in the off-target which is
        numerically unstable.

    Args:
        temp_fg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the on-target
        forward and reverse strand, to any of the primers in primer set, where the index corresponds to the h5py path prefix.
        temp_bg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the off-target
        forward and reverse strand, to any of the primers in primer set, where the index corresponds to the h5py path prefix.
        fg_seq_lengths: The lengths of the on-target genomes corresponding to the h5py files (in the same order as fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes corresponding to the h5py files (in the same order as fname_prefixes).
        return_X: Boolean indicating whether to return the feature matrix or not.
        fg_circular: Boolean variable indicating whether the on-target genome is circular. Set to true by default.
        bg_circular: Boolean variable indicating whether the off-target genome is circular. Set to true by default.

    Returns:
        score: The predicted score of he primer set.
        X (optional): The feature matrix to be optionally returned.
    """
    warnings.filterwarnings('ignore') 

    fg_features = get_features_simplified(temp_fg_fname_to_positions, fg_circular, fg_seq_lengths, target=True)
    bg_features = get_features_simplified(temp_bg_fname_to_positions, bg_circular, bg_seq_lengths, target=False)

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

    X = all_features[['ratio', 'agnostic_mean_gap_ratio', 'on_gap_gini', 'off_gap_gini', 'within_mean_gap_ratio']]
    
    fo = open(os.path.join(src.parameter.src_dir, 'ratio_agnostic_mean_gap_ratio_all_ginis_within_mean_gap_ratio.p'),
        'rb')
    clf = pickle.load(fo)
    fo.close()


    score = clf.predict(X)[0]

    if return_X:
        return score, X

    return score

def evaluate_pairs_helper(new_primer, curr_fg_fname_to_positions=None, curr_bg_fname_to_positions=None, fg_fname_prefixes=None, bg_fname_prefixes=None, fg_seq_lengths=None, bg_seq_lengths=None, fg_circular=True, bg_circular=False):
    """
     Computes the score resulting from adding a new primer to the primer set.

     Args:
        new_primer: The new primer to incorporate into the primer set.
        curr_fg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand of the on-target genome, to any of the primers in primer set, where the index corresponds
        to the h5py path prefix.
        curr_bg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand of the off-target genome, to any of the primers in primer set, where the index corresponds
        to the h5py path prefix.
        fg_prefixes: The path prefix to the h5py files of the on-target genome, basically the path minus '_6mer_positions.h5' where k = 6.
        bg_prefixes: The path prefix to the h5py files of the off-targe genome, basically the path minus '_6mer_positions.h5' where k = 6.
        fg_seq_lengths: The lengths of the on-target genomes corresponding to the h5py files (in the same order as fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes corresponding to the h5py files (in the same order as fname_prefixes).

     Returns:
        score: The predicted amplification score for the primer set and the new primer.
        temp_fg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand, respectively, of the on-target genome to any of the primers in primer set, where the index corresponds to the h5py path prefix.
        temp_bg_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand, respectively, of the off-target genome to any of the primers in primer set, where the index corresponds to the h5py path prefix.
     """
    temp_fg_fname_to_positions, _ = incorporate_positions_row(curr_fg_fname_to_positions, fg_fname_prefixes, new_primer)
    temp_bg_fname_to_positions, _ = incorporate_positions_row(curr_bg_fname_to_positions, bg_fname_prefixes, new_primer)
    score = evaluate(temp_fg_fname_to_positions, temp_bg_fname_to_positions, fg_seq_lengths=fg_seq_lengths, bg_seq_lengths=bg_seq_lengths, fg_circular=fg_circular, bg_circular=bg_circular)
    return score, temp_fg_fname_to_positions, temp_bg_fname_to_positions

def build_fname_to_positions(primer_set, fname_prefixes=None, seq_lengths=None, circular=True):
    """
    This is a helper function which builds a dictionary of h5py path prefixes to the positions of exact binding sites
    between any primer in primer_set and the genomes corresponding to the h5py files.

    Args:
        primer_set: The primer set of interest.
        fname_prefixes: The path prefix to the h5py files, basically the path minus '_6mer_positions.h5' where k = 6.
        seq_lengths: The lengths of the genomes corresponding to the h5py files (in the same order as fname_prefixes).
        circular: Boolean variable indicating whether the genome is circular. Set to true by default.

    Returns:
        temp_fname_to_positions: List of tuples of sets containing positions of exact binding locations in the forward
        and reverse strand, respectively, to any of the primers in primer set, where the index corresponds to the h5py path prefix.
    """
    temp_fname_to_positions = initialize_fname_to_positions_dict(fname_prefixes, seq_lengths)
    for primer in primer_set:
        temp_fname_to_positions, total_matches = incorporate_positions_row(temp_fname_to_positions, fname_prefixes, primer)
    return temp_fname_to_positions

def evaluate_wrapper(fname_to_positions, fg_seq_lengths=None, bg_seq_lengths=None, fg_circular=True, bg_circular=False):
    """
    Simple helper function that calls evaluate for a particular primer set.

    Args:
        fname_to_positions: List of tuples of positions for the foreground and background, respectively.
        fg_seq_lengths: The lengths of the on-target genomes corresponding to the h5py files (in the same order as fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes corresponding to the h5py files (in the same order as fname_prefixes).

    Returns: returns the output of evaluate.

    """
    return evaluate(fname_to_positions[0], fname_to_positions[1], fg_seq_lengths=fg_seq_lengths,bg_seq_lengths=bg_seq_lengths, fg_circular=fg_circular, bg_circular=bg_circular)

def stand_alone_score_primer_sets(primer_sets, fg_prefixes, bg_prefixes, fg_seq_lengths, bg_seq_lengths, fg_circular, bg_circular):
    """
    This is a function to quickly evaluate a list of candidate primer sets.

    Args:
        primer_sets: The primer sets to evaluate.
        fg_prefixes: The path prefix to the h5py files of the on-target genome, basically the path minus '_6mer_positions.h5' where k = 6.
        bg_prefixes: The path prefix to the h5py files of the off-targe genome, basically the path minus '_6mer_positions.h5' where k = 6.
        fg_seq_lengths: The lengths of the on-target genomes corresponding to the h5py files (in the same order as fname_prefixes).
        bg_seq_lengths: The lengths of the off-target genomes corresponding to the h5py files (in the same order as fname_prefixes).
        fg_circular: Boolean variable indicating whether the on-target genome is circular. Set to true by default.
        bg_circular: Boolean variable indicating whether the off-target genome is circular. Set to true by default.

    Returns:
        top_scores: The scores of the evaluated primer sets.
    """
    fg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=fg_prefixes,
                                   seq_lengths=fg_seq_lengths, circular=fg_circular)
    initial_fg_fname_to_positions = src.utility.create_pool(fg_seq_initialized_f, primer_sets, cpus=src.parameter.cpus)

    bg_seq_initialized_f = partial(build_fname_to_positions, fname_prefixes=bg_prefixes,
                                   seq_lengths=bg_seq_lengths, circular=bg_circular)
    initial_bg_fname_to_positions = src.utility.create_pool(bg_seq_initialized_f, primer_sets, cpus=src.parameter.cpus)

    partial_initialized_f = partial(evaluate_wrapper, fg_seq_lengths=fg_seq_lengths,
                                    bg_seq_lengths=bg_seq_lengths, fg_circular=fg_circular, bg_circular=bg_circular)
    top_scores = src.utility.create_pool(partial_initialized_f,
                                     list(zip(initial_fg_fname_to_positions, initial_bg_fname_to_positions)), cpus=src.parameter.cpus)
    print(primer_sets)
    print(top_scores)
    return top_scores

if __name__ == "__main__":
    get_position_features([3, 8], True, 9)
