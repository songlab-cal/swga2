import warnings
from optparse import OptionParser
import argparse
import src.parameter
import src.utility
import multiprocessing
import json
import sys
import pandas as pd
import src.kmer
import src.rf_preprocessing
import src.utility
import pickle
import src.optimize
import numpy as np
import src.filter
import src.string_search
import os
import warnings


defaults = {
    "min_fg_freq": float(1 / 100000),
    "max_bg_freq": float(1 / 200000),
    "max_gini": 0.6,
    "max_primer": 500,
    "min_amp_pred": 10,
    "min_tm": 15,
    "max_tm": 45,
    "max_dimer_bp": 3,
    "max_self_dimer_bp": 4,
    "selection_metric": "deterministic",
    "iterations": 8,
    "top_set_count": 10,
    "retries": 5,
    "max_sets": 5,
    "fg_circular": True,
    "bg_circular": False,
    "drop_iterations": 5,
    "verbose": False,
    "cpus": int(multiprocessing.cpu_count()),
}

parser = OptionParser()
parser.add_option("-k", "--kmer_fore", help="")
parser.add_option("-l", "--kmer_back", help="")
parser.add_option("-x", "--fasta_fore", help="")
parser.add_option("-y", "--fasta_back", help="")
parser.add_option("-j", "--json-file", type=str, help="")
parser.add_option("-z", "--data_dir", type=str, help="")
parser.add_option("-!", "--src_dir", type=str, help="")
parser.add_option("-p", "--min_fg_freq", type=float, help="")
parser.add_option("-q", "--max_bg_freq", type=float, help="")
parser.add_option("-g", "--max_gini", type=float, help="")
parser.add_option("-e", "--max_primer", type=int, help="")
parser.add_option("-a", "--min_amp_pred", type=float, help="")
parser.add_option("-m", "--min_tm", type=int, help="")
parser.add_option("-n", "--max_tm", type=int, help="")
parser.add_option("-t", "--max_dimer_bp", type=int, help="")
parser.add_option("-u", "--max_self_dimer_bp", type=int, help="")
parser.add_option("-s", "--selection_metric", type=str, help="")
parser.add_option("-i", "--iterations", type=int, help="")
parser.add_option("-#", "--top_set_count", type=int, help="")
parser.add_option("-r", "--retries", type=int, help="")
parser.add_option("-w", "--max_sets", type=int, help="")
parser.add_option("-f", "--fg_circular", help="")
parser.add_option("-b", "--bg_circular", help="")
parser.add_option("-d", "--drop_iterations", type=int, help="")
parser.add_option("-v", "--verbose", help="")
parser.add_option("-c", "--cpus", type=int, help="")
(options, args) = parser.parse_args()
params = src.parameter.get_params(options)

fg_prefixes = params["fg_prefixes"]
bg_prefixes = params["bg_prefixes"]

fg_genomes = params["fg_genomes"]
bg_genomes = params["bg_genomes"]

fg_seq_lengths = params["fg_seq_lengths"]
bg_seq_lengths = params["bg_seq_lengths"]

fg_circular = params["fg_circular"]
bg_circular = params["bg_circular"]


def main():

    arg_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)

    arg_parser.add_argument("command", type=str, choices=["step1", "step2", "step3", "step4"])

    args, remaining = arg_parser.parse_known_args()

    with open(params["json_file"], "r") as infile:
        json_params = json.load(infile)

    # command line params will override json params
    for param, param_val in params.items():
        if param_val is not None:
            json_params[param] = param_val

    # default params will be added if not included in the json params or the command line params
    for default_param, default_val in defaults.items():
        if default_param not in json_params or json_params[default_param] is None:
            json_params[default_param] = default_val

    with open(params["json_file"], "w+") as outfile:
        json.dump(json_params, outfile, indent=4)

    print("======================================================================================================")
    print("Parameters:")
    for param_key, param_val in params.items():
        print(param_key + ": " + str(param_val))
    print("======================================================================================================")

    if args.command == "step1":
        step1()
    elif args.command == "step2":
        step2()
    elif args.command == "step3":
        step3()
    elif args.command == "step4":
        step4()
    else:
        print("Please input a command: step1, step2, step3, or step4.")


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

    print("Done running jellyfish")

    return fg_seq_lengths, bg_seq_lengths


def step2(all_primers=None):
    """
    Filters all candidate primers according to primer design principles (http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html)
    and minimum foreground frequency, maximum background frequency, and maximum Gini index of distances between binding sites.

    Args:
        all_primers: The list of candidate primers to consider. Defaults to all k-mers from 6 to 12 of the target genome if not input.
        circular: Boolean variable indicating whether the genome is circular. Set to true by default.

    Returns:
        filtered_gini_df: Pandas dataframe containing sequences which pass the foreground frequency, background frequency, and Gini index filters.
    """
    kwargs = {
        "fg_prefixes": fg_prefixes,
        "bg_prefixes": bg_prefixes,
        "fg_total_length": sum(fg_seq_lengths),
        "bg_total_length": sum(bg_seq_lengths),
    }
    if all_primers is None:
        all_primers = src.kmer.get_primer_list_from_kmers(fg_prefixes)
    print("Computing foreground and background rates...")
    rate_df = src.filter.get_all_rates(all_primers, **kwargs)
    # print(rate_df)
    filtered_rate_df = rate_df[(rate_df["fg_bool"]) & (rate_df["bg_bool"])]
    # print(filtered_rate_df)
    filtered_rate_df = filtered_rate_df.drop(["fg_bool", "bg_bool"], axis=1)
    print(
        "Filtered "
        + str(len(rate_df) - len(filtered_rate_df))
        + " number of primers based on foreground/background rate."
    )

    if src.parameter.verbose:
        print(rate_df)

    print("Computing Gini index...")
    # print(src.parameter.fg_circular)
    gini_df = src.filter.get_gini(fg_prefixes, fg_genomes, fg_seq_lengths, filtered_rate_df, fg_circular)
    print("Filtered " + str(len(filtered_rate_df) - len(gini_df)) + " number of primers based on gini rate.")
    gini_df["ratio"] = gini_df["bg_count"] / gini_df["fg_count"]
    filtered_gini_df = gini_df.sort_values(by=["ratio"], ascending=False)[: src.parameter.max_primer]

    # pickle.dump(filtered_gini_df, open(os.path.join(src.parameter.data_dir, 'step2_df.p'), 'wb'))
    filtered_gini_df.to_csv(os.path.join(src.parameter.data_dir, "step2_df.csv"))
    src.string_search.get_positions(
        filtered_gini_df["primer"], fg_prefixes, fg_genomes, circular=src.parameter.fg_circular
    )
    src.string_search.get_positions(
        filtered_gini_df["primer"], bg_prefixes, bg_genomes, circular=src.parameter.bg_circular
    )
    print("Number of remaining primers: " + str(len(filtered_gini_df["primer"])))
    # print(filtered_gini_df['primer'])
    return filtered_gini_df


# RANK BY RANDOM FOREST
def step3():
    """
    Filters primers according to primer efficacy. To adjust the threshold, use option -a or --min_amp_pred.

    Returns:
        joined_step3_df: Pandas dataframe of sequences passing step 3.
    """
    step2_df = pd.read_csv(os.path.join(src.parameter.data_dir, "step2_df.csv"))
    # step2_df = pickle.load(open(os.path.join(src.parameter.data_dir, 'step2_df.p'), 'rb'))

    primer_list = step2_df["primer"]
    print("Number of primers beginning of step 3: " + str(len(primer_list)))
    fg_scale = sum(fg_seq_lengths) / 6200

    print("Computing random forest features...")
    df_pred = src.rf_preprocessing.create_augmented_df(fg_prefixes, primer_list)
    df_pred = src.rf_preprocessing.scale_delta_Gs(df_pred, on_scale=fg_scale)
    df_pred["molarity"] = 2.5

    print("Computing amplification efficacy...")
    results = src.rf_preprocessing.predict_new_primers(df_pred)
    results.sort_values(by=["on.target.pred"], ascending=[False], inplace=True)

    step3_df = results[results["on.target.pred"] >= src.parameter.min_amp_pred]

    step2_df = step2_df.set_index("primer")
    step3_df = step3_df.rename({"sequence": "primer"}, axis="columns")
    step3_df = step3_df.set_index("primer")

    joined_step3_df = step3_df.join(step2_df[["ratio", "gini", "fg_count", "bg_count"]], how="left").sort_values(
        by="gini"
    )

    # pickle.dump(joined_step3_df, open(os.path.join(src.parameter.data_dir, 'step3_df.p'), "wb"))
    joined_step3_df.to_csv(os.path.join(src.parameter.data_dir, "step3_df.csv"))

    print("Filtered " + str(step2_df.shape[0] - joined_step3_df.shape[0]) + " number of primers based on efficacy.")

    if src.parameter.verbose:
        print(joined_step3_df)

    return joined_step3_df


# src.optimize AND SEARCH
def step4(primer_list=None, scores=None, initial_primer_sets=None):
    """Searches for primer sets using candidate primers from step 3.

    Args:
        primer_list: The list of primers from which to build sets. Otherwise, the lists defaults to all primers from step  3.
        scores: The scores or probabilities to use for initial randomization. Defaults the primer efficacy scores from step 3.
        initial_primer_sets: The list of initialized primer sets. Defaults to random selection according to normalized scores.
    """
    if primer_list is None:
        # step3_df = pickle.load(open(os.path.join(src.parameter.data_dir, 'step3_df.p'), 'rb'))
        step3_df = pd.read_csv(os.path.join(src.parameter.data_dir, "step3_df.csv"))
        primer_list = step3_df["primer"]
        scores = step3_df["on.target.pred"]

    all_results = []
    all_scores = []
    cache = {}
    banned_primers = []

    for i in range(src.parameter.retries):
        print("Repeat #: " + str(i + 1))

        primer_list = list(set(primer_list) - set(banned_primers))

        if initial_primer_sets is not None:
            primer_list.extend(src.utility.flatten(initial_primer_sets))
        else:
            initial_primer_sets = src.optimize.random_initial_start(
                np.asarray(primer_list).reshape((-1, 1)), scores, max_sets=src.parameter.max_sets
            )

        print("Initial primers: ")
        print(initial_primer_sets)

        primer_list = sorted(list(set(primer_list)))
        random_primers = src.utility.flatten(initial_primer_sets)
        combined_primer_list = list(set(primer_list + random_primers))

        print("Banned primers: " + str(banned_primers))
        results, scores, cache = src.optimize.bfs(
            combined_primer_list,
            fg_prefixes,
            bg_prefixes,
            fg_seq_lengths,
            bg_seq_lengths,
            initial_primer_sets=initial_primer_sets,
            iterations=src.parameter.iterations,
            max_sets=src.parameter.max_sets,
            selection_method=src.parameter.selection_metric,
            banned_primers=banned_primers,
            fg_circular=src.parameter.fg_circular,
            bg_circular=src.parameter.bg_circular,
            cache=cache,
            drop_indices=src.parameter.drop_iterations,
        )

        all_results.append(results)
        all_scores.append(scores)

        banned_primers.append(src.utility.most_frequent(src.utility.flatten(results)))
        banned_primers.append(src.utility.most_frequent(src.utility.flatten(results)))
        banned_primers = list(set(banned_primers))

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

    print("FINAL TOP SETS FROM ALL RESTARTS: " + str(src.parameter.top_set_count))
    for i, final_set in enumerate(final_sets[: src.parameter.top_set_count]):
        print("[" + ", ".join(map(str, final_set)) + "], " + str(round(final_scores[i], 2)))


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
    # step4()
    # print(parameter.min_fg_freq)

    src.optimize.stand_alone_score_primer_sets(
        [["AAAGTATGG", "ATTACAGCA", "GATAAGAGA", "GTGATGAAA", "TCAACTATA"]],
        fg_prefixes,
        bg_prefixes,
        fg_seq_lengths,
        bg_seq_lengths,
    )
