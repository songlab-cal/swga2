import src.thermo_estimation as te
import src.utility
import pandas as pd
from functools import partial
import src.parameter
import numpy as np
import pickle
import melting
import os
import src.kmer
import sklearn
import warnings

# Names of the features not including the thermodynamicfeatures
base_features = ['molarity', 'sequence.length', 'number.of.A', 'proportion.of.A', 'number.of.T', 'proportion.of.T', 'number.of.G', 'proportion.of.G', 'number.of.C', 'proportion.of.C', 'GC.content', 'melting_tm', 'GC.clamp', 'longest.A.repeat', 'longest.T.repeat', 'longest.G.repeat', 'longest.C.repeat', 'AA repeat', 'CC repeat', 'TT repeat', 'GG repeat', '3.end.first.base', '3.end.second.base', '3.end.third.base', '3.end.fourth.base', '3.end.fifth.base']

# Creates the thermodynamic features for the on target (d_g_on_features)
bins = [-20, -18, -16, -14, -12, -10, -9, -8, -7, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
histogram_upper_bound=bins[-1]
histogram_lower_bound=bins[0]
features = list(base_features)
delta_g_on_features = ['on_target_' + str(bin) for bin in bins[1:]]

# List of all features in the regression
regression_features = features + delta_g_on_features

def get_features(primer, target=None, molarity=2.5):
    """
    Function called by create_base_feature_matrix pool and used to make feature vector for each primer.

    Args:
        primer: The primer for which to make the feature vector for.
        target: Boolean variable--true if the primer exists in the on-target genome (dropped in the regression).
        molarity: The molarity of the plasmid experiment concentration.

    Returns:
        result: A list of all the base feature values, the primer sequence, and whether it is in the on-target genome.
    """
    feature_dict = {}
    feature_dict['molarity'] = molarity
    feature_dict['sequence.length'] = len(primer)

    feature_dict['number.of.A'] = primer.count('A')
    feature_dict['number.of.G'] = primer.count('G')
    feature_dict['number.of.T'] = primer.count('T')
    feature_dict['number.of.C'] = primer.count('C')

    feature_dict['proportion.of.A'] = feature_dict['number.of.A']/feature_dict['sequence.length']
    feature_dict['proportion.of.G'] = feature_dict['number.of.G'] / feature_dict['sequence.length']
    feature_dict['proportion.of.T'] = feature_dict['number.of.T'] / feature_dict['sequence.length']
    feature_dict['proportion.of.C'] = feature_dict['number.of.C'] / feature_dict['sequence.length']

    feature_dict['GC.content'] = (feature_dict['number.of.G'] + feature_dict['number.of.C'])/ feature_dict['sequence.length']

    feature_dict['longest.A.repeat'] =  src.utility.longest_char_repeat(primer,'A')
    feature_dict['longest.G.repeat'] =  src.utility.longest_char_repeat(primer, 'G')
    feature_dict['longest.T.repeat'] =  src.utility.longest_char_repeat(primer, 'T')
    feature_dict['longest.C.repeat'] =  src.utility.longest_char_repeat(primer, 'C')

    feature_dict['AA repeat'] = primer.count('AA')
    feature_dict['GG repeat'] = primer.count('GG')
    feature_dict['CC repeat'] = primer.count('CC')
    feature_dict['TT repeat'] = primer.count('TT')

    feature_dict['target'] = target

    feature_dict['melting_tm'] = melting.temp(primer)
    last_five_end = primer[:-5]
    feature_dict['GC.clamp'] = last_five_end.count('C') + last_five_end.count('G')

    feature_dict['3.end.first.base'] =  src.utility.char_to_int_dict[primer[-1]]
    feature_dict['3.end.second.base'] =  src.utility.char_to_int_dict[primer[-2]]
    feature_dict['3.end.third.base'] =  src.utility.char_to_int_dict[primer[-3]]
    feature_dict['3.end.fourth.base'] =  src.utility.char_to_int_dict[primer[-4]]
    feature_dict['3.end.fifth.base'] =  src.utility.char_to_int_dict[primer[-5]]

    result = [primer, target]
    for feature in base_features:
        result.append(feature_dict[feature])
    return result

def create_base_feature_matrix(primer_list, molarity):
    """
    Creates the feature matrix, not including the delta_G features.

    Args:
        primer_list: The list of primers to create the feature matrix for.
        molarity: The molarity of the plasmid experiment experiments.

    Returns:
        df: The pandas dataframe including the features for the base regression features, the sequence of the primer,
        and whether it exists in the target genome (where the last two features are dropped in regression).
    """
    get_features_partial = partial(get_features, molarity = molarity)
    results =  src.utility.create_pool(get_features_partial, primer_list, src.parameter.cpus)
    f = ['sequence','target']
    f.extend(base_features)
    df = pd.DataFrame(results, columns=f)
    return df

def get_all_predicted_delta_G_for_all_files_transformed(primer_list, target, fnames, f_out_name=None):
    """
    Gets all the predicted thermodynamic values for all the candidate primers (primer_list) in all the genome files.

    Args:
        primer_list: The list of candidate primers which we are evaluating.
        target: Boolean variable--true if the primer exists in the on-target genome (dropped in the regression).
        fnames: A list of the paths to the genome files.
        f_out_name: If this is not None, this saves the feature matrix using this as the file path.

    Returns:
        df: A pandas dataframe representing the feature matrix of the basic features (defined at the top) for the
        primers in primer_list.

    """
    seq_initialized_f = partial(get_all_predicted_delta_G_per_primer_transformed, fnames=fnames)

    arr_transformed =  src.utility.create_pool(seq_initialized_f, primer_list, src.parameter.cpus)           #this is the results, its in the form of list of lists of length 2 containing

    if target:
        columns = ["on_target_" + str(i) for i in bins[1:]]
    else:
        columns = ["off_target_" + str(i) for i in bins[1:]]

    df = pd.DataFrame(arr_transformed, columns=columns)

    if f_out_name:
        pickle.dump(df, open(f_out_name, "wb"))

    return df

def get_all_predicted_delta_G_per_primer_transformed(primer, fnames=None, penalty=4):
    """
    This function takes a primer and computes a vector of the binned thermodynamic values. It takes account both the
    forward and the reverse strand. For computing potential binding to k-mers in the forward strand, we want the
    thermodynamic values between the k-mer in the genome written 5' to 3' and the primer (written 3' to 5').
    For computing potential binding to k-mers in the reverse strand, we want the thermodynamic values between
    the primer and the complement of the k-mer in the genome.

    Args:
        primer: The primer for which to compute the predicted thermodynamic values.
        fnames: A list of the genomes to search and compute the predicted thermodynamic values binding with the primer.
        penalty: A value to add to the thermodynamic sum if a pair does not exist in Santa-Lucia's table.
        This usually means everything is a mismatch and it makes the binding much less likely.

    Returns:
        all_delta_G_vals: Vector of the histogram frequencies from -20 to 3 of the thermodynamic values resulting
        between the primer and all possible kmers in the genome files.
    """
    k = len(primer)

    all_delta_G_vals = [0 for i in range(len(bins)-1)]

    for fname_prefix in fnames:
        if os.path.exists(fname_prefix + '_' + str(k) + 'mer_all.txt'):
            for kmer_i, count in src.kmer.get_kmer_to_count_dict(fname_prefix + '_' + str(k) + 'mer_all.txt').items():
                delta_G_val_forward = te.compute_free_energy_for_two_strings(kmer_i,  src.utility.reverse(primer), penalty)
                delta_G_val_reverse = te.compute_free_energy_for_two_strings(primer,  src.utility.complement(kmer_i), penalty)  # this is to search the backward strand, but we don't need to reverse it

                if delta_G_val_forward <= histogram_upper_bound:
                    if delta_G_val_forward == histogram_upper_bound:
                        i = len(bins) - 1
                    else:
                        i = np.digitize(delta_G_val_forward, bins)
                    all_delta_G_vals[i-1] += count

                if delta_G_val_reverse <= histogram_upper_bound:
                    if delta_G_val_reverse == histogram_upper_bound:
                        i = len(bins) - 1
                    else:
                        i = np.digitize(delta_G_val_reverse, bins)
                    all_delta_G_vals[i-1] += count
    return all_delta_G_vals

def create_augmented_df(fg_fnames, primer_list=None, feature_matrix=None, delta_G_matrix=None, molarity=None):
    """
    Creates the feature matrix including the thermodynamic features.

    Args:
        fg_fnames: A list of paths to the genome files.
        primer_list: The list of primers to create the feature matrix for. If None the feature_matrix cannot be None.
        feature_matrix: The feature matrix of the base features built using the above function create_base_feature_matrix.
        delta_G_matrix: A matrix containing the histogram frequences of the thermodynamic features, assuming it was
        computed already.
        molarity: The molarity of the plasmid experiment experiments.

    Returns:
        df: The pandas dataframe of the regression features (including the thermodynamic features),
        the sequence of the primer, and whether it exists in the target genome (where the last two features are
        dropped in regression).
    """
    if feature_matrix is None:
        print("feature matrix was none")
        feature_matrix = create_base_feature_matrix(primer_list, molarity=molarity)
    feature_matrix = feature_matrix.reset_index()
    if delta_G_matrix is None:
        print("delta G matrix was none")
        if primer_list is None:
            primer_list = feature_matrix['sequence']

        delta_G_matrix = get_all_predicted_delta_G_for_all_files_transformed(primer_list, True, fg_fnames)
        print("Done with fg")

    df = pd.concat([feature_matrix, delta_G_matrix], axis=1)
    return df

def scale_delta_Gs(df_pred, on_scale=4000):
    """
    This scales the histogram frequencies of the thermodynamic values. This is necessary because often the off-target
    genome is much larger than the on-target genome. Using the ratio of the on-target length to off-target length as
    the on_scale is recommended.

    Args:
        df_pred: The dataframe containing the histogram frequencies of the thermodynamic values.
        on_scale: The value by which to normalize each value.

    Returns:
        df_pred: The rescaled dataframe.

    """
    # print(df_pred.columns)
    df_pred[delta_g_on_features] = df_pred[delta_g_on_features]/(on_scale)
    return df_pred

def predict_new_primers(df):
    """
    Runs the random forest regression on the input features to predict amplification scores.

    Args:
        df: The input features.

    Returns:
        output_df: The predicted amplification scores.

    """
    warnings.filterwarnings('ignore') 
    
    X_test = df[regression_features]
    X_test = X_test.dropna(axis='columns')

    clf = pickle.load(open(os.path.join(src.parameter.src_dir, 'random_forest_filter.p'),'rb'))
    y_pred = clf.predict(X_test)
    output_df = pd.DataFrame(y_pred, columns=['on.target.pred'])
    output_df = pd.concat([output_df, df], axis=1)
    return output_df

if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    myco_prefixes = [src.parameter.data_dir + 'kmer_files/myco']

    print(create_augmented_df(myco_prefixes, ['ACAACC','TACGTCA'], molarity=2.5))