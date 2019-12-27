import thermo_estimation as te
import utility
import pandas as pd
from functools import partial
import parameter
import numpy as np
import pickle
import melting
import os
import kmer
import sklearn

base_features = ['molarity', 'sequence.length', 'number.of.A', 'proportion.of.A', 'number.of.T', 'proportion.of.T', 'number.of.G', 'proportion.of.G', 'number.of.C', 'proportion.of.C', 'GC.content', 'melting_tm', 'GC.clamp', 'longest.A.repeat', 'longest.T.repeat', 'longest.G.repeat', 'longest.C.repeat', 'AA repeat', 'CC repeat', 'TT repeat', 'GG repeat', '3.end.first.base', '3.end.second.base', '3.end.third.base', '3.end.fourth.base', '3.end.fifth.base']
bins = [-20, -18, -16, -14, -12, -10, -9, -8, -7, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
histogram_upper_bound=bins[-1]
histogram_lower_bound=bins[0]

features = list(base_features)
delta_g_on_features = ['on_target_' + str(bin) for bin in bins[1:]]
regression_features = features+delta_g_on_features

def get_features(primer, target=None, molarity=2.5):
    '''Called by create_base_feature_matrix pool and used to make feature vector for each primer.'''
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

    feature_dict['longest.A.repeat'] = utility.longest_char_repeat(primer,'A')
    feature_dict['longest.G.repeat'] = utility.longest_char_repeat(primer, 'G')
    feature_dict['longest.T.repeat'] = utility.longest_char_repeat(primer, 'T')
    feature_dict['longest.C.repeat'] = utility.longest_char_repeat(primer, 'C')

    feature_dict['AA repeat'] = primer.count('AA')
    feature_dict['GG repeat'] = primer.count('GG')
    feature_dict['CC repeat'] = primer.count('CC')
    feature_dict['TT repeat'] = primer.count('TT')

    feature_dict['target'] = target

    feature_dict['melting_tm'] = melting.temp(primer)
    last_five_end = primer[:-5]
    feature_dict['GC.clamp'] = last_five_end.count('C') + last_five_end.count('G')

    feature_dict['3.end.first.base'] = utility.char_to_int_dict[primer[-1]]
    feature_dict['3.end.second.base'] = utility.char_to_int_dict[primer[-2]]
    feature_dict['3.end.third.base'] = utility.char_to_int_dict[primer[-3]]
    feature_dict['3.end.fourth.base'] = utility.char_to_int_dict[primer[-4]]
    feature_dict['3.end.fifth.base'] = utility.char_to_int_dict[primer[-5]]

    result = [primer, target]
    for feature in base_features:
        result.append(feature_dict[feature])
    return result

def create_base_feature_matrix(primer_list, molarity):
    '''Creates the feature matrix, not including the delta_G features.'''
    get_features_partial = partial(get_features, molarity = molarity)
    results = utility.create_pool(get_features_partial, primer_list, parameter.cpus)
    f = ['sequence','target']
    f.extend(base_features)
    df = pd.DataFrame(results, columns=f)
    return df

def get_all_predicted_delta_G_for_all_files_transformed(primer_list, target, fnames, f_out_name=None):

    seq_initialized_f = partial(get_all_predicted_delta_G_per_primer_transformed, fnames=fnames)

    arr_transformed = utility.create_pool(seq_initialized_f, primer_list, parameter.cpus)           #this is the results, its in the form of list of lists of length 2 containing

    if target:
        columns = ["on_target_" + str(i) for i in bins[1:]]
    else:
        columns = ["off_target_" + str(i) for i in bins[1:]]

    df = pd.DataFrame(arr_transformed, columns=columns)

    if f_out_name:
        pickle.dump(df, open(f_out_name, "wb"))

    return df

def get_all_predicted_delta_G_per_primer_transformed(primer, fnames=None, penalty=parameter.mismatch_penalty):
    k = len(primer)

    all_delta_G_vals = [0 for i in range(len(bins)-1)]

    for fname_prefix in fnames:
        if os.path.exists(fname_prefix + '_' + str(k) + 'mer_all.txt'):
            for kmer_i, count in kmer.get_kmer_to_count_dict(fname_prefix + '_' + str(k) + 'mer_all.txt').items():
                delta_G_val_forward = te.compute_free_energy_for_two_strings(kmer_i, utility.reverse(primer), penalty)
                delta_G_val_reverse = te.compute_free_energy_for_two_strings(primer, utility.complement(kmer_i), penalty)  # this is to search the backward strand, but we don't need to reverse it

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
    print(df_pred.columns)
    df_pred[delta_g_on_features] = df_pred[delta_g_on_features]/(on_scale)
    return df_pred

def predict_new_primers(df):
    X_test = df[regression_features]
    X_test = X_test.dropna(axis='columns')

    clf = pickle.load(open(os.path.join(parameter.src_dir, 'random_forest_filter.p'),'rb'))
    y_pred = clf.predict(X_test)
    output_df = pd.DataFrame(y_pred, columns=['on.target.pred'])
    output_df = pd.concat([output_df, df], axis=1)
    return output_df

if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    myco_prefixes = [parameter.data_dir + 'kmer_files/myco']

    print(create_augmented_df(myco_prefixes, ['ACAACC','TACGTCA'], molarity=2.5))