import utility
import parameter
import numpy as np

def compatible_set(dimer_mat, selected_primers, primer_to_index_dict):
    for i, primer in enumerate(selected_primers):
        for j in range(i+1, len(selected_primers)):
            if dimer_mat[primer_to_index_dict[primer]][primer_to_index_dict[selected_primers[j]]] == 1:
                print(primer)
                print(selected_primers[j])
                return False
    return True

def compatible(dimer_mat, selected_primers, primer, primer_to_index_dict):
    for selected_primer in selected_primers:
        if dimer_mat[primer_to_index_dict[primer]][primer_to_index_dict[selected_primer]] == 1:
            return False
    return True

def is_dimer(seq_1, seq_2, max_dimer_bp=parameter.max_dimer_bp):
    binding_len = utility.longest_common_substring(seq_1, utility.reverse_complement(seq_2))
    return(binding_len > max_dimer_bp)

def heterodimer_matrix(primer_list, max_dimer_bp=parameter.max_dimer_bp, het_type='boolean'):

    n = len(primer_list)
    het_matrix = np.zeros((n, n))

    # Pairwise comparison, note this can be precomputed
    for i in range(n):
        het_matrix[i, i] = is_dimer(primer_list[i], primer_list[i], max_dimer_bp=parameter.max_self_dimer_bp)
        for j in range(i+1, n):  # Also check internal hairpin (binds with itself)
            het_matrix[i, j] = is_dimer(primer_list[i], primer_list[j], max_dimer_bp=max_dimer_bp)
            het_matrix[j, i] = het_matrix[i, j]
    return het_matrix

def is_compatible_set(primer_set, max_dimer_bp=parameter.max_dimer_bp):
    het_mat = np.matrix(heterodimer_matrix(primer_set, max_dimer_bp=max_dimer_bp))
    if het_mat.sum() > 0:
        return False
    return True

if __name__ == "__main__":
    primer_set = ["ATCGACAAC","CGAATCCG","CGTTACGG","CTACGACG","GACGATCG","GATCGACTC","TCGACGAA"]

    dimer_mat = heterodimer_matrix(primer_set)
    print(dimer_mat)
    print(is_compatible_set(primer_set))
    print(is_dimer(primer_set[5], primer_set[6]))
