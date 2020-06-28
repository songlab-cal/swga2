import src.utility
import src.parameter
import numpy as np


def compatible_set(dimer_mat, selected_primers, primer_to_index_dict):
    """
    Checks if the primer set selected_primers has a risk of forming heterodimers.

    Args:
        dimer_mat: A 2D binary, symmetric array where a 1 in row i and column j indicates primer i and primer j may form a heterodimer.
        selected_primers: The list of primers to be evaluated, written in the 5' to 3' direction.
        primer_to_index_dict: A dictionary from each primer to its corresponding index in dimer_mat.

    Returns:
        compatible_bool: True if no two primers in selected_primers are likely to form a heterodimer.
    """
    for i, primer in enumerate(selected_primers):
        for j in range(i + 1, len(selected_primers)):
            if dimer_mat[primer_to_index_dict[primer]][primer_to_index_dict[selected_primers[j]]] == 1:
                print(primer)
                print(selected_primers[j])
                return False
    return True


def compatible(dimer_mat, selected_primers, primer, primer_to_index_dict):
    """
    Checks if adding primer to primer set selected_primers has a risk of forming heterodimers.

    Args:
        dimer_mat: A 2D binary, symmetric array where a 1 in row i and column j indicates primer i and primer j may form a heterodimer.
        selected_primers: The already constructed primer set, written in the 5' to 3' direction.
        primer: Candidate primer being considered for addition to the primer set selected_primers, written in the 5' to 3' direction.
        primer_to_index_dict: A dictionary from each primer to its corresponding index in dimer_mat.

    Returns:
        compatible_bool: True if no two primers in selected_primers are likely to form a heterodimer.
    """
    for selected_primer in selected_primers:
        if dimer_mat[primer_to_index_dict[primer]][primer_to_index_dict[selected_primer]] == 1:
            return False
    return True


def is_dimer(seq_1, seq_2, max_dimer_bp=3):
    """
    Checks if two primers may form a heterodimer, according to if the length of the longest common substring is greater than max_dimer_bp. Adjust max_dimer_bp using options -t or --max_dimer_bp or the optional function argument.

    Args:
        seq_1: One of the primers, written in the 5' to 3' direction.
        seq_2: The other primer, written in the 5' to 3' direction.
        max_dimer: The maximum length of the longest common substring permitted. By default is set to 3.

    Returns:
        heterodimer bool: True if the primers may form a heterodimer.
    """
    binding_len = src.utility.longest_common_substring(seq_1, src.utility.reverse_complement(seq_2))
    return (binding_len > max_dimer_bp)


def heterodimer_matrix(primer_list, max_dimer_bp=3):
    """Computes a 2D binary, symmetric array where a 1 in row i and column j indicates primer i and primer j may form a heterodimer. Adjust max_dimer_bp using options -t or --max_dimer_bp or the optional function argument.

    Args:
        primer_list: The list of primers to evaluate.
        max_dimer_bp: The maximum length of the longest common substring permitted. By default is set to 3.

    Returns:
        het_matrix: 2D binary, symmetric array where a 1 in row i and column j indicates primer i and primer j may form a heterodimer and primer i corresponds to the ith element in primer_list.

    """
    n = len(primer_list)
    het_matrix = np.zeros((n, n))

    # Pairwise comparison, note this can be precomputed
    for i in range(n):
        het_matrix[i, i] = is_dimer(primer_list[i], primer_list[i], max_dimer_bp=src.parameter.max_self_dimer_bp)
        for j in range(i + 1, n):  # Also check internal hairpin (binds with itself)
            het_matrix[i, j] = is_dimer(primer_list[i], primer_list[j], max_dimer_bp=max_dimer_bp)
            het_matrix[j, i] = het_matrix[i, j]
    return het_matrix


def is_compatible_set(primer_set, max_dimer_bp=3):
    """
    Checks if primer_set has a risk of forming heterodimers. Adjust max_dimer_bp using options -t or --max_dimer_bp or the optional function argument.

    Args:
        primer_set: The list of primers to be evaluated, written in the 5' to 3' direction.
        max_dimer_bp: The maximum length of the longest common substring permitted. By default is set to 3.

    Returns:
        compatible_bool: True if no two primers in primer_set are likely to form a heterodimer.
    """
    het_mat = np.matrix(heterodimer_matrix(primer_set, max_dimer_bp=max_dimer_bp))
    if het_mat.sum() > 0:
        return False
    return True


if __name__ == "_main_":
    primer_set = ["ATCGACAAC", "CGAATCCG", "CGTTACGG", "CTACGACG", "GACGATCG", "GATCGACTC", "TCGACGAA"]

    dimer_mat = heterodimer_matrix(primer_set)
    print(dimer_mat)
    print(is_compatible_set(primer_set))
    print(is_dimer(primer_set[5], primer_set[6]))