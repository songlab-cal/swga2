import src.utility

#Global values for the thermodynamic nearest-neighbor model
delta_G_vals = {'GA/CA': 0.17,
                'GA/CC': 0.81,
                'GA/CG': -0.25,
                'GA/CT': -1.30,
                'GC/CA': 0.47,
                'GC/CC': 0.79,
                'GC/CG': -2.24,
                'GC/CT': 0.62,
                'GG/CA': -0.52,
                'GG/CC': -1.84,
                'GG/CG': -1.11,
                'GG/CT': 0.08,
                'GT/CA': -1.44,
                'GT/CC': 0.98,
                'GT/CG': -0.59,
                'GT/CT': 0.45,
                'CA/GA': 0.43,
                'CA/GC': 0.75,
                'CA/GG': 0.03,
                'CA/GT': -1.45,
                'CC/GA': 0.79,
                'CC/GC': 0.70,
                'CC/GG': -1.84,
                'CC/GT': 0.62,
                'CG/GA': 0.11,
                'CG/GC': -2.17,
                'CG/GG': -0.11,
                'CG/GT': -0.47,
                'CT/GA': -1.28,
                'CT/GC': 0.40,
                'CT/GG': -0.32,
                'CT/GT': -0.12,
                'AA/TA': 0.61,
                'AA/TC': 0.88,
                'AA/TG': 0.14,
                'AA/TT': -1.00,
                'AC/TA': 0.77,
                'AC/TC': 1.33,
                'AC/TG': -1.44,
                'AC/TT': 0.64,
                'AG/TA': 0.02,
                'AG/TC': -1.28,
                'AG/TG': -0.13,
                'AG/TT': 0.71,
                'AT/TA': -0.88,
                'AT/TC': 0.73,
                'AT/TG': 0.07,
                'AT/TT': 0.69,
                'TA/AA': 0.69,
                'TA/AC': 0.92,
                'TA/AG': 0.42,
                'TA/AT': -0.58,
                'TC/AA': 1.33,
                'TC/AC': 1.05,
                'TC/AG': -1.30,
                'TC/AT': 0.97,
                'TG/AA': 0.74,
                'TG/AC': -1.45,
                'TG/AG': 0.44,
                'TG/AT': 0.43,
                'TT/AA': -1.00,
                'TT/AC': 0.75,
                'TT/AG': 0.34,
                'TT/AT': 0.68
                }

nn_init_corrections = {'G/C': 0.98, 'A/T': 1.03, 'C/G': 0.98, 'T/A': 1.03}

def compute_free_energy_per_nn(doublet_1, doublet_2, penalty):
    """Computes the delta G values for two doublets.

    Args:
        doublet_1 (str): The 5' to 3' doublet.
        doublet_2 (str): The 3' to 5' doublet.
        penalty: The penalty for mismatches that are not specified (terminal or more than one mismatch).

    Returns:
        delta_G: The delta G value computed aligning doublet_1 and doublet_2.
    """
    nn = doublet_1 + '/' + doublet_2

    if nn in delta_G_vals:
        return delta_G_vals[nn]
    else:
        return penalty

# ASSUMES THAT THEY ARE OF THE FORM IN SANTA LUCIA PAPER.
# x is assumed 5' to 3' and y is assumed 3' to 5'
def compute_free_energy_for_two_strings(x, y, penalty=4):
    """Computes the delta G values for between two DNA sequences using the thermodynamic nearest-neighbor model.

    Args:
        x (str): The 5' to 3' sequence
        y (str): The 3' to 5' sequence
        penalty (double): The penalty for mismatches that are not specified (terminal or more than one mismatch).
            Default value is 4.

    Returns:
        delta_G: The delta G value computed aligning x and y.
    """
    x = x.upper()
    y = y.upper()
    delta_G = 0

    if x[0] == src.utility.complement(y[0]):
        delta_G = nn_init_corrections[x[0]+'/'+y[0]]

    if x[-1] == src.utility.complement(y[-1]):
        delta_G += nn_init_corrections[x[-1]+'/'+y[-1]]

    for i in range(1, len(x)):
        delta_G += compute_free_energy_per_nn(x[i - 1:i + 1], y[i - 1:i + 1], penalty)
        if delta_G > penalty * 10:
            break
    if src.utility.complement(x) == src.utility.reverse(y):
        delta_G += 0.43
    return delta_G

if __name__ == "__main__":
    primer_x = 'CGTTGA'
    primer_y = 'GCAACT'
    print(compute_free_energy_for_two_strings(primer_x, primer_y))