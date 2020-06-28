import os
import src.parameter
import melting

def run_jellyfish(genome_fname=None, output_prefix=None):
    """
    Runs jellyfish program using the output_prefix and transfroms the kmer count information txt files. Count k-mers from 6 to 12.

    Args:
        genome_fname: The fasta file used to count kmers.
        output_prefix: The output path prefix for the output files. Resulting output files will be suffixed by _kmer_all.txt for k from 6 to 12 inclusive.
    """
    for k in range(6, 13, 1):
        if not os.path.exists(output_prefix+'_'+str(k)+'mer_all.txt'):
            os.system("jellyfish count -m "+str(k) + " -s 1000000 -t " + str(src.parameter.cpus) + " " + genome_fname + " -o " + output_prefix+'_'+str(k)+'mer_all.jf')
            os.system("jellyfish dump -c " + output_prefix+'_'+str(k)+'mer_all.jf' + " > " + output_prefix+'_'+str(k)+'mer_all.txt')
        if os.path.exists(output_prefix+'_'+str(k)+'mer_all.jf'):
            os.system("rm " + output_prefix+'_'+str(k)+'mer_all.jf')

def get_kmer_to_count_dict(f_in_name):
    """
    Computes the counts of all kmers in the 5' to 3' direction.

    Args: The path to the fasta file which has the genome.

    Returns:
        primer_to_count_dict: Dictionary mapping kmer to count in the 5' to 3' direction.
    """
    primer_to_count_dict = {}

    with open(f_in_name, 'r') as f_in:
        for line in f_in:
            primer = line.split(" ")[0]
            count = int(line.split(" ")[1].rstrip())
            primer_to_count_dict[primer] = count

    return primer_to_count_dict

def get_primer_list_from_kmers(prefixes, kmer_lengths=None):
    """
    Gets all the kmers from the jellyfish output files.

    Args:
        prefixes: The prefix path that all the jellyfish output files share.

    Returns:
        primer_list: List of all the kmers that occur at least and satisfy the temperature conditions.
    """
    primer_list = []

    if not kmer_lengths:
        kmer_lengths = range(6,13,1)

    for prefix in prefixes:
        for k in kmer_lengths:
            with open(prefix+'_'+str(k)+'mer_all.txt', 'r') as f_in:
                for line in f_in:
                    curr_kmer = line.split(" ")[0]
                    tm = melting.temp(curr_kmer)
                    if tm < src.parameter.max_tm and tm > src.parameter.min_tm:
                        primer_list.append(curr_kmer)
    return primer_list
