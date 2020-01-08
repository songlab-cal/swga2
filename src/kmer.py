import os
import src.parameter
import melting

def run_jellyfish(genome_fname=None, output_prefix=None):
    for k in range(6, 13, 1):
        if not os.path.exists(output_prefix+'_'+str(k)+'mer_all.txt'):
            os.system("jellyfish count -m "+str(k) + " -s 1000000 -t " + str(src.parameter.cpus) + " " + genome_fname + " -o " + output_prefix+'_'+str(k)+'mer_all.jf')
            os.system("jellyfish dump -c " + output_prefix+'_'+str(k)+'mer_all.jf' + " > " + output_prefix+'_'+str(k)+'mer_all.txt')
        if os.path.exists(output_prefix+'_'+str(k)+'mer_all.jf'):
            os.system("rm " + output_prefix+'_'+str(k)+'mer_all.jf')

def get_kmer_to_count_dict(f_in_name):
    primer_to_count_dict = {}

    with open(f_in_name, 'r') as f_in:
        for line in f_in:
            primer = line.split(" ")[0]
            count = int(line.split(" ")[1].rstrip())
            primer_to_count_dict[primer] = count

    return primer_to_count_dict

def get_primer_list_from_kmers(prefixes, kmer_lengths=None):
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
