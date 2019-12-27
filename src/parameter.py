import os
import multiprocessing
import json
import utility

src_dir=os.path.dirname(os.path.abspath(__file__))

def get_all_files(prefix_or_file):
    if type(prefix_or_file) is str:
        if 'fasta' == prefix_or_file.split('.')[-1] or 'fa' == prefix_or_file.split('.')[-1]:
            print(prefix_or_file.split('.')[-1])
            return [prefix_or_file]
        else:
            prefix_or_file = [prefix_or_file]
    total_files = []
    for individual_prefix in prefix_or_file:
        list_of_files = os.listdir(os.path.dirname(individual_prefix))  # list of files in the current directory
        for each_file in list_of_files:
            if each_file.startswith(individual_prefix.split('/')[-1]):
                if each_file.endswith('fasta') or each_file.endswith('fa'):
                    total_files.append(os.path.join(os.path.dirname(individual_prefix), each_file))
    # print(total_files)
    return total_files

def get_value_or_default(arg_value, data, key, default):
    # print(arg_value)
    # print(data)
    # print(key)
    # print(default)
    if arg_value is not None:
        return arg_value
    if not key in data:
        return default
    return data[key]

def write_args_to_json(args, out_fname=None):
    global min_fg_freq
    global max_bg_freq
    global min_tm
    global max_tm
    global max_gini
    global max_primer
    global min_amp_pred
    global cpus
    global max_dimer_bp
    global max_self_dimer_bp
    global mismatch_penalty
    global data_dir
    global selection_metric
    global iterations
    global max_sets
    global fg_circular
    global bg_circular
    global drop_iterations

    data = {}

    data_dir = data['data_dir'] = args.data_dir
    # data['json_file'] = get_value_or_default(args.json_file, data, 'json_file', os.path.join(data['data_dir'], 'params')) if out_fname is None else out_fname

    if args.json_file is not None:
        data['json_file'] = args.json_file
        if not os.path.exists(os.path.dirname(data['json_file'])):
            os.makedirs(os.path.dirname(data['json_file']))

        data_extra = read_args_from_json(args.json_file)

        for k, v in data_extra.items():
            if k not in data:
                data[k] = v

    min_fg_freq = data['min_fg_freq'] = args.min_fg_freq
    max_bg_freq = data['max_bg_freq'] = args.max_bg_freq
    min_tm = data['min_tm'] = args.min_tm
    max_tm = data['max_tm'] = args.min_tm
    max_gini = data['max_gini'] = args.max_gini
    max_primer = data['max_primer'] = args.max_primer
    min_amp_pred = data['min_amp_pred'] = args.min_amp_pred
    cpus = data['cpus'] = args.cpus
    max_dimer_bp = data['max_dimer_bp'] = args.max_dimer_bp
    max_self_dimer_bp = data['max_self_dimer_bp'] = args.max_self_dimer_bp
    mismatch_penalty = data['mismatch_penalty'] = args.mismatch_penalty

    if args.fasta_fore is not None:
        data['fg_genomes'] = get_all_files(args.fasta_fore)
    if args.fasta_back is not None:
        data['bg_genomes'] = get_all_files(args.fasta_back)

    if args.kmer_fore is not None:
        if len(data['fg_genomes']) > 1:
            data['fg_prefixes'] = [args.kmer_fore + os.path.splitext(os.path.basename(genome_fname))[0] for genome_fname
                                   in data['fg_genomes']]
        else:
            data['fg_prefixes'] = [args.kmer_fore]
    if args.kmer_back is not None:
        if len(data['bg_genomes']) > 1:
            data['bg_prefixes'] = [args.kmer_back + '_' + os.path.splitext(os.path.basename(genome_fname))[0] for
                                   genome_fname in data['bg_genomes']]
        else:
            data['bg_prefixes'] = [args.kmer_back]

    if 'fg_seq_lengths' not in data or len(data['fg_seq_lengths']) != len(data['fg_genomes']):
        data['fg_seq_lengths'] = utility.get_all_seq_lengths(fname_genomes=data['fg_genomes'], cpus=data['cpus'])
    if 'bg_seq_lengths' not in data or len(data['bg_seq_lengths']) != len(data['bg_genomes']):
        data['bg_seq_lengths'] = utility.get_all_seq_lengths(fname_genomes=data['bg_genomes'], cpus=data['cpus'])

    with open(data['json_file'], 'w+') as outfile:
        json.dump(data, outfile, indent=4)

    return data

def read_args_from_json(in_fname):
    with open(in_fname) as json_file:
        data = json.load(json_file)
    return data

if __name__ == "__main__":
    print(cpus)