import os
import json
import src.utility

src_dir=os.path.dirname(os.path.abspath(__file__))

def get_all_files(prefix_or_file):
    """
    If passed a directory path, this returns all the fasta files in a directory. Otherwise, if given a file it
    passes back the file if it is a fasta file.

    Args:
        prefix_or_file: Directory path or path to a fasta file.

    Returns:
        total_files: List of fasta file(s).
    """
    if type(prefix_or_file) is str:
        if 'fasta' == prefix_or_file.split('.')[-1] or 'fa' == prefix_or_file.split('.')[-1]:
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
    return total_files

def get_value_or_default(arg_value, data, key):
    """
    Helper function for getting the input values. If a command line argument value was given, just return arg_value.
    Else we grab the value from the dictionary parsed from the json. This is to ensure command line inputs overwrite
    json inputs.

    Args:
        arg_value: Argument value to be returned
        data: dictionary of parameter to parameter values read from the input json file.
        key: The name of the parameter.

    Returns:
        val: Returns the command line input value or the json input value if command line one wasn't entered.
    """
    if arg_value is not None:
        return arg_value
    if key in data:
        return data[key]
    else:
        print("Please input values for " + str(key) + ".")

def get_params(args):
    """
    Writes the arguments of a pipeline instance to a json file for future use.

    Args:
        args: All the command line parameter inputs.

    Returns:
        data: Returns all the parameter inputs in dictionary form.
    """
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
    global verbose
    global top_set_count
    global retries
    global src_dir


    data = {}

    # data['json_file'] = get_value_or_default(args.json_file, data, 'json_file', os.path.join(data['data_dir'], 'params')) if out_fname is None else out_fname

    if args.json_file is not None:
        data['json_file'] = args.json_file
        if not os.path.exists(os.path.dirname(data['json_file'])) or not os.path.isfile(data['json_file']):
            if not os.path.exists(os.path.dirname(data['json_file'])):
                os.makedirs(os.path.dirname(data['json_file']))
            # if not os.path.isfile(data['json_file']):
            #     open(data['json_file'], 'w+').close()
        else:
            data_extra = read_args_from_json(args.json_file)

            for k, v in data_extra.items():
                if k not in data:
                    data[k] = v

    data_dir = data['data_dir'] = get_value_or_default(args.data_dir, data, 'data_dir')
    src_dir = data['src_dir'] = get_value_or_default(args.src_dir, data, 'src_dir')
    min_fg_freq = data['min_fg_freq'] = get_value_or_default(args.min_fg_freq, data, 'min_fg_freq')
    max_bg_freq = data['max_bg_freq'] = get_value_or_default(args.max_bg_freq, data, 'max_bg_freq')

    min_tm = data['min_tm'] = get_value_or_default(args.min_tm, data, 'min_tm')
    max_tm = data['max_tm'] = get_value_or_default(args.max_tm, data, 'max_tm')



    max_gini = data['max_gini'] = get_value_or_default(args.max_gini, data, 'max_gini')
    max_primer = data['max_primer'] = get_value_or_default(args.max_primer, data, 'max_primer')
    min_amp_pred = data['min_amp_pred'] = get_value_or_default(args.min_amp_pred, data, 'min_amp_pred')
    cpus = data['cpus'] = get_value_or_default(args.cpus, data, 'cpus')
    max_dimer_bp = data['max_dimer_bp'] = get_value_or_default(args.max_dimer_bp, data, 'max_dimer_bp')
    max_self_dimer_bp = data['max_self_dimer_bp'] = get_value_or_default(args.max_self_dimer_bp, data, 'max_self_dimer_bp')
    verbose = data['verbose'] = get_value_or_default(args.verbose, data, 'verbose')

    drop_iterations = data['drop_iterations'] = get_value_or_default(args.drop_iterations, data, 'drop_iterations')
    iterations = data['iterations'] = get_value_or_default(args.iterations, data, 'iterations')
    top_set_count = data['top_set_count'] = get_value_or_default(args.top_set_count, data, 'top_set_count')
    retries = data['retries'] = get_value_or_default(args.retries, data, 'retries')
    max_sets = data['max_sets'] = get_value_or_default(args.max_sets, data, 'max_sets')
    selection_metric = data['selection_metric'] = get_value_or_default(args.selection_metric, data, 'selection_metric')
    fg_circular = data['fg_circular'] = get_value_or_default(args.fg_circular, data, 'fg_circular')
    bg_circular = data['bg_circular'] = get_value_or_default(args.bg_circular, data, 'bg_circular')

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
        data['fg_seq_lengths'] = src.utility.get_all_seq_lengths(fname_genomes=data['fg_genomes'], cpus=data['cpus'])
    if 'bg_seq_lengths' not in data or len(data['bg_seq_lengths']) != len(data['bg_genomes']):
        data['bg_seq_lengths'] = src.utility.get_all_seq_lengths(fname_genomes=data['bg_genomes'], cpus=data['cpus'])

    return data

def read_args_from_json(in_fname):
    """
    Returns the parameters inputs in json form as a dictionary.

    Args:
        in_fname: Path to the json file.

    Returns:
        data: Parameter inputs in dictionary form.

    """
    with open(in_fname, 'r') as json_file:
        data = json.load(json_file)
    return data

if __name__ == "__main__":
    print(get_all_files("/Users/janeyu/Desktop/primer_data_dir/genomes/"))