# SOAPswga (**S**elective **O**timal **A**mplifying **P**rimers for **s**elective **w**hole-**g**enome **a**mplification)

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Workflow](#workflow)
* [All parameters](#parameters)

## Introduction 

SOAPswga is a command-line tool for choosing   sets of primers for selective whole-genome amplification. The overall pipeline for proposing primer sets takes in a set of genomes which we would like to amplify (target genome(s)) and would not like to amplify (off-target genome(s)), where the goal is to find sets of primers which will effectively and evenly amplify the target genome(s) and not the off-target genome(s). 

## Installation

We recommend that you install `soapswga` into a python [virtual environment](https://virtualenv.pypa.io/en/latest/).

For example, you can create the virtual environment `soapswga_venv` using

```bash
$ virtualenv soapswga_venv
```

Enter the virtual environment using 

```bash
$ source soapswga_venv/bin/activate
```

Install dependencies in the virtual environment using

```bash
$ pip install -r requirements.txt
```

To install `soapswga`

 ```bash
$ pip install soapswga
```

## Workflow

The general workflow has four main stages: 1) preprocessing of locations in the target and off-target genome of all motifs in the target genome, 2) filtering all motifs in the target genome based on individual primer properties and frequencies in the genomes, 3) scoring the remaining primers for amplification efficacy using a machine learning model, and 4) search and evaluate aggregations of primers as candidate primer sets.

### Step 1: K-mer preprocessing 

The primary step in the program identifies the k-mers of length 6 to 12 in the target genome, which serve as possible candidate primers for downstream steps. The counts of these k-mers in the target and off-target genome(s) are computed using \verb!jellyfish! (\cite{marccais2011fast}), a fast, parallel k-mer counter for DNA. This entire preprocessing step is parallelized, and we have provided pre-computed files for \textit{Mycobacterium tuberculosis} as well as human. If parameters are changed, this step will not need to be re-run. 

The k-mer files for the target and off-target gnomes will be output to a file with a path prefix determined by command line parameters `-k` or `--kmer-fore` and `-j` or `--kmer-back`, respectively. 

For example, if you have a ready json_file (see example/params.json) you can run 

```bash
$ step1 -j example.json
```
and for each prefix in `fg_prefixes` and `bg_prefixes`, a file with suffixes `_Xmer_all.txt` for X = 6 to 12 containing all k-mers of length X. The program will use the fasta file in `fg_genomes` or `bg_genomes` at the corresponding index. Thus, `fg_prefixes` should have the same length as `fg_genomes` and the same for `bg_prefixes` and `bg_genomes`. It they do not have the same length, the name of the fasta file will be used as the file prefix but at the path specified by `fg_prefixes` or `bg_prefixes` if it exists.

Alternatively, if there is no pre-existing json file, the k-mer files for the target and off-target gnomes will be output to a file with a path prefix determined by command line parameters `-k` or `--kmer-fore` and `-j` or `--kmer-back`, respectively. 

The genome files will be read from the `fg_genomes` and `bg_genomes` entries in the json file or by using all files with the path prefix determined by command line parameters `-f` or `--fasta-fore` and `g` or `--fasta-back`. 

Example:

```bash
$ step1 -j ../example/params.json
```
or if a json file does not exist/you want to overwrite parameters in the json:

```bash
$ step1 --kmer_fore ../example/kmer_files/myco --kmer_back ../example/kmer_files/human
 --fasta_fore ../example/genomes/MTBH37RV --fasta_back ../example/genomes/chr --json_file  ../example/params.json --data_dir ../example/
```
This would produce the given `.txt` files in `example/kmer_files/` and it would use all the fasta files with the path prefix `../example/genomes/MTBH37RV` and `../example/genomes/chr`.

#### Step 1 relevant parameters

| Short option | Long option | Default value | Description |
| ------------- | ------------- | ------------- | ------------- |
| -j | --json_file | None | path of json file, either existing or to be written |
| -k | --kmer_fore | None | path prefix for the kmer files of the target genomes | 
| -l | --kmer_back | None | path prefix for the kmer files of the off-target genomes |
| -x | --fasta_fore | None | path or path prefix to the fasta files of the on-target genomes | 
| -y | --fasta_back | None | path or path prefix to the fasta files of the off-target genomes |
| -z | --data_dir | soapswga/project/ | the project directory where metadata files will be stored |

### Step 2: Candidate primer filtering 

In this step, we filter out candidate primers from the set of all motifs in the target genome based on having certain properties of each primer (as suggested by `http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html`). Additionally, we filter by

- Binding frequency: binding frequency is computed as the number of exact matches of a primer in the genome, normalized by the total genome length. Primers that bind too sparsely to the target genome (lower than parameter `min_fg_freq`) or too frequently to the off-target (higher than `max_bg_freq`) are removed. 
- Binding evenness: the evenness of binding is calculated by finding the Gini index of the distances between each primer binding site on the target, and primers with Gini indices higher than `max_gini` are removed. 

Finally, primers are ranked by the ratio of the binding frequency in the target genome(s) to the binding frequency in the off-target genome(s) and those primers with the highest ratio are identified for downstream use (by default, this currently identifies the top 500 primers, and is modifiable via the `max_primers` parameter). 

A h5py files are then created for storing primers and their respective locations for each genome (or chromosome) in a same location as the kmer files.

Example:

```bash
$ step2 -j ../example/params.json
```

#### Step 2 relevant parameters
| Short option | Long option | Default value | Description |
| ------------- | ------------- | ------------- | ------------- |
| -j | --json_file | None | path of json file, either existing or to be written |
| -k | --kmer_fore | None | path prefix for the kmer files of the target genomes | 
| -l | --kmer_back | None | path prefix for the kmer files of the off-target genomes |
| -x | --fasta_fore | None | path or path prefix to the fasta files of the on-target genomes | 
| -y | --fasta_back | None | path or path prefix to the fasta files of the off-target genomes |
| -u | --min_fg_freq | 1e-5 | minimum normalized frequency of occurrences of the candidate primer in the foreground genomes |
| -v | --max_bg_freq | 5e-5 | maximum normalized frequency of occurrences of the candidate primer in the background genomes |
| -g | --max_gini | 0.6 | maximum Gini index of gap distances |
| -m | --min_tm | 15 | minimum predicted melting temperature |
| -n | --max_tm | 45 | maximum predicted melting temperature |
| -e | --max_self_dimer_bp | 4 | maximum number of self-complementary base pairs |
| -c | --cpus | all cpus | number of cpus to use for multi-processed tasks |
| -z | --data_dir | soapswga/project/ | the project directory where metadata files will be stored |

<!-- ```bash
$ step2 --min_fg_freq 1e-05 --max_fg_freq 5e-06 --max_gini 0.6 --max_primer 500 --min_amp_pred 5
``` -->

### Step 3: Amplification efficacy scoring

In this step, before we optimize over the combinatorial space of primer sets, we predict the on-target amplification value for each individual candidate primer from the previous step. Afterwards, we filter out according to the minimum predicted on-target amplification threshold (`min_amp_pred`) parameter, which by default is set to 5. This step significantly reduces search computation by weeding out low-amplification primers. 

Example:

```bash
$ step3 -j ../example/params.json
```
or if a json file does not exist/you want to overwrite parameters in the json:

<!-- ```bash
$ step3 --min_amp_pred 5
``` -->
#### Step 3 relevant parameters
| Short option | Long option | Default value | Description |
| ------------- | ------------- | ------------- | ------------- |
| -j | --json_file | None | path of json file, either existing or to be written |
| -a | --min_amp_pred | 5 | minimum amplification score from random forest regressor |
| -c | --cpus | all cpus | number of cpus to use for multi-processed tasks |
| -z | --data_dir | soapswga/project/ | the project directory where metadata files will be stored |

### Step 4: Primer set search and evaluation

Using the filtered list of primers from the previous step, SOAPswga searches for primer sets using a machine-learning guided scoring function and a breadth-first, greedy approach. At a high level, `max_sets` number of top primer sets are built-in parallel, primer by primer, by adding primers which increases the evaluations score the most. More specifically, we run the following algorithm. 


```bash
$ step4 -j ../example/params.json
```

<!-- ```bash
$ step4 --max-sets 5 --drop_iterations [4]
``` -->

#### Step 4 relevant parameters
| Short option | Long option | Default value | Description |
| ------------- | ------------- | ------------- | ------------- |
| -j | --json_file | None | path of json file, either existing or to be written |
| -r | --max_dimer_bp | 4 | maximum number of complementary base pairs |
| -s | --selection_method | 'deterministic' | selection method for choosing the next top primer sets |
| -d | --drop_iterations | [5] | the iterations which will have drop out |
| -i | --iterations | 10 | number of iterations to run the primer set search; the maximum length of the resulting primer sets will be the number of iterations minus the number of drop iterations |
| -c | --cpus | all cpus | number of cpus to use for multi-processed tasks |
| -z | --data_dir | soapswga/project/ | the project directory where metadata files will be stored |

## All parameters
| Short option | Long option | Default value | Description |
| ------------- | ------------- | ------------- | ------------- |
| -j | --json_file | None | path of json file, either existing or to be written |
| -k | --kmer_fore | None | path prefix for the kmer files of the target genomes | 
| -l | --kmer_back | None | path prefix for the kmer files of the off-target genomes |
| -x | --fasta_fore | None | path or path prefix to the fasta files of the on-target genomes | 
| -y | --fasta_back | None | path or path prefix to the fasta files of the off-target genomes |
| -c | --cpus | all cpus | number of cpus to use for multi-processed tasks |
| -z | --data_dir | soapswga/project/ | the project directory where metadata files will be stored |
| -u | --min_fg_freq | 1e-5 | minimum normalized frequency of occurrences of the candidate primer in the foreground genomes |
| -v | --max_bg_freq | 5e-5 | maximum normalized frequency of occurrences of the candidate primer in the background genomes |
| -g | --max_gini | 0.6 | maximum Gini index of gap distances |
| -m | --min_tm | 15 | minimum predicted melting temperature |
| -n | --max_tm | 45 | maximum predicted melting temperature |
| -e | --max_self_dimer_bp | 4 | maximum number of self-complementary base pairs |
| -a | --min_amp_pred | 5 | minimum amplification score from random forest regressor |
| -r | --max_dimer_bp | 4 | maximum number of complementary base pairs |
| -s | --selection_method | 'deterministic' | selection method for choosing the next top primer sets |
| -d | --drop_iterations | [5] | the iterations which will have drop out |
| -i | --iterations | 10 | number of iterations to run the primer set search; the maximum length of the resulting primer sets will be the number of iterations minus the number of drop iterations |
