# SOAPswga (**S**elective **O**timal **A**mplifying **P**rimers for **s**elective **w**hole-**g**enome **a**mplification)

## Contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Workflow](#workflow)
* [Parameters](#parameters)

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

The k-mer files will be output into the directory

### Step 2: Candidate primer filtering 

In this step, we filter out candidate primers from the set of all motifs in the target genome based on having certain properties of each primer (as suggested by \cite{biosoft}). These are as follows:

\begin{enumerate}
    \item Melting temperature: the melting temperature must be within the standard temperature range of \verb!min_tm! (default $15^{\circ}$C) and \verb!max_tm! (default $45^{\circ}$C) as established in \cite{Leichty2014}. These temperatures are computed based on \cite{allawi1997}. 
    
    \item Self-dimer: candidate primers that could possibly form self-dimers are eliminated at this stage. This is estimated by the longest common subsequence between the candidate primer and its reverse complement being greater than \verb!default_max_self_dimer_bp! (default 4).
    
    \item Number of consecutive runs of a single base pair: primers with runs of 5 of five or longer are eliminated. 
    
    \item GC content: The GC content for each primer is maintained to be within \verb!min_GC! and \verb!max_GC! (default 0.375 and 0.625 respectively). Three or more G or C's are avoided in the last five and three base pairs of the 3'-end of the primer.
    
    \item Di-nucleotide repeats: di-nucleotide repeats of 5 or more are avoided.
    
    \item Binding frequency: binding frequency is computed as the number of exact matches of a primer in the genome, normalized by the total genome length. Primers that bind too sparsely to the target genome (lower than parameter \verb!min_fg_freq!) or too frequently to the off-target (higher than \verb!max_bg_freq!) are removed. 
    
    \item Binding evenness: the evenness of binding is calculated by finding the Gini index of the distances between each primer binding site on the target, and primers with Gini indices higher than \verb!max_gini! are removed. 
    
\end{enumerate}

While the melting temperature filter is more a requirement of the experimental procedure, most of these filters are done in an effort to reduce mis-priming (binding to an undesirable location, e.g. the off-target genome or another primer). The computing of the frequencies of these candidate primers in the genomes is critical for estimating the number of possible binding positions in the target and off-target genomes, the former of which we would like to maximize and the latter of which we would like to minimize. Binding evenness as measured by the Gini index of the gap distances is important for downstream analysis such as short-read assembly and copy number estimation.

Finally, primers are ranked by the ratio of the binding frequency in the target genome(s) to the binding frequency in the off-target genome(s) and those primers with the highest ratio are identified for downstream use (by default, this currently identifies the top 500 primers, and is modifiable via the \verb!max_primers! parameter). 

An h5py file is then created for storing primers and their respective locations for each genome (or chromosome). 

\subsubsection*{[Step 3] Amplification efficacy scoring: predicting individual potential strength to amplify} %DB - this section needs additional clarity JY - I'm not sure if I addressed your comments sufficiently. I really just intended to explain what this step does, just as the other steps do, and I was hoping to defer explanation of the random forest to the plasmid experiment section, which is how the random forest is trained. 
In this step, before we optimize over the combinatorial space of primer sets, we score each individual candidate primer from the previous step. To do this we use a random forest regressor model trained from prior experimentation (discussed in \nameref{plasmidexp}). In short, this non-linear regression model is trained on plasmid experiments conducted using sets of a single primer and a single plasmid genome. The goal of this regressor is to predict amplification efficacy from various properties of the primer, including computed thermodynamically-principled features estimating a primer's binding affinity for the target genome (see \nameref{sec:tp_features}). %DB - not clear what this means. The features are not being generated. %DB - this is confusing. The random forest regressor has not been discussed and it is not clear yet where it comes from. It can be made clear here that optimal therodynamic features have been identified experimentally using random forest models. 
After predicting an on-target amplification value for each of the candidate primers from the previous step, the primers are filtered out according to the minimum predicted on-target amplification threshold (\verb!min_amp_pred!) parameter, which by default is set to 5. This step significantly reduces search computation by weeding out low-amplification primers. 

\subsubsection*{[Step 4] Primer set search and evaluation}

Using the filtered list of primers from the previous step, SOAPswga searches for primer sets using a machine-learning guided scoring function and a breadth-first, greedy approach. At a high level, \verb!max_sets! number of top primer sets are built-in parallel, primer by primer, by adding primers which increases the evaluations score the most. More specifically, we run the following algorithm. \\ % DB-some additional clarity here too. For example, it may not be clear to the reader what the filtered list of is. 
\begin{algorithm}[H] % DB-I am not sure this will work for readers. A schematic could work but probably not computer language. Maybe in a supplement TODO
    \SetKwFunction{Evaluate}{Evaluate}
    \SetKwFunction{ChooseMaxSets}{ChooseMaxSets}
    \SetKwFunction{Compatible}{Compatible}
    \SetKwFunction{RandomInitialStart}{RandomInitialStart}
    \SetKwFunction{DropOut}{DropOut}
    \KwIn{primer\_list (list of candidate primers)}
    \additionalinput{max\_sets (maximum number of sets to explore at each step)}
    \additionalinput{primer\_efficacy\_scores (amplification efficacy predicted from random forest model)}
     \additionalinput{drop\_indices (defines which iterations are drop out layers)}
    \KwOut{top\_sets}
    
    top\_sets $\leftarrow$ \RandomInitialStart(primer\_list, primer\_efficacy\_scores, max\_sets);
    
    top\_scores $\leftarrow$ \Evaluate(top\_sets);
    
    curr\_sets = $\leftarrow$ []\;
    curr\_scores $\leftarrow$ []\;
    
    \For {$i = 1$ to $max\_iterations$}{
        \For {top\_set in top\_sets}{ 
            \For {primer in S}{
                \If{\Compatible(top\_set $\cup$ [primer])}{
                new\_set $\leftarrow$ [top\_set $\cup$ [primer]]\; 
                score $\leftarrow$ \Evaluate(new\_set)\;
                curr\_scores $\leftarrow$ curr\_scores $\cup$ [score]\;
                curr\_sets $\leftarrow$ curr\_sets $\cup$ new\_set;
                }
            }
        }
        % \If{max(curr\_scores) $<$ min(top\_scores)}{
        %     \Return top\_sets
        % }
        \If{i is in drop\_indices}{
            top\_sets, top\_scores = \DropOut(curr\_scores, curr\_sets, max\_sets)
        }
        \Else{
        top\_sets, top\_scores = \ChooseMaxSets(curr\_scores, curr\_sets, max\_sets)}
    }
    \Return top\_sets
    \caption{Primer Set Search}
\end{algorithm}
The function $\RandomInitialStart$ randomly chooses the first primer in each of the \verb!max_sets! number of primers where probabilities are the normalized amplification efficacy scores from the previous step. This allows for randomization of the search initialization and permits better exploration of the search space. 
The function $\Compatible$ checks that no two primers have subsequences longer than \verb!default_max_dimer_bp! that are complementary. This is done in an effort to reduce the risk of primers binding with each other rather than the target genome. %DB - it is not explained why heterodimers are bad. I do not know what "no risk" means.
The function $\DropOut$ allows each of the highest-scoring primer sets to drop one primer. This is particularly useful, for example, if a primer set has an evaluation score much higher than the primer set excluding the initially chosen primer. ``Dropping'' a primer also allows for the possibility of adding a primer that would otherwise be barred because of its risk in forming dimers with the dropped primer.
$\ChooseMaxSets$ simply chooses up to \verb!max_sets! number of best sets based on their score.

$\Evaluate$ evaluates the primer set using the following equation:
\begin{align*}
\text{Score} = 
  \beta_0 + 
  \beta_1\,\text{freq\_ratio} + 
  \beta_2\,\text{mean\_gap\_ratio} + 
  \beta_3\,\text{coverage\_ratio} +
  \beta_4\,\text{on\_gap\_gini} + 
  \beta_5\,\text{off\_gap\_gini}
 \end{align*} \label{eqn:first}
where Table~\ref{tab:regression} contains the coefficient values and the score is a prediction of the 1x coverage per base pair sequenced. This step can be re-run multiple times and includes the option of withholding primers too frequently used in previous primer sets until after the dropout layer.
