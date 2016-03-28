##############################################################################
###############              HELPER FUNCTIONS                #################
##############################################################################

import requests
import numpy as np
import pandas as pd
import re
from itertools import groupby
import pickle
import itertools
from itertools import product

def find_regions(annot,data):
    """This finds all the genes, non-coding regions,
    reverse regions in the genome and returns them as a list
    annot is the annotated genome with c's as forward genes,
    n's as noncoding regions and r's as genes on the reverse complement strand.
    data is the genome comprised of a's,t's,g's,c's.
    Both annot and data are strings"""
    stateGroups = groupby(enumerate(annot), lambda value: value[1])
    genome_regions = [(state, [x for x, _ in iterator]) for state, iterator in stateGroups]
    regions = [(region[0],min(region[1]),max(region[1])+1,
                data[min(region[1]):max(region[1])+1]) for region in genome_regions]
    return regions


def load_file(url,genome=True):
    """input: url for genome,
    output: pulled data without header and as a single txt string"""
    data = ''.join(requests.get(url).text.split("\n")[1::])
    return data

def load_file2(url,genome=True):
    """input: url for genome,
    output: pulled data without header and as a single txt string"""
    data = requests.get(url).text.split("\n")
    title = data[0][1::]
    genome = ''.join(data[1::])
    return [title,genome]


## load all Data
def all_data(staph_url,staph_annot):
    d = {}
    for i in range(len(staph_url)):
        title, genome = load_file2(staph_url[i])
        _, annot = load_file2(staph_annot[i])
        d[i] = [load_file(staph_url[i]),load_file(staph_annot[i]),title]
        print(title)
    return d

def find_regions_df(annot,data):
    region = find_regions(annot,data)
    regions = pd.DataFrame(region)
    regions.columns = [["state","start","end","seq"]]
#     print("Built table of sequences for each state in genome...")
    return regions


def segment_genome(data_string,split_size):
    """input: genome data string (data_string) and then size we should split the data (split_size)
    output: array split by specified split_size"""
    chunks = len(data_string)
    data = [data_string[i:i+split_size] for i in range(0,chunks,split_size)]
    return data

def gene_label(annotation, min_gene=10):
    """Input: annotation
    ouput: the label array for data (true if C or R else false for > 50% of line)"""
    size = len(annotation[0])
    print(size)
    label = [True if len(line.strip("C")) < size/2 or len(line.strip("R")) < size/2 else False
             for line in annotation]
    return label

def gene_label_dir(annotation,direction="C"):
    """Input: annotation and direction is C (forward) unless specified as R (reverse)
    ouput: the label array for data"""
    size = len(annotation[0])
    label = [True if len(line.strip(direction)) < size/2 else False
             for line in annotation]
    return label

def nucleotide_frequency(seq):
    '''Count the occurrences of characters in "seq".'''
    counts = {'A':0,'C':0,'G':0,'T':0}
    for c in seq:
        counts[c] +=1
    total = sum(list(counts.values()))
    for key, value in counts.items():
        counts[key] = value / total
    return counts

def ORF_finds(seq):
    """input: A sequence string (the genome) is taken and all the ATG's in the sequence are found
    output: a list of start indices
    *** on this reverse strand... these are from the end of the sequence!!"""
    starts = [m.start() for m in re.finditer('ATG', seq)]
    return starts

def sequence_list2(genome,ends,start,direction):
    """input: """
    frame = start % 3
    sequences = [(start, start+end+3, 'taa', frame, direction) for end in ends[0]]
    sequences += [(start, start+end+3, 'tag', frame, direction) for end in ends[1]]
    sequences += [(start, start+end+3, 'tga', frame, direction) for end in ends[2]]
    return sequences

def geneList(genome,direction):
    """input: genome (character string), which direction we are looking at (the forward or reverse strand)
    output: a list of sequences on that particular strand of interest (forward or reverse complement)"""
    starts = ORF_finds(genome)
    sequences = []
    for start in starts:
        ## m.start()+1 so that when we sample the sequence we go from
        ends = [[m.start() for m in re.finditer(x, genome[start:start+2000]) if m.start() % 3 == 0] for x in ["TAA", "TAG", "TGA"]]
        sequences.append(sequence_list2(genome,ends,start,direction))
    return sequences

def fixed(possible_reverse_sequences, genome_length):
    fixed_reversed_sequences = []
    for possible_reverse in possible_reverse_sequences:
        f = [(genome_length - end , genome_length - start,end_codon,frame,strand) for start,end,end_codon,frame,strand in possible_reverse]
        fixed_reversed_sequences.append(f)
    return fixed_reversed_sequences

def potentialGenes(genome):
    """input: genome which is a string of a's, t's, g's, c's
    ouput: the number of possible genes that start with ATG and
    end with taa,tag, or tga on the forward and reverse complement DNA strand"""
    forwards = geneList(genome,'F')
    reversed_genome = reverse_complement(str(genome))
    reverse = geneList(reversed_genome,'R')
    reverses = fixed(reverse, len(genome))
    combos = forwards + reverses
    return combos

def reverse_complement(genome):
    """input: genome, output: reverse complement"""
    switch = {"A":"T","G":"C","T":"A","C":"G"}
    rc = ''.join([switch[letter] for letter in genome])[::-1]
    return rc

def tri_split(seq):
    return ''.join([char+' ' if (i+1) % 3 == 0 else char for i,char in enumerate(seq)]).split()

def codonToAA(codon_list):
    aa_seq = ''.join([codon_to_aa[codon] for codon in codon_list])
    return aa_seq

def aa_frequency(seq):
    """Count the occurrences of amino acids in "seq"."""
    zeros = [0]*len(aa)
    counts = dict(zip(aa, zeros))
    for c in seq:
        counts[c] +=1

    total = sum(list(counts.values()))
    for key, value in counts.items():
        counts[key] = value / total
    return counts

def no_prior_stops(aa_seq):
    prior = aa_seq[0:len(aa_seq)-1]
    if 'x' in prior:
        return False
    else:
        return True

def find_all_occurrences(pattern, string):
    """finds all occurrences of pattern in string,
    returning a list of positions for each occurrence"""
    occurence_indices = [m.start() for m in re.finditer(pattern, string)]
    return occurence_indices

def stop_indices(seq):
    patterns = ["TAA","TAG","TGA"]
    stop_locations = [find_all_occurrences(pattern,seq) for pattern in patterns]
    stop_locations = list(itertools.chain(*stop_locations))
    return stop_locations

def start_indices(seq):
    pattern = "ATG"
    start_locations = find_all_occurrences(pattern,seq)
    return start_locations

def start_stops(genome):
    """input: output:"""
    stop_locations = sorted(stop_indices(genome))
    start_locations = start_indices(genome)
    return [start_locations, stop_locations]

def split_into_frames(starts,stops):
    """input: starts and stops are a list of indices
    where starts: (ATGs) or stops: (taa,tag,tga) where found in the genome.
    output: a list of tuples in the form of [...(start, stop, frame),...]"""
    res = []
    for i in range(3):
        s = np.array([-1] + [start for start in starts if (start - i) % 3 == 0])
        e = [end for end in stops if (end - i) % 3 == 0]
        indices = np.searchsorted(e,s)
        res += [(s[index],e[ind]+3,i) for index,ind in enumerate(indices[:-1])]
    return res

def process_genome(genome):
    """input:
    output: """
    starts, stops = start_stops(genome)
    pairs = [pair for pair in split_into_frames(starts,stops) if pair[0] !=-1]
    pairs = [pair for pair in pairs if pair[1] - pair[0] > 500]
    return pairs


def genome_end_genes(genome):
    """input:
    output: """
    starts, stops = closest_inframe_start(genome)
    pairs = [pair for pair in closest_start_end(starts,stops) if pair[0] !=-1]
    return pairs
