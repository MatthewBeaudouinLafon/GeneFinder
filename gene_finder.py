# -*- coding: utf-8 -*-
"""
YOUR HEAD HERE

@author: Matthew Beaudouin-Lafon

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # I cooould make a bunch of if statements, but that's no fun. Let's golf this shit.

    ATGC = 'ATGC'
    TACG = 'TACG'   # Appropriate map
    return TACG[ ATGC.index(nucleotide) ] # Find the index in ATGC, then return the element at that index in the map


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("C")
    'G'
    """
    dnaSize = len(dna) # Handy to avoid recomputation

    # Traverse the dna backwards, get the complement before appending
    reverse_compl = ''
    for i in range(dnaSize):
        reverse_compl = reverse_compl + get_complement(dna[dnaSize-i-1])

    return reverse_compl

    
def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGCGATCATAGGTTAG")
    'ATGCGATCA'
    """

    iter = 3 # Start after the start codon

    while iter < len(dna):

        daCodon = dna[iter:iter+3]
        if daCodon == 'TAG' or daCodon == 'TAA' or daCodon == 'TGA':
            return dna[:iter]

        iter += 3 # Traverse looking at codons rather than individual nucleotides

    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TACATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    demORFs = []
   

    # Find the start codon in that frame
    Iter = 0
    codon = dna[:Iter+3]

    # Find the rest of the sequence
    while Iter < len(dna):
        if dna[Iter:Iter+3] == 'ATG':
            demORFs.append(rest_of_ORF(dna[Iter:])) # Add next ORF
            Iter = Iter + len(demORFs[-1]) + 3  # Move by how long the ORF was, add three to skip the stop codon
        else:
            Iter += 3

    return demORFs



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    demORFs = []
    for i in range(3):
        demORFs.extend(find_all_ORFs_oneframe(dna[i:]))

    return demORFs 
    


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    allORFs = []
    allORFs.extend(find_all_ORFs(dna))  # Find ORFs in the dna
    allORFs.extend(find_all_ORFs(get_reverse_complement(dna)))

    return allORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string. Returns empty string if there are no ORFs
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("GATATGCGAATGTAGCATCAAA") # Extra codon in front
    'ATGCTACATTCGCAT'
    >>> longest_ORF("GATGCGAATGTAGCATCAAA") # Just one nucleotide in front
    'ATGCTACATTCGCAT'
    """
    allORFs = find_all_ORFs_both_strands(dna)
    largestORF = ''
    for anORF in allORFs:
        if len(anORF) > len(largestORF):
            largestORF = anORF

    return largestORF

    


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    uberLongestORF = 0
    for n in range(num_trials):    
        shuffledDNA = shuffle_string(dna)                   # Shuffle the DNA
        longestORFinShuffle = len(longest_ORF(shuffledDNA))      # Find the longest ORF in there
        if uberLongestORF < longestORFinShuffle:  # Make it the new uberLongestORF if it is longer 
            uberLongestORF = longestORFinShuffle

    return uberLongestORF


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGCCCGCTT")
        'MPA'
    """
    acids = ''
    Iter = 0
    while Iter+2 < len(dna):                            # Sketchy decision here
        acids = acids + aa_table[dna[Iter:Iter+3]]

        Iter += 3 # Move to the next Codon

    return acids


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    allORFs = find_all_ORFs_both_strands(dna)
    demAminos = [coding_strand_to_AA(ORF) for ORF in allORFs if len(ORF) > threshold] # Oh yeah comprehend that list

    return  demAminos


if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    print "Let's get started"
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print "Imported DNA, let's get computing"
    print gene_finder(dna)
