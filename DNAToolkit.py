# DNA Toolkit file
import collections
import random
import re
import itertools
import math
import numpy

Nucleotides = ["A", "C", "G", "T"]
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# Creating a random DNA sequence for testing:
randDNAStr = ''.join([random.choice(Nucleotides)
                      for nuc in range(100)])


def validateSeq(dna_seq):
    """Check the sequence to make sure it is a DNA String"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(seq):
    """Counts the frequency of types of nucleotides"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    # "tmpFreqDict" is a Dictionary first thing in the ratio is the term : second thing is the value
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict


def transcription(seq):
    """ DNA to RNA transcription """
    return seq.replace("T", "U")


def reversecomplement(seq):
    return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]  # returns 5' to 3'


def GCcontent(seq):
    """Counts the frequency of types of G and C base pairs"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    GC = 0
    # "tmpFreqDict" is a Dictionary first thing in the ratio is the term : second thing is the value
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    GC += tmpFreqDict["G"]
    GC += tmpFreqDict["C"]
    return GC / len(seq)


def SubsectionGCcontent(seq, k):
    """Counts the frequency of types of G and C base pairs within subsection of "k" length"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    GC = 0
    subsections = []
    # "tmpFreqDict" is a Dictionary first thing in the ratio is the term : second thing is the value
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        subsections.append(GCcontent(subseq))
    return subsections


def HighestGCcontent(file):
    """Takes a fasta file returns the highest GC content and its index in the list provided by the function read_one_fasta IN PERCENT"""
    (idnum, ID, seqlist) = read_one_fasta(file)  # defining multiple variables from one function return
    GClist = []  # stores a list of GC content for each unique string found
    X = len(seqlist)  # length of the list for the following loop
    hold = 0
    holdindex = 0  # stores index of the highest found GC content
    for i in range(X):
        GClist.append(GCcontent(seqlist[i]))  # translates the list of sequences into a list of GC contents
        if GClist[i] > hold:
            hold = GCcontent(seqlist[i])
            holdindex = i

    # print(GClist) #, print(ID[holdindex])
    return ID[holdindex], (hold * 100)


def read_one_fasta(file):
    # function to read one sequence in FASTA format from a textfile.
    # input: the name (complete path) of the file that contains the sequence to read
    # output: strings containing the Identifier and the Sequence
    # this function will read ONE sequence only.  If more that one sequence is in the file, only first is read

    F = open(file)
    idnum = 0  # will keep track of (index) how many ID lines have been found
    ID = []
    seq = ""  # will hold Identifier (name) of sequence
    seqlist = []  # will hold Sequences

    for line in F:
        line = line.rstrip("\n")  # removes new line characters (returns) from line
        m = re.search(r"^>(\S+)", line)  # uses a Regular Expression to check if line contains ">"
        if m and idnum > 0:
            ID.append(line)  # remember Identifier
            seqlist.append(seq)
            idnum = idnum + 1  # remember ID#
            seq = ""
        elif m and idnum == 0:
            ID.append(line)  # remember ID
            seq = ""  # make sure seq is empty
            idnum = idnum + 1
        else:
            seq += line  # IF not ID line AND IF if in 1st seqeunce, add line to sequence
            # print("found sequence info")
    seqlist.append(seq)
    F.close()
    return idnum, ID, seqlist


def Fibonacci(n, k):
    """Takes n months and k number of offspring per breeding pair nad returns number of breeding adults"""
    babies = 1
    adults = 0
    newbabies = 0
    for i in range(n):
        newbabies = adults * k
        adults = adults + babies
        babies = newbabies
        newbabies = 0
    return (adults)


def LifespanFibonacci(n, l):
    """Takes n months and a lifespan of l, assumes each breeding pair has one viable offspring"""
    babies = 1
    adults = 0
    newbabies = 0
    reaper = [0]  # this lists how many adults that will die, each month is a different spot on the list
    for i in range(n):

        print("month " + str(i))
        print("Killing next number in list")
        print(reaper)
        adults = adults - reaper.pop(
            0)  # takes the current number of adults that are supposed to die this month and removes it from the list
        print("remaining list")
        print(reaper)
        newbabies = adults
        if i < 1:
            for x in range(l - 1):  # at the begining of the simulation, you need to add how many months will pass before the first adult dies of old age
                reaper.append(0)
        reaper.append(babies)  # notes how many deaths will occur after lifespan ends
        print("added babies to death list")
        print(reaper)
        adults = adults + babies
        babies = newbabies
        newbabies = 0

    return adults


def PointMutationCounter(seq1, seq2):
    """Takes two strings checks if theres a mutation at each position"""
    mutationcount = 0
    print("captured seq")
    print(str(len(seq1)))

    for y in range(len(seq1)):
        print("comparing at position " + str(y))
        print(seq1[y])
        print(seq2[y])
        if seq1[y] != seq2[y]:
            mutationcount = mutationcount + 1
    return (mutationcount)


def Permutations(n):
    print(math.factorial(n))
    perm = itertools.permutations(list(range(1, n + 1)))
    for i, j in enumerate(list(perm)):
        permutation = ''
        for item in j:
            permutation += str(item) + ' '
        print(permutation)





DNAstr1 = "ATGCACTCCGCATGCCTATTATATAATGCTGTTTGATACCGGTCTTCTGTAGGGGTTCGTTTAGTTCGTGCCGGTATCCGCTCATCCGGTCACCCGCAAGAAACCTCAGGAGTGTAACGTAGCCCCTGGACGGATGATCACGACAACACAGCCCATGGTCTGGGCACTTGAGATGCTAGGGTGACGCCGTAACTTATGGCCTTTCTACTCCCTGATGGAGATTACACATTCTTCTATGCTTGTCACTAATAGTATCACGTCCCAAGGGCTCTCAAGCAAATCGCTGTTGTCGTATCTGTTCGTGAGTGCGGCAAGGAGTGTCTGATGGGTCTACAACGGAGTTCGCTCAAAATTGAACGTGTGATGCAGGAGATTATTCAGCCACGCGCGTTACGGAGTTGATTATATGAACAGTCAGCTTTATAGCGCCATGAGACTCACGTCGAATAATCCACCCGCGGCAATCCCTAGCAAGACGCACTTTCGAGCCAATCACGTTCATGCCAGCTCGGTCCTTAACTTGGGAGCGGCACATATATGCCAGTGCTGGTGAAGGCTTAATAGCAACATTATACCAACACTGCGGGAAATCTTTGGTAAGGAGAATGTCCAGTGATTGGCATTGAGTGTGAGACTGGTGAATCTACGGTGCCAATGACACCTAGGCCGGGGGTTGATGTGCTTCGCCCTGAAACATATTCTGTTTGGCCGATTGCGCTTTTCTATCTAATAATCTAGATTGTAGTAGGGATCTTCGCCCTTCCCTTAACCGAGACGTTCTATTCCCAGGCCCAAAGCGAGGTGGGTCGCGCTGCTTGGACGGCCACAGGTTTCCCCATTGTTTTAAGTACGTTACAGATTTCTTCATGCATACCGCCTAGGTGCAGCACAGCACAGATTATGGGGCCCAAAGCACGGTCTTACGAAAACTCTCGCCGTGCTATTTACGTGTCTACAACAACATCTAAACGGATTGACAAGAA"
DNAstr2 = "CCGGCACAGGCCTTACAATTATACCTCACTGATTTATTCTTGGCTTCCGTAGAAGCTCGCTCAGTTCTAGCTAGTAACCGCTTGACCGGGTAGACTATAGAGACCGTGCGAGTGGCACTTAGGGTCTTGTCTGATAATCACCTTGTCAGAACCGGATAGCGAGTAACTCCATGGTCAAGAGTTATGCTCGCATTTATCGCCTGACGTTACAGCGCTGGGTAATCCCTATTGGCTTAGGCAGCTGACGGAGTGTAGCTAGTCACAGCCCGGCGCTTGCAATTACACCTTCCAATCCGGGATTCTAAATCCAACAGGTTTGGCGTGAACCTCGCATGGGGGAAAACTCCTCATAGGGAAGGACATATGCATAAGTTTTTGCTGACAAATGGGAGAGGAATTACATATACTGCAGCGTCCGCATGATTGCACCATGACACTCTAATAAATTATTGCGCGACACGTAATGCCCTACCTGGCGGACACGCTACACTATTCCTAAGGGGACACAGCCGTCATTCAACCCGGAGCGGCTGCAAATTCGCCGTAAAGTCCAACTCTGATTATCGATGTCAAACCTCTAGTGCACGAGGTCTTGGGGCACGTATAACTCTAATCTGGCGCGAATGCTGTCAGAACCCACGGACAGCTGAACATTAAGCATCCAGGCCGGAATTGGGGGTGCATCACCAGGAGGAAGTTTCTAAGGGACCGCGTTCTCGATACAATTAGGTTGGGCATCTCGGAGAACCGCCCTACGGACTGCCATCAGTGGGCATGATTTGCCCTCAGAGGATATGGTACGGCGCGCACCCACGTACGGCGTCATGAGCTACATCTATTGTTTTTAGGAGTATAGGGATTTTGGTTGATATCCAGCCCAAAGACGTTAGAGGCCGAGTTGAGAGGAACTACCCGCGGACCACTGTCAAGTGGGTCTCTTGTGGATGTCAATTACAAGTTAATGCTTACTCAATTCATTAGCC"

# if validateSeq(randDNAStr) != True:
#     print(randDNAStr)
#     print(countNucFrequency(randDNAStr))
#     print(transcription(randDNAStr))
#     print(randDNAStr)
#     print(reversecomplement(randDNAStr))
#     print(GCcontent(randDNAStr))
#     print(SubsectionGCcontent(randDNAStr, 10))
#     print(f" {read_one_fasta('Amel_DSCAM.FASTA')} ")
#     print(f" {LifespanFibonacci(94, 19)} ")
#     print(HighestGCcontent('input.txt'))
# print(str(PointMutationCounter(DNAstr1, DNAstr2)))
print(Permutations(3))
