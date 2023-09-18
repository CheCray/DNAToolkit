# DNA Toolkit file
import collections
import random
import re
import itertools
import math
from collections import Counter
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

RawNucleotides = ["A", "C", "G", "T", "N"] # allows for raw read anlysis with bad basecalls
Nucleotides = ["A", "C", "G", "T"]
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N':'N'} # N allows for raw read analysis
# Creating a random DNA sequence for testing:
randDNAStr = ''.join([random.choice(Nucleotides)
                      for nuc in range(100)])

DNACodonTable = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W'}

RNACodonTable = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L',  # key: Value
                 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
                 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 'UCU': 'S', 'CCU': 'P',
                 'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 'CCG': 'P',
                 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
                 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': '_', 'CAA': 'Q',
                 'AAA': 'K', 'GAA': 'E', 'UAG': '_', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                 'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R',
                 'AGC': 'S', 'GGC': 'G', 'UGA': '_', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

ProteinCodonTable = {'F': ['UUU', 'UUC'], 'L': ['CUU', 'CUC', 'UUA', 'CUA', 'UUG', 'CUG'], 'I': ['AUU', 'AUC', 'AUA'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                     'M': ['AUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'],
                     'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'Y': ['UAU', 'UAC'], 'H': ['CAU', 'CAC'], 'N': ['AAU', 'AAC'], 'D': ['GAU', 'GAC'], '_': ['UAA', 'UAG', 'UGA'],
                     'Q': ['CAA', 'CAG'], 'K': ['AAA', 'AAG'], 'E': ['GAA', 'GAG'], 'C': ['UGU', 'UGC'], 'R': ['CGU', 'CGC', 'CGA', 'AGA', 'CGG', 'AGG'],
                     'G': ['GGU', 'GGC', 'GGA', 'GGG'], 'W': ['UGG']}

Monoisotopic_Mass = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694,
                     'E': 129.04259, 'F': 147.06841, 'G': 57.02146,
                     'H': 137.05891, 'I': 113.08406, 'K': 128.09496,
                     'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
                     'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
                     'S': 87.03203, 'T': 101.04768, 'V': 99.06841,
                     'W': 186.07931, 'Y': 163.06333}


def mRNAinference(AAseq): # WIP
    permutationcount = 3
    for AA in AAseq:
       permutationcount = permutationcount * len(ProteinCodonTable[AA])
    return permutationcount % 1000000

def reversedictionarycreator(seq): # Thanks Kevin
    """ Takes a dictionary and swaps keys with values """
    proteinToCodon = {}
    for key in seq:
        if seq[key] in proteinToCodon:
            proteinToCodon[seq[key]].append(key)
        else:
            proteinToCodon[seq[key]] = [key]
    print(proteinToCodon)

# reversedictionarycreator(RNACodonTable)

def motifFinder(seq, motif):
    """accepts a sequence string and a motif string and finds all the indexes that the motif is located in a list"""
    motifposition = []
    for position in range(len(seq)):
        if seq[position: position + len(motif)] == motif:
            motifposition.append(position+1)
    return motifposition

def listopenreadingframes(file):
    _, _, seq = read_one_fasta(file)
    RNAseq = transcription(seq[0])
    rcRNAseq = transcription(reversecomplement(seq[0]))
    orflist = []
    openreadingframe = ''
    for pos in range(len(RNAseq)):
        if RNAseq[pos:3 + pos] in RNACodonTable and RNACodonTable[RNAseq[pos:3 + pos]] == 'M':
            if translation(RNAseq[pos:len(RNAseq)]) == False:
                openreadingframe = ''
            else:
                openreadingframe += translation(RNAseq[pos:len(RNAseq)])
                orflist.append(openreadingframe)
                print(openreadingframe)
                openreadingframe = ''
    for pos in range(len(rcRNAseq)):
        if rcRNAseq[pos:3 + pos] in RNACodonTable and RNACodonTable[rcRNAseq[pos:3 + pos]] == 'M':
            if translation(rcRNAseq[pos:len(rcRNAseq)]) == False:
                openreadingframe = ''
            elif translation(rcRNAseq[pos:len(rcRNAseq)]) in orflist:
                openreadingframe = ''
            else:
                openreadingframe += translation(rcRNAseq[pos:len(rcRNAseq)])
                orflist.append(openreadingframe)
                print(openreadingframe)
                openreadingframe = ''
    print(f" Found {len(orflist)} open reading frames. Copied to 'Output.txt'")
    file = open('output.txt', 'w')
    for i in range(len(orflist)):
        file.write(orflist[i] + '\n')
    file.close()
    return orflist

def translation(seq) -> str:
    """Translates DNA or RNA seq into a protein stopping at the first '_' codon."""
    transcriptedseq = "" # Stores provided sequence in RNA form
    position = 0 # Stores position in the coding strand of DNA
    CodonSeq = "" # Stores translated sequence
    if validateDNASeq(seq) == False and validateDNASeq(reversetranscription(seq)) == False:  # Check the sequence to make sure it is a valid DNA/RNA String
        return print("Invalid Sequence Provided")
    else:
        transcriptedseq = transcription(seq)
        while position < len(transcriptedseq):
            block = transcriptedseq[position:3 + position]
            # print(block)
            if len(block) < 3 or RNACodonTable[block] == "_":
                #print('end block found')
                return CodonSeq
            else:
                CodonSeq += RNACodonTable[block]
                # print(gencode[block])
                position = position + 3
    return False

def seqtoprotein(seq) -> str:
    """Fully translates DNA or RNA seq into a protein."""
    transcriptedseq = "" # Stores provided sequence in RNA form
    position = 0 # Stores position in the coding strand of DNA
    CodonSeq = "" # Stores translated sequence
    if validateDNASeq(seq) == False and validateDNASeq(reversetranscription(seq)) == False:  # Check the sequence to make sure it is a valid DNA/RNA String
        return print("Invalid Sequence Provided")
    else:
        transcriptedseq = transcription(seq)
        while position < len(transcriptedseq):
            block = transcriptedseq[position:3 + position]
            # print(block)
            if len(block) < 3 or position > len(transcriptedseq):
                #print('end block found')
                return CodonSeq
            else:
                CodonSeq += RNACodonTable[block]
                # print(gencode[block])
                position = position + 3
    return CodonSeq

def randDNAstr(len):
    DNA = ''.join([random.choice(Nucleotides)
             for i in range(len)])
    return DNA

def validateDNASeq(dna_seq):
    """Check the sequence to make sure it is a DNA String. Considers "N" a valid base call"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return True

def validateRawReadSeq(dna_seq):
    """Check the sequence to make sure it is a DNA String"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in RawNucleotides:
            return False
    return True

def countNucFrequency(seq):
    """Counts the frequency of types of nucleotides"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    # "tmpFreqDict" is a Dictionary first thing in the ratio is the term : second thing is the value
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

def transcription(DNAseq):
    """ DNA to RNA transcription (RNA copy of CODING strand) """
    return DNAseq.replace("T", "U")

def reversetranscription(seq):
    """ RNA to DNA transcription """
    return seq.replace("U", "T")

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

def proteinmass(seq):
    weight = 0
    for AA in range(len(seq)): ## AA = Amino Acid
        peptide = seq[AA]
        weight += Monoisotopic_Mass[peptide]
    return weight

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
            for x in range(
                    l - 1):  # at the begining of the simulation, you need to add how many months will pass before the first adult dies of old age
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
    return

def search_phage_reads(file):
    """checks FASTA/Q file for reads then checks for the vector sequence"""

    # function to read one sequence in FASTA/Q format from a textfile.
    # input: the name (complete path) of the file that contains the sequence to read
    # output: strings containing the Identifier and the Sequence
    # this function will read ONE sequence only.  If more that one sequence is in the file, only first is read

    F = open(file)
    readlist = []  # will hold all reads
    ampliconseq = [] # holds valid vectors
    vectorseq = [] # translated inserts between vector ends
    startvector = "TTCGCAATTCCTTTAGTTGTTCCTT"
    endvector = "AGCAAAACCTCATACAGAAAATTCA"
    readcount = 0
    seqfreq = {}
    pepseq = []

    for line in F:
        line = line.rstrip("\n")  # removes new line characters (returns) from line

        if startvector and endvector in line: # if the sequence has both vectors proceed to removing flanks from sequence
            start = line.find(startvector)  # finds index for the start vector
            stop = line.find(endvector)  # finds index for the end vector
            print(line)
            readlist.append(line)
            ampliconseq.append(line[start:(stop + len(endvector))])

            removeEnd = line[(start):(stop + len(endvector))]
            pepseq.append(ampliconseqtoprotein(removeEnd))
            print(line[start:stop])
            vectorseq.append(ampliconseqtoprotein(line[start:stop+len(endvector)]))
        elif reversecomplement(startvector) and reversecomplement(endvector) in line:
            start = line.find(reversecomplement(startvector))  # finds index for the start vector
            stop = line.find(reversecomplement(endvector))  # finds index for the end vector
            readlist.append(line)
            print("reverse found")
            ampliconseq.append(reversecomplement(line[start:(stop + len(endvector))]))
            removeEnd = reversecomplement(line[(start):(stop + len(endvector))])
            pepseq.append(ampliconseqtoprotein(removeEnd))
            print(line[start:stop])
            vectorseq.append(ampliconseqtoprotein(line[start:stop + len(endvector)]))

        readcount += 1

    for item in pepseq:
        if (item in seqfreq):
            seqfreq[item] += 1
        else:
            seqfreq[item] = 1

   # for key, value in seqfreq.items():
      #  print(key, value)
    F.close()
   # print((readcount / 4))
    readcount = (readcount / 4) # takes the lines from the fastq and devides by 4 lines for each read to get the read count
    df = pd.DataFrame(seqfreq.items())
    return readlist, ampliconseq, vectorseq, readcount, seqfreq, df

def search_cluster_peptide(file, startvector, endvector):
    """take fastaptamer FASTA "-count" reads the clusters and their associated info, finds vectors."""
    # function to read fastaptamer cluster FASTA output
    # input: the name (complete path) of the file that contains the sequence to read
    # output: dicts, containing the cluster, sequence and other details about the sequence

    F = open(file)
    idcount = 0
    ampliconseq = []
    """ holds valid nucleotide vectors (inclusive) in a list"""
    clusternum = []
    """index of the cluster of duplicate seqs found"""
    #startvector = "GTTCCTTTCTATTCTCACTCGGCCGACGGGGCT"
    #endvector = "GCTGGGGCCCGCCGTGCTGGGGCCGAA"
    readcount = []
    """ Stores the number of reads found for a particular seq """
    RPM = []
    """stores the Reads per million"""

    rawread = []
    """ stores any line that is made of entirely valid nucleotides and is longer than 5 characters long """
    fastatrack = []
    """ Stores the entire ID line of VALID PEPTIDE SEQUENCES from a fasta format """
    fastaline = []
    """ Stores ALL ID lines for every read in a fasta format file """

    for line in F:
        line = line.rstrip("\n")  # removes new line characters (returns) from line

        if ">" and "-" in line:
            fastaline.append(line.upper())  # stores the "seq id" from Fasta format seq
            cluster = line[line.find(">") + 1:line.find("-")] # storing the count rank of the seq
            clusternum.append(cluster)
            # print(line[line.find("-") + 1:line.rfind("-")])
            count = int(line[line.find("-") + 1:line.rfind("-")]) # storing the number of reads of that seq
            print(count)
            readcount.append(count)
            print(count)
            ReadperM = line[line.rfind("-") + 1:len(line)]
            RPM.append(float(ReadperM)) # reads per million storage
            print(ReadperM)

        if validateRawReadSeq(line.upper()) == True and len(line) > 5: #checks for a seq line, has to be of adequate length and made of only nucleotides
            idcount +=1
            rawread.append(line.upper())

    pepseq = []
    rawindex = []
    ampToPep = {}
    ampcount = []
    rcfreq = {}
    pepcount = {}
    ampRPM = []

    for index in range(len(rawread)):
        if startvector and endvector in rawread[index]:
            start = rawread[index].find(startvector)  # finds index for the start vector
            stop = rawread[index].find(endvector)  # finds index for the end vector
            removeEnd = rawread[index][(start):(stop + len(endvector))] # isolating the nuc vector
            removeEndpep = ampliconseqtoprotein((removeEnd)) # isolating the peptide vector
            if removeEndpep in pepcount: # putting the
                # print(pepcount[removeEndpep])
                # print(readcount[index])
                pepcount[removeEndpep] = pepcount[removeEndpep] + readcount[index]
            else:
                pepcount[removeEndpep] = readcount[index]

            pepseq.append(ampliconseqtoprotein(removeEnd))
            ampliconseq.append(removeEnd)# do you need this?
            rawindex.append(index)
            ampcount.append(readcount[index]) #storing the count from fastaptamer
            ampRPM.append(RPM[index])
            fastatrack.append(fastaline[index])

        elif startvector and endvector in reversecomplement(rawread[index]): # collecting relevant sequences
            reverse = reversecomplement(rawread[index])
            start = reverse.find(startvector)  # finds index for the start vector
            stop = reverse.find(endvector)  # finds index for the end vector
            # print("reverse compliment found")
            removeEnd = reverse[(start):(stop + len(endvector))]

            try:
                removeEndpep = ampliconseqtoprotein(removeEnd)
            except KeyError:
                print("break")
                continue

            if removeEndpep in rcfreq:
                rcfreq[removeEndpep] = rcfreq[removeEndpep] + 1
            else:
                rcfreq[removeEndpep] = 1

            if removeEndpep in pepcount:
                # print(pepcount[removeEndpep])
                # print(readcount[index])
                pepcount[removeEndpep] = pepcount[removeEndpep] + readcount[index]
            else:
                pepcount[removeEndpep] = readcount[index]
            pepseq.append(removeEndpep)
            ampliconseq.append(removeEnd)  # do you need this?
            rawindex.append(index)
            ampcount.append(readcount[index])
            ampRPM.append(RPM[index])
            fastatrack.append(fastaline[index])


    pepfreq = {} #stores the freq of each peptide
    cumRPM = {} #stores the RPM
    fasta = {}
    final = []
    for index in range(len(pepseq)): #compiling stats for each seq
        ampToPep[pepseq[index]] = ampliconseq[index]
        fasta[pepseq[index]] = fastatrack[index]

        if pepseq[index] in pepfreq:
            # print(pepfreq[pepseq[index]])
            pepfreq[pepseq[index]] = (pepfreq[pepseq[index]] + ampcount[index])
            # print(pepfreq[pepseq[index]])
            cumRPM[pepseq[index]] = (cumRPM[pepseq[index]] + ampRPM[index])
        else:
            pepfreq[pepseq[index]] = ampcount[index]
            cumRPM[pepseq[index]] = ampRPM[index]

    for item in pepseq: # remove seq with less than 10 reads

        if item in pepfreq and pepfreq[item] < 3:
            del pepfreq[item]
            del cumRPM[item]
            del ampToPep[item]

    for index in range(len(pepfreq)):
        if len(pepseq[index]) > 2:
            final.append(fastatrack[index])
            final.append(pepseq[index])
        else:
            continue



    print(pepfreq)
    df = pd.DataFrame(pepfreq.items())

    df.sort_values
    df.columns = ['Peptide Sequence', 'Read Count']
    df["RPM"] = cumRPM.values()
    return df, pepfreq, final, cumRPM

def build_negative_peptide(file):
    """take fastaptamer FASTA "-count" reads the clusters and their associated info, finds vectors."""
    # function to read fastaptamer cluster FASTA output
    # input: the name (complete path) of the file that contains the sequence to read
    # output: dicts, containing the cluster, sequence and other details about the sequence

    F = open(file)
    idcount = 0
    ampliconseq = [] # holds valid vectors as
    clusternum = [] #index of the cluster of duplicate seqs found
    startvector = "TTCGCAATTCCTTTAGTTGTTCCTT"
    endvector = "AGCAAAACCTCATACAGAAAATTCA"
    readcount = []
    RPM = []
    rawread = []
    clustertrack = []
    fastaline = []
    keep = []

    for line in F:
        line = line.rstrip("\n")  # removes new line characters (returns) from line

        if ">" and "-" in line:
            count = int(line[line.find("-") + 1:line.rfind("-")])  # storing the number of reads with seq turned into an integer
            print(count)
            if count < 10:
                keep.append(False)
                continue
            else:
                keep.append(True)
                readcount.append(count)
            fastaline.append(line)
            cluster = line[line.find(">") + 1:line.find("-")] # storing the index of the cluster
            clusternum.append(cluster)
           # print(cluster)
            ReadperM = line[line.rfind("-") + 1:len(line)]
            RPM.append(float(ReadperM)) # reads per million storage
            print(ReadperM)

        if validateRawReadSeq(line.upper()) == True and len(line) > 5:
            idcount +=1
            rawread.append(line.upper())

    pepseq = []
    rawindex = []
    ampToPep = {}
    ampcount = []
    rcfreq = {}
    pepcount = {}
    ampRPM = []

    for index in range(len(rawread)):
        if keep[index] == True:
            if startvector and endvector in rawread[index]:
                start = rawread[index].find(startvector)  # finds index for the start vector
                stop = rawread[index].find(endvector)  # finds index for the end vector
                removeEnd = rawread[index][(start):(stop + len(endvector))]
                removeEndpep = ampliconseqtoprotein((removeEnd))
                if removeEndpep in pepcount:
                    # print(pepcount[removeEndpep])
                    # print(readcount[index])
                    pepcount[removeEndpep] = pepcount[removeEndpep] + readcount[index]
                else:
                    pepcount[removeEndpep] = readcount[index]

                pepseq.append(ampliconseqtoprotein(removeEnd))
                ampliconseq.append(removeEnd)# do you need this?
                rawindex.append(index)
                ampcount.append(readcount[index]) #storing the count from fastaptamer
                ampRPM.append(RPM[index])
                clustertrack.append(clusternum[index])
            elif startvector and endvector in reversecomplement(rawread[index]):
                reverse = reversecomplement(rawread[index])
                start = reverse.find(startvector)  # finds index for the start vector
                stop = reverse.find(endvector)  # finds index for the end vector
                # print("reverse compliment found")
                removeEnd = reverse[(start):(stop + len(endvector))]

                try:
                    removeEndpep = ampliconseqtoprotein(removeEnd)
                except KeyError:
                    print("break")
                    continue

                if removeEndpep in rcfreq:
                    rcfreq[removeEndpep] = rcfreq[removeEndpep] + 1
                else:
                    rcfreq[removeEndpep] = 1

                if removeEndpep in pepcount:
                    # print(pepcount[removeEndpep])
                    # print(readcount[index])
                    pepcount[removeEndpep] = pepcount[removeEndpep] + readcount[index]
                else:
                    pepcount[removeEndpep] = readcount[index]
                pepseq.append(removeEndpep)
                ampliconseq.append(removeEnd)  # do you need this?
                rawindex.append(index)
                ampcount.append(readcount[index])
                ampRPM.append(RPM[index])
                clustertrack.append(clusternum[index])

    pepfreq = {}
    cumRPM = {}

    for index in range(len(pepseq)):
        ampToPep[pepseq[index]] = ampliconseq[index]

        if pepseq[index] in pepfreq:
            # print(pepfreq[pepseq[index]])
            pepfreq[pepseq[index]] = (pepfreq[pepseq[index]] + ampcount[index])
            # print(pepfreq[pepseq[index]])
            cumRPM[pepseq[index]] = (cumRPM[pepseq[index]] + ampRPM[index])
        else:
            pepfreq[pepseq[index]] = ampcount[index]
            cumRPM[pepseq[index]] = ampRPM[index]


    print(pepfreq)
    peptidefreqDF = pd.DataFrame(pepfreq.items())

    peptidefreqDF.sort_values
    # df.columns = ['Peptide Sequence', 'Read Count']
    peptidefreqDF["RPM"] = cumRPM.values()
    backround_peptides = pepfreq.keys()
    return pepfreq, peptidefreqDF

def peptide_backround_removal(positive_selection_dict, posRPM, negative_selection_dict, pos_fastalist):
    ""
    cumRPM = posRPM
    posdict = positive_selection_dict
    negdict = negative_selection_dict

    for seq in negdict.keys():
        if seq in posdict:
            del posdict[seq]
            del cumRPM[seq]

    fasta = []

    for index in range(len(pos_fastalist)):
        if (index + 1) % 2 == 0:
            if pos_fastalist[index] in posdict:
                fasta.append(pos_fastalist[index-1])
                fasta.append(pos_fastalist[index])




    filteredFreqDf = pd.DataFrame(posdict.items())
    filteredFreqDf.columns = ['Peptide Sequence', 'Read Count']
    filteredFreqDf["RPM"] = cumRPM.values()
    return filteredFreqDf, fasta

def ampliconseqtoprotein(seq) -> str:
    """Fully translates DNA or RNA seq into a protein."""

    position = 0 # Stores position in the coding strand of DNA
    CodonSeq = "" # Stores translated sequence

    while position < len(seq):
        block = seq[position:3 + position]
        # print(block)
        if len(block) < 3 or position > len(seq):
            #print('end block found')
            return CodonSeq
        else:
            try:
                CodonSeq += DNACodonTable[block]
                # print(gencode[block])
                position = position + 3
            except KeyError:
                # print("found bad block")
                CodonSeq += "{" + block.lower() + "}"
                # print("{" + block.lower() + "}")
                # print(CodonSeq)
                position = position + 3

    return CodonSeq

def peptide_quality_check(peptideSeqList):
    trimPeptide = []
    for item in peptideSeqList:
        while (item.rfind("AGA") - item.find("GA")) > 0: # Check for valid reading frame
            if item.find("AGA") != item.rfind("AGA"):
                peptideSeqList.append(item[item.find("GA"):item.rfind("AGA")])
                item = item[item.find("GA"):item.rfind("AGA")]
            else:
                trimPeptide.append(item[:(item.rfind("AGA")+3)])
                item = item[:item.rfind("AGA")]
    return trimPeptide


    startvector = "GTTCCTTTCTATTCTCACTCGGCCGACGGGGCT"
    endvector = "GCTGGGGCCCGCCGTGCTGGGGCCGAA"


### First set of primers
startvector = "TTCGCAATTCCTTTAGTTGTTCCTT"  # original "TTCGCAATTCCTTTAGTTGTTCCTT"
# endvector = "AGCAAAACCTCATACAGAAAATTCA"

(df1, pepfreq, fasta, cumRPM) = search_cluster_peptide("15mer_lin_f3tr1_S1_L001/15mer_lin_f3tr1_S1_L001_R2_001_count_merged_trimmed.fasta", "GTTCCTTTCTATTCTCACTCGGCCGACGGGGCT", "GCTGGGGCCCGCCGTGCTGGGGCCGAA")
(df2, pepfreq2, fasta2, cumRPM2) = search_cluster_peptide("Count_15mer_linear_library_S1_L001_001_trimmed.fasta", "TTCGCAATTCCTTTAGTTGTTCCTT", "AGCAAAACCTCATACAGAAAATTCA")
vectordf = pd.concat([df1, df2], ignore_index=True)
vectordf = vectordf.sort_values(by='Read Count', ascending=False)

print(fasta)

vectordf.to_csv ('Combined_old.csv', index = False, header=True)


vectordf['Rank'] = vectordf['Read Count'].rank(ascending=False)
#vectordf.to_csv ('CountData_For_Super_Library.csv', index=False, header=True)

# Filter data for ranks 1 to 100
filtered_df = vectordf.iloc[2:]  #[vectordf['Rank'] <= 100]

# Plotting
plt.figure(figsize=(14, 9))
plt.plot(filtered_df['Rank'], filtered_df['Read Count'], marker='o', linestyle='-', color='b')
plt.yscale('log')  # Set y-axis to log scale
plt.xscale('log')
plt.xlabel('Rank')
plt.ylabel('RPM')
plt.title('RPM vs Rank (Log Scale) for 15mer_lin_f3tr1  ')
plt.grid(True)
plt.show()


#(backround, broundDF) = build_negative_peptide("Negative_Control_S1_L001_Trimmed_Merged.fasta-count.fasta")

#broundDF.to_csv('Negative_backround_peptide.csv', index = False, header=True)
#(filtertest, fasta) = peptide_backround_removal(pepfreq, cumRPM, backround, fasta)
# df = search_cluster_peptide("R1-_merged_001-count.fasta")
# # print(ampliconseqtoprotein("TNCGCAATTCCTTTAGTTGTTCCTTTCTATTCTCACTCGGCCGACGGGGCTCGTTGGTGTATATGGGTTCCGACTTCTTGGATGTGTCAGAAAGCTGGGGCC"))
# (readlist, ampliconseq, vectorseq, readcount, seqfreq, df) = search_phage_reads("R1-_merged_001.fastq")
# (seqstatdict, pepfreq, ampliconseq, ampfreq) = search_cluster_peptide("R1-_merged_001-count.fasta")
# # print(vectorseq)
# #
# # print(f"peptides:reads post QC: {len(peptide_quality_check(vectorseq)) / len(readlist)}")
# # postQC = peptide_quality_check(vectorseq)
#
#
# # print(c.values())
# # print(c.keys())
# # print(c.most_common(50))
#
# # print("Total reads: " + str(len(readlist)))
# # print("Reads with intact vector ends: " + str(len(ampliconseq)))
# # print("Total peptides found: " + str(len(vectorseq)))
# # print("Ratio of peptide:reads : " + str((len(vectorseq) / readcount)))
# print(ampfreq)
# c = collections.Counter(seqstatdict)
# print(c.most_common(50))
#

#filtertest.to_csv ('Filtered_Count_N4_Trypsin_only_S2_L001_001.csv', index = False, header=True)
#print(fasta)

#with open(r'Filtered_Count_N4_Trypsin_only_S2_L001_001.fasta', 'w') as fp:
#    fp.write('\n'.join(fasta))



# print(df)
#
#
# # c = Counter( peptideseq )
# #
# # # print(c.items())
# #
# # DNAstr1 = "ACCCAGTCTAGTTTGGTCTAGTAGTCTAGTGCGAAAACTCTCTGCGTCCGTCTAGTGTGAGTCTAGTCGTCTAGTGCCACGGTCTAGTTAGGTCTAGTCGCCGCGGTCTAGTAGTATCGTCTAGTGTCTAGTCGTCTAGTACGACTGTCTAGTTCCCGTCTAGTTGTGTCTAGTAGTCTAGTAGTCTAGTAGTCTAGTTGTCTAGTTTTATCGTTGTCTAGTGTGTCTAGTAGTCTAGTGTCTAGTACGTCTAGTAATCGTCTAGTGGCGCGCGGTCTAGTATCGGTGCGCGTCTAGTGTCTAGTAGTCTAGTGTCTAGTTGCAACTTTAGTCTAGTGAACGTGTCTAGTGTCTAGTGAAGGGTCTAGTTTGTCTAGTGTCTAGTGTCTAGTCGTTCGTCTAGTGGTCTAGTTGTCTAGTCTTGAATTAGGTCTAGTGTCTAGTGGTATGGGGTCTAGTGTCTAGTCTGGTCTAGTGGCGTCTAGTGTCTAGTGGTCTAGTGGAGTCTAGTGTCTAGTCGTCTAGTGGTCTAGTGTCCAAGTCTAGTGTCTAGTGTCTAGTGTCTAGTCGCTACCACCAGCGGTCTAGTTGTACCAGTCTAGTGTCATGCGGTCTAGTCGTCTAGTAGTCTAGTGGTCTAGTGTAGTCTAGTGTCTAGTGTCTAGTACGTCTAGTGTTGTCTAGTTTATTCGTCTAGTAAGTCTAGTCTATGTCTAGTGGGTCTAGTCTCAAGGTCTAGTTCGTCTAGTTACGATTGTCTAGTCAGTCTAGTAGTCTAGTGTCTAGTGTACGGTCTAGTAAGTCTAGTCAGTCTAGTGTCTAGTAGTGTCTAGTGTCTAGTTGTCTAGTGGTCTAGTGTCTAGTTCTGGGTCTAGTAATCGATACGAGTCTAGTGGTCTAGTATTCTGTCTAGTGTCTAGTGTCTAGTTGGGGTCTAGTTAAGTCTAGTGTCTAGT"
# # DNAstr2 = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
# # NewDNA = randDNAstr(40)
# # protein = 'MGGVFQWLVVGDAWQGGDFYCADHVGQWPEDYIWMMGMMIKWMFNMCGPNYMDVPEFCKSVWYYSVMRHYLRNLEPSACGWFRTMLVGWRGPPGHYHEEWDKVGRPSNTRMYPIHWFPENHPFYAWQSGIAHFWSTLEYNYAIPFFHDKSWWGMIWRGDGEQLCNRYKPFPLYAAAEDDPARNDYCRRSSAWKCRFCAIRHARMKWWAEVCIGMAEFWHSWMMTTRVNQESVADCAISCDGPQVVTWREQNVLAVDCFRIDESTNARSEQDCKITFYLPQWEERNAFFGHHDQRNCNKGMGFHPFAPCYTNKELMWVAFWLILCFCWQFTDYHMSPDKMIICFQFRIPQKSVFQDCIRYHKAFFRPPQKWMIRIQHRLNRHYVVENWSINSTNEDESTLHRGGLWHSEWMLLQMLLVDNDMKTHGEASRALKQCRGDNYGVNSHDSSNNPKITMYLKEFTPWHFIFPFNHMSNTRHMWHGRHLMSGMFILMLSKYWMHVSDNYCNEFNYGRMMDLRSLFMMAHKFVRIQDYYWLYKRHIWYPENGRTGAYEKYNIPAEMSKGQPVQCWTIGHKNYIALQVCNDEDSEGFIQDFWIFHFITKVMVSVFYMYISRMGIWYDPGHYTDIQHECCGCVAVPVFFLPRLILDRRKTNIDAKTGAPKFSDGFEILIRLNRSIPILIGLQFECDRKQGLPSEKGNTIGWIEVGGRDTPYHIQPNDNYQSWAAIPECTFRYKTEHWQRTCFAHHVRDTEQRGAMTQADDLEVCHVIFKDFRSHFVIYFKHEIYVLDCAKLHPHTSWLRWLGHDCIRWCLQYFSGPLDLFPEGYFANCYNNYGCRFPYWVTDHWRVMYQCFDHHSGVQQYCDQAPWKVPWVGMSNEGAPSVQSTMYREWWLHAYYYPIFFALERVVIQQYATSLLYRAYCCPRYFIQPQVRKIWMAHRGYDDAYIFRCWNNQRWVCVEGHGHALFCCDQRWWAGVADKTNNFWMYMLEWLWPANHSW'
# # exampleRNA = 'TAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
# # motif = 'GTCTAGTGT'
# # if validateDNASeq(NewDNA) == True:
# #     print("Printing Sequence, Nuc Frequency, DNA Transcript, Reverse Compliment, GC Content, and GC Content of 10 Nucleotides. ")
# #     print(f" Sequence: {NewDNA}")
# #     print(f" Nucleotide Frequency: {countNucFrequency(NewDNA)}")
# #     print(f" RNA Transcript: {transcription(NewDNA)}")
# #     print(f" Reverse Compliment: {reversecomplement(NewDNA)}")
# #     print(f" Total GC Content: {GCcontent(NewDNA)}")
# #     print(f" GC Content of Subsections of 10 Nucleotides: {SubsectionGCcontent(NewDNA, 10)}")
# #     print(f" Protein Mass of protein: {proteinmass(protein)}")
# #     print(f" Translated exampleRNA: {translation(DNAstr2 )}")
# #     print(f" Protein mRNA Inference: {mRNAinference(protein)}")
# #     print(f" Looking for motif {motif} in DNAStr1: {motifFinder(DNAstr1, motif)}")
# #     print(f" Looking for Open Reading Frames in file: {listopenreadingframes('rosalind_orf.txt')}")
# #     print(randDNAStr)
# #     print(SubsectionGCcontent(NewDNA, 10))
# print(f" {search_cluster_peptide('R1-_merged_001-cluster.fasta')} ")
# print(pepfreq)
# c = collections.Counter(ampliconseq)
# print(c.most_common(50))
# #     print(f" {LifespanFibonacci(94, 19)} ")
# #     print(HighestGCcontent('input.txt'))
# # print(str(PointMutationCounter(DNAstr1, DNAstr2)))
# # print(Permutations(6))
