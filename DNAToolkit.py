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
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'Stop', 'TGG': 'W'}

RNACodonTable = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L',  # key: Value
                 'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
                 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 'UCU': 'S', 'CCU': 'P',
                 'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 'CCG': 'P',
                 'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
                 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': 'Stop', 'CAA': 'Q',
                 'AAA': 'K', 'GAA': 'E', 'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                 'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R',
                 'AGC': 'S', 'GGC': 'G', 'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

ProteinCodonTable = {'F': ['UUU', 'UUC'], 'L': ['CUU', 'CUC', 'UUA', 'CUA', 'UUG', 'CUG'], 'I': ['AUU', 'AUC', 'AUA'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                     'M': ['AUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'],
                     'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'Y': ['UAU', 'UAC'], 'H': ['CAU', 'CAC'], 'N': ['AAU', 'AAC'], 'D': ['GAU', 'GAC'], 'Stop': ['UAA', 'UAG', 'UGA'],
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
    """Translates DNA or RNA seq into a protein stopping at the first 'Stop' codon."""
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
            if len(block) < 3 or RNACodonTable[block] == "Stop":
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
    """Check the sequence to make sure it is a DNA String"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
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

DNAstr1 = "ACCCAGTCTAGTTTGGTCTAGTAGTCTAGTGCGAAAACTCTCTGCGTCCGTCTAGTGTGAGTCTAGTCGTCTAGTGCCACGGTCTAGTTAGGTCTAGTCGCCGCGGTCTAGTAGTATCGTCTAGTGTCTAGTCGTCTAGTACGACTGTCTAGTTCCCGTCTAGTTGTGTCTAGTAGTCTAGTAGTCTAGTAGTCTAGTTGTCTAGTTTTATCGTTGTCTAGTGTGTCTAGTAGTCTAGTGTCTAGTACGTCTAGTAATCGTCTAGTGGCGCGCGGTCTAGTATCGGTGCGCGTCTAGTGTCTAGTAGTCTAGTGTCTAGTTGCAACTTTAGTCTAGTGAACGTGTCTAGTGTCTAGTGAAGGGTCTAGTTTGTCTAGTGTCTAGTGTCTAGTCGTTCGTCTAGTGGTCTAGTTGTCTAGTCTTGAATTAGGTCTAGTGTCTAGTGGTATGGGGTCTAGTGTCTAGTCTGGTCTAGTGGCGTCTAGTGTCTAGTGGTCTAGTGGAGTCTAGTGTCTAGTCGTCTAGTGGTCTAGTGTCCAAGTCTAGTGTCTAGTGTCTAGTGTCTAGTCGCTACCACCAGCGGTCTAGTTGTACCAGTCTAGTGTCATGCGGTCTAGTCGTCTAGTAGTCTAGTGGTCTAGTGTAGTCTAGTGTCTAGTGTCTAGTACGTCTAGTGTTGTCTAGTTTATTCGTCTAGTAAGTCTAGTCTATGTCTAGTGGGTCTAGTCTCAAGGTCTAGTTCGTCTAGTTACGATTGTCTAGTCAGTCTAGTAGTCTAGTGTCTAGTGTACGGTCTAGTAAGTCTAGTCAGTCTAGTGTCTAGTAGTGTCTAGTGTCTAGTTGTCTAGTGGTCTAGTGTCTAGTTCTGGGTCTAGTAATCGATACGAGTCTAGTGGTCTAGTATTCTGTCTAGTGTCTAGTGTCTAGTTGGGGTCTAGTTAAGTCTAGTGTCTAGT"
DNAstr2 = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
NewDNA = randDNAstr(40)
protein = 'MGGVFQWLVVGDAWQGGDFYCADHVGQWPEDYIWMMGMMIKWMFNMCGPNYMDVPEFCKSVWYYSVMRHYLRNLEPSACGWFRTMLVGWRGPPGHYHEEWDKVGRPSNTRMYPIHWFPENHPFYAWQSGIAHFWSTLEYNYAIPFFHDKSWWGMIWRGDGEQLCNRYKPFPLYAAAEDDPARNDYCRRSSAWKCRFCAIRHARMKWWAEVCIGMAEFWHSWMMTTRVNQESVADCAISCDGPQVVTWREQNVLAVDCFRIDESTNARSEQDCKITFYLPQWEERNAFFGHHDQRNCNKGMGFHPFAPCYTNKELMWVAFWLILCFCWQFTDYHMSPDKMIICFQFRIPQKSVFQDCIRYHKAFFRPPQKWMIRIQHRLNRHYVVENWSINSTNEDESTLHRGGLWHSEWMLLQMLLVDNDMKTHGEASRALKQCRGDNYGVNSHDSSNNPKITMYLKEFTPWHFIFPFNHMSNTRHMWHGRHLMSGMFILMLSKYWMHVSDNYCNEFNYGRMMDLRSLFMMAHKFVRIQDYYWLYKRHIWYPENGRTGAYEKYNIPAEMSKGQPVQCWTIGHKNYIALQVCNDEDSEGFIQDFWIFHFITKVMVSVFYMYISRMGIWYDPGHYTDIQHECCGCVAVPVFFLPRLILDRRKTNIDAKTGAPKFSDGFEILIRLNRSIPILIGLQFECDRKQGLPSEKGNTIGWIEVGGRDTPYHIQPNDNYQSWAAIPECTFRYKTEHWQRTCFAHHVRDTEQRGAMTQADDLEVCHVIFKDFRSHFVIYFKHEIYVLDCAKLHPHTSWLRWLGHDCIRWCLQYFSGPLDLFPEGYFANCYNNYGCRFPYWVTDHWRVMYQCFDHHSGVQQYCDQAPWKVPWVGMSNEGAPSVQSTMYREWWLHAYYYPIFFALERVVIQQYATSLLYRAYCCPRYFIQPQVRKIWMAHRGYDDAYIFRCWNNQRWVCVEGHGHALFCCDQRWWAGVADKTNNFWMYMLEWLWPANHSW'
exampleRNA = 'TAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
motif = 'GTCTAGTGT'
if validateDNASeq(NewDNA) == True:
    print("Printing Sequence, Nuc Frequency, DNA Transcript, Reverse Compliment, GC Content, and GC Content of 10 Nucleotides. ")
    print(f" Sequence: {NewDNA}")
    print(f" Nucleotide Frequency: {countNucFrequency(NewDNA)}")
    print(f" RNA Transcript: {transcription(NewDNA)}")
    print(f" Reverse Compliment: {reversecomplement(NewDNA)}")
    print(f" Total GC Content: {GCcontent(NewDNA)}")
    print(f" GC Content of Subsections of 10 Nucleotides: {SubsectionGCcontent(NewDNA, 10)}")
    print(f" Protein Mass of protein: {proteinmass(protein)}")
    print(f" Translated exampleRNA: {translation(DNAstr2 )}")
    print(f" Protein mRNA Inference: {mRNAinference(protein)}")
    print(f" Looking for motif {motif} in DNAStr1: {motifFinder(DNAstr1, motif)}")
    print(f" Looking for Open Reading Frames in file: {listopenreadingframes('rosalind_orf.txt')}")
#     print(randDNAStr)
#     print(SubsectionGCcontent(NewDNA, 10))
#     print(f" {read_one_fasta('Amel_DSCAM.FASTA')} ")
#     print(f" {LifespanFibonacci(94, 19)} ")
#     print(HighestGCcontent('input.txt'))
# print(str(PointMutationCounter(DNAstr1, DNAstr2)))
# print(Permutations(6))
