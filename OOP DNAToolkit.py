# DNA Toolkit file
import collections
import random
import re
import itertools
import math
from collections import Counter
import csv
import pandas as pd

Nucleotides = ["A", "C", "G", "T"]
RNAnuc = ["A", "C", "G", "U"]
DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# Creating a random DNA sequence for testing:
randDNAStr = ''.join([random.choice(Nucleotides)
                      for nuc in range(100)])

RNACodonTable = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L',
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
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '{Stop}', 'TAG': '{Stop}',
            'TGC': 'C', 'TGT': 'C', 'TGA': '{Stop}', 'TGG': 'W'}



class fasta():
    """takes one FASTA file as an input for several types of analysis"""
    filename = ""

    def __init__(self, file):
        self.filename = file
        self.idnum = 0 # tracks of how many ID lines have been found
        self.seqID = []  # holds Identifier (name) of sequence
        self.seq = ""
        self.seqlist = []  # holds Sequences
        self.percentGC = []  # holds percent GC
        self.nucfreq = []  # holds a dictionary with the count of every nucleotide in a sequence

        F = open(file)
        idnum = 0
        ID = []
        seq = ""
        seqlist = []

        for line in F:
            line = line.rstrip("\n")  # removes new line characters (returns) from line
            m = re.search(r"^>(\S+)", line)  #  checks if line has '>'
            if m and idnum > 0:
                ID.append(line)  # remember ID Name
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
        F.close()
        self.idnum = idnum
        self.seqID = ID
        self.seqlist.append(seq)
        self.seq = seq


        # This will get GC content for each seq in the FASTA file
        for seq in seqlist:
            GCFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
            GC = 0
            for nuc in seq:
                GCFreqDict[nuc] += 1
            GC += GCFreqDict["G"]
            GC += GCFreqDict["C"]
            self.percentGC.append(GC / len(seq))

        # Counts the frequency of types of nucleotides
        for seq in self.seqlist:
            nucFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
            for nuc in seq:
                nucFreqDict[nuc] += 1
            self.percentGC.append(nucFreqDict)

    def validateSeq(self):
        """Check the sequence to make sure it is a DNA String"""
        for nuc in self.seq:
            if nuc in Nucleotides:
                pass
            else:
                return False
        return True

    def GCcontent(self):
        """Gives percent of G and C base pairs sequence"""
        if self.validateSeq() == True:
            tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
            GC = 0
            for nuc in self.seq:
                tmpFreqDict[nuc] += 1
            GC = tmpFreqDict["G"] + tmpFreqDict["C"]
            print(tmpFreqDict["G"])
            print(tmpFreqDict["C"])
            GC = GC / len(self.seq)
            return GC

    def transcription(self):
        """ DNA to RNA transcription (RNA copy of CODING strand) """
        self.RNASeq = self.seq.replace("T", "U")
        return self.RNASeq

    def translation(self):
        """Translates DNA or RNA seq into a protein starting at the first 'M' and stopping at the first 'Stop' codon."""
        transcriptedseq = ""  # Stores provided sequence in RNA form
        position = 0  # Stores position in the coding strand of DNA
        CodonSeq = ""  # Stores translated sequence
        encode = True
        transcriptedseq = self.transcription()
        while position < len(transcriptedseq):
            block = transcriptedseq[position:3 + position]
            # print(block)
            if RNACodonTable[block] == 'M':
                encode = True
            if encode == True:
                if len(block) < 3 or RNACodonTable[block] == "Stop":
                    # print('end block found')
                    self.translatedseq = CodonSeq
                    return CodonSeq
                else:
                    CodonSeq += RNACodonTable[block]
                    # print(gencode[block])
                    position = position + 3

        return CodonSeq

class fastaptamer():
    """takes one FASTQ file as an input for several types of analysis"""
    filename = ""


    def __init__(self, file):
        """collects attributes and begins to parse a fastaptamer count file"""
        self.filename = file
        self.idnum = 0 # tracks of how many ID lines have been found
        self.seqlist = []  # holds Sequences
        self.percentGC = []  # holds percent GC
        self.nucfreq = []  # holds a dictionary with the count of every nucleotide in a sequence
        self.rawread =[]
        self.clusternum = []
        self.readcount = []
        """occurance count for each sequence found"""
        self.RPM = []
        """Reads Per Million for each sequence found"""
        self.RawNucleotides = ["A", "C", "G", "T", "N"]
        self.Nucleotides = ["A", "C", "G", "T"]
        self.DNA_ReverseComplement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        self.DNACodonTable = {
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
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '{Stop}', 'TAG': '{Stop}',
            'TGC': 'C', 'TGT': 'C', 'TGA': '{Stop}', 'TGG': 'W'}
        self.RNACodonTable = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 'UUC': 'F', 'CUC': 'L',  # key: Value
                         'AUC': 'I', 'GUC': 'V', 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
                         'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 'UCU': 'S', 'CCU': 'P',
                         'ACU': 'T', 'GCU': 'A', 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                         'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'UCG': 'S', 'CCG': 'P',
                         'ACG': 'T', 'GCG': 'A', 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
                         'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'UAA': '{Stop}', 'CAA': 'Q',
                         'AAA': 'K', 'GAA': 'E', 'UAG': '{Stop}', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                         'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 'UGC': 'C', 'CGC': 'R',
                         'AGC': 'S', 'GGC': 'G', 'UGA': '{Stop}', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                         'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
        self.ProteinCodonTable = {'F': ['UUU', 'UUC'], 'L': ['CUU', 'CUC', 'UUA', 'CUA', 'UUG', 'CUG'],
                             'I': ['AUU', 'AUC', 'AUA'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                             'M': ['AUG'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
                             'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'],
                             'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'Y': ['UAU', 'UAC'], 'H': ['CAU', 'CAC'],
                             'N': ['AAU', 'AAC'], 'D': ['GAU', 'GAC'], '{Stop}': ['UAA', 'UAG', 'UGA'],
                             'Q': ['CAA', 'CAG'], 'K': ['AAA', 'AAG'], 'E': ['GAA', 'GAG'], 'C': ['UGU', 'UGC'],
                             'R': ['CGU', 'CGC', 'CGA', 'AGA', 'CGG', 'AGG'],
                             'G': ['GGU', 'GGC', 'GGA', 'GGG'], 'W': ['UGG']}
        F = open(file)
        index = 0

        for line in F:
            line = line.rstrip("\n")  # removes new line characters (returns) from line
            if ">" and "-" in line:  #  checks if line has '>'
                cluster = line[line.find(">") + 1:line.find("-")]  # storing the index of the cluster
                self.clusternum.append(cluster)
                # print(cluster)
                count = int(
                    line[line.find("-") + 1:line.find("-", line.find("-") + 1)])  # storing the number of reads with seq
                self.readcount.append(count)
                #   print(count)
                ReadperM = line[
                           line.find("-", line.find("-") + 2):line.find("-", line.find("-", line.find("-", line.find(
                               "-") + 1)) + 1)]
                self.RPM.append(abs(float(ReadperM)))  # reads per million storage

            if self.validateRawReadSeq(line) == True and len(line) > 2:  # if the sequence is a valid DNA sequence and the length is valid
                self.rawread.append(line)


            index +=1
        F.close()
        self.clusterstatdict = {'Cluster Rank': self.clusternum,
                                'Read Frequency': self.readcount,
                                'Reads Per Million': self.RPM,
                                'Read': self.rawread}


        self.RawReadStats = pd.DataFrame.from_dict(self.clusterstatdict, orient='index').transpose()
        return print(self.RawReadStats)


    def validateRawReadSeq(self,seq):
        """Check the sequence to make sure it is a DNA String"""
        tmpseq = seq.upper()
        for nuc in tmpseq:
            if nuc not in self.RawNucleotides:
                return False
        return True

    def reversecomplement(self,seq):
        return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]  # returns 5' to 3'

    def ampliconseqtoprotein(self,seq) -> str:
        """Fully translates DNA or RNA seq into a protein."""

        position = 0  # Stores position in the coding strand of DNA
        CodonSeq = ""  # Stores translated sequence

        while position < len(seq):
            block = seq[position:3 + position]
            # print(block)
            if len(block) < 3 or position > len(seq):
                # print('end block found')
                return CodonSeq
            else:
                try:
                    CodonSeq += DNACodonTable[block]
                    # print(gencode[block])
                    position = position + 3
                except KeyError:
                    #print("found bad block")
                    CodonSeq += "{" + block.lower() + "}"
                    # print("{" + block.lower() + "}")
                    # print(CodonSeq)
                    position = position + 3
        return CodonSeq

    def search_cluster_peptide(self):
        """take fastaptamer FASTA "-count" reads the clusters and their associated info, finds vectors."""

        # function to read fastaptamer cluster FASTA output
        # input: the name (complete path) of the file that contains the sequence to read
        # output: dataframe with the peptide sequence, Readcount, and vectors found for it

        # next thing to do is use the dataframe made of all the raw stats and see if you can add a column next to it with the vector info

        readlist = []
        """holds all the reads from the fastq """

        ampliconseq = []
        """ holds the confirmed vector sequence from a read"""
        vectorcount = 0
        startvector = "TCCTTTCTATTCTCACTCGGCCGACGGGGCT"
        endvector = "GCTGGGGCCCGCCGTGCTGGGGCCGAA"
        vectorlist=[]
        ampcount = {}
        ampliconcount = []
        vectordict = {}
        RPM={}
        """tracks vectors that map to a peptide"""

        for index in range(len(self.rawread)):

            line = self.rawread[index].rstrip("\n")  # removes new line characters (returns) from line
            start = line.find(startvector)  # finds index for the start vector
            stop = line.find(endvector)  # finds index for the end vector
            rcstart = line.find(self.reversecomplement(startvector))  # finds index for the start vector
            rcstop = line.find(self.reversecomplement(endvector))


            if startvector and endvector in line:  # if the sequence has both vectors proceed to removing flanks from sequence
                # print(line)

                removeEnd = line[(start):(stop + len(endvector))]
                if "N" in removeEnd:
                    continue
                readlist.append(line)
                vectorlist.append(removeEnd)
                amp = self.ampliconseqtoprotein(removeEnd)
                ampliconseq.append(self.ampliconseqtoprotein(removeEnd))

                # print(removeEnd)
                if (removeEnd in ampcount) == False:
                    ampcount[removeEnd] = self.readcount[index]
                    vectordict[amp]=[] #tracks vectors that map to a peptide
                    vectordict[amp].append(removeEnd)
                    RPM[amp] = self.RPM[index]
                else:
                    ampcount[removeEnd] = ampcount[removeEnd] + self.readcount[index]
                    RPM[amp] = RPM[amp] + self.RPM[index]
                    vectordict[amp].append(removeEnd)


                vectorcount += 1
                ampliconcount.append(ampcount[removeEnd])
            elif self.reversecomplement(startvector) and self.reversecomplement(endvector) in line:
                # print(line)
                removeEnd = line[(rcstart + len(startvector)):(rcstop + len(endvector))]
                if "N" in removeEnd:
                    continue

                readlist.append(line)
                amp = self.ampliconseqtoprotein(self.reversecomplement(removeEnd))
                ampliconseq.append(
                    self.ampliconseqtoprotein(line[rcstart:(rcstop + len(endvector))]))  # do you need this?

                vectorlist.append(removeEnd)
                # print(removeEnd)
                if (self.reversecomplement(removeEnd) in ampcount) == False:
                    ampcount[self.reversecomplement(removeEnd)] = self.readcount[index]
                    vectordict[amp]=[] #tracks vectors that map to a peptide
                    vectordict[amp].append(removeEnd)
                    RPM[amp] = self.RPM[index]
                else:
                    ampcount[self.reversecomplement(removeEnd)] = ampcount[self.reversecomplement(removeEnd)] + self.readcount[index]
                    vectordict[amp].append(removeEnd)
                    RPM[amp] = self.RPM[index] + RPM[amp]
                vectorcount += 1
                ampliconcount.append(ampcount[self.reversecomplement(removeEnd)])

        ampfreq = {}

        vectordf = pd.DataFrame()
        for peptide in vectordict:
            tempveclist = list(set(vectordict[peptide]))
            #print(tempveclist)
            vectors = pd.DataFrame([tempveclist], columns=[f'Column {i + 1}' for i in range(len(tempveclist))])
            newcol = pd.DataFrame([peptide], columns=['peptide'])
            finalvectors = pd.concat([newcol, vectors], axis=1)
            vectordf = pd.concat([finalvectors, vectordf], axis=0)

        vectordf.columns = ['Peptide Seq', 'Seq 1', 'Seq 2', 'Seq 3']
        print(vectordf)
        for index in range(len(ampliconseq)):
            ampfreq[ampliconseq[index]] = ampliconcount[index]
        df = pd.DataFrame(ampfreq.items())
        RPMdf = pd.DataFrame(RPM.items())
        RPMdf.columns = ['Peptide Seq', 'RPM']
        df.columns = ['Peptide Seq', 'Read Count']
        dfmerge = pd.merge(vectordf, df, on='Peptide Seq') # Merged DF that contains the peptide, all of the vectors found to map to it, and the count containing the occurances
        dfmerge = pd.merge(dfmerge, RPMdf, on='Peptide Seq')
        print(dfmerge)
        self.peptides = dfmerge.sort_values('Read Count', ascending=False)
        return self.peptides


    def backround_peptide(self,num):
        """take fastaptamer FASTA "-count" runs a search for valid peptides. Only keeps the samples a read number larger than the input you give. Outputs a pandas dataframe"""

        # function to read fastaptamer cluster FASTA output
        # input: the name (complete path) of the file that contains the sequence to read
        # output: dataframe with the peptide sequence, Readcount, and vectors found for it


        readlist = []
        """holds all the reads from the fastq """

        ampliconseq = []
        """ holds the confirmed vector sequence from a read"""
        vectorcount = 0
        startvector = "TTCGCAATTCCTTTAGTTGTTCCTT"
        endvector = "AGCAAAACCTCATACAGAAAATTCA"
        vectorlist=[]
        ampcount = {}
        ampliconcount = []
        vectordict = {}
        RPM={}
        """tracks vectors that map to a peptide"""

        for index in range(len(self.rawread)):

            line = self.rawread[index].rstrip("\n")  # removes new line characters (returns) from line
            start = line.find(startvector)  # finds index for the start vector
            stop = line.find(endvector)  # finds index for the end vector
            rcstart = line.find(self.reversecomplement(startvector))  # finds index for the start vector
            rcstop = line.find(self.reversecomplement(endvector))


            if startvector and endvector in line:  # if the sequence has both vectors proceed to removing flanks from sequence
                # print(line)

                removeEnd = line[(start):(stop + len(endvector))]
                if "N" in removeEnd:
                    continue
                readlist.append(line)
                vectorlist.append(removeEnd)
                amp = self.ampliconseqtoprotein(removeEnd)
                ampliconseq.append(self.ampliconseqtoprotein(removeEnd))

                # print(removeEnd)
                if (removeEnd in ampcount) == False:
                    ampcount[removeEnd] = self.readcount[index]
                    vectordict[amp]=[] #tracks vectors that map to a peptide
                    vectordict[amp].append(removeEnd)
                    RPM[amp] = self.RPM[index]
                else:
                    ampcount[removeEnd] = ampcount[removeEnd] + self.readcount[index]
                    RPM[amp] = RPM[amp] + self.RPM[index]
                    vectordict[amp].append(removeEnd)


                vectorcount += 1
                ampliconcount.append(ampcount[removeEnd])
            elif self.reversecomplement(startvector) and self.reversecomplement(endvector) in line:
                # print(line)
                removeEnd = line[(rcstart + len(startvector)):(rcstop + len(endvector))]
                if "N" in removeEnd:
                    continue

                readlist.append(line)
                amp = self.ampliconseqtoprotein(self.reversecomplement(removeEnd))
                ampliconseq.append(
                    self.ampliconseqtoprotein(line[rcstart:(rcstop + len(endvector))]))  # do you need this?

                vectorlist.append(removeEnd)
                # print(removeEnd)
                if (self.reversecomplement(removeEnd) in ampcount) == False:
                    ampcount[self.reversecomplement(removeEnd)] = self.readcount[index]
                    vectordict[amp]=[] #tracks vectors that map to a peptide
                    vectordict[amp].append(removeEnd)
                    RPM[amp] = self.RPM[index]
                else:
                    ampcount[self.reversecomplement(removeEnd)] = ampcount[self.reversecomplement(removeEnd)] + self.readcount[index]
                    vectordict[amp].append(removeEnd)
                    RPM[amp] = self.RPM[index] + RPM[amp]
                vectorcount += 1
                ampliconcount.append(ampcount[self.reversecomplement(removeEnd)])

        ampfreq = {}

        vectordf = pd.DataFrame()
        for peptide in vectordict:
            tempveclist = list(set(vectordict[peptide]))
            #print(tempveclist)
            vectors = pd.DataFrame([tempveclist], columns=[f'Column {i + 1}' for i in range(len(tempveclist))])
            newcol = pd.DataFrame([peptide], columns=['peptide'])
            finalvectors = pd.concat([newcol, vectors], axis=1)
            vectordf = pd.concat([finalvectors, vectordf], axis=0)

        vectordf.columns = ['Peptide Seq', 'Seq 1', 'Seq 2', 'Seq 3']
        print(vectordf)
        for index in range(len(ampliconseq)):
            ampfreq[ampliconseq[index]] = ampliconcount[index]
        df = pd.DataFrame(ampfreq.items())
        RPMdf = pd.DataFrame(RPM.items())
        RPMdf.columns = ['Peptide Seq', 'RPM']
        df.columns = ['Peptide Seq', 'Read Count']
        dfmerge = pd.merge(vectordf, df, on='Peptide Seq') # Merged DF that contains the peptide, all of the vectors found to map to it, and the count containing the occurances
        dfmerge = pd.merge(dfmerge, RPMdf, on='Peptide Seq')
        print(dfmerge)
        dfmerge = dfmerge.sort_values('Read Count', ascending=False)

        for index, row in dfmerge.iterrows():
            if row['Read Count'] < num:
                dfmerge.drop(index, inplace=True)

        return dfmerge

    def backround_filter(self,Nselection):
        """Takes a negative selection dataframe and uses it as a backround for removing peptides from a dataframe of peptides """

        df = self.peptides
        backrounddict = {}
        for index, row in Nselection.iterrows():
            df = df.drop(df['Peptide Seq'] == row['Peptide Seq'])
        print(df)

        return df


merged_001 = fastaptamer("R1-_merged_001-count.fasta")
vectordf = merged_001.search_cluster_peptide()
merged_001.backround_peptide(20)
vectordf.to_csv ('vectordf.csv', index=False, header=True)
