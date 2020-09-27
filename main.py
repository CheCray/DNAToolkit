# DNA Toolset/Code testing file
from DNAToolkit import *
import random

# Creating a random DNA sequence for testing:
randDNAStr = ''.join([random.choice(Nucleotides)
                      for nuc in range(20)])

DNAStr = randDNAStr   #validateSeq(randDNAStr)


print(f"[1.] Sequence {DNAStr}\n")
print(f'[2.] Nucleotide Frequency: {countNucFrequency(DNAStr)}\n')
print(f"[3.] Transcribed Sequence { (DNAStr)}\n")

print(f"[4.] Reverse Complement Sequence: \n\n  Original Seq: 5' {DNAStr} '3")
print(f"                   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"Complement Seq: 3' {reversecomplement(DNAStr)[::-1]} '5")
print(f"Complement Seq: 5' {reversecomplement(DNAStr)} '3")
print(f"[5.] GC Content:  {GCcontent(DNAStr)}")
