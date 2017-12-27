# RUN SCRIPT ON PYTHON 2.7
"""
kmerDNA.py screens a continguous DNA sequence for 20-mer sgRNA that cuts 
at two distinct low copy repeats within the sequence to introduce a microdeletion
(see scoreGuideDesign method below for more information).
It takes in a single DNA sequence in FASTA format and outputs a CSV file containing
the candidate gRNA targets, their the location in the sequence, and size of deletion.

EXAMPLE:
scoreGuideDesign("22q11.21 hg38 UCSC.fasta")
    Requires a single DNA sequence in FASTA format and prints a csv file with
    possible sgRNAs for dual cutting to introduce microdeletions (Tai et al.
    Nature neuroscience 2016)

EXAMPLE:
scoreGuideDesign("22q11.21 hg19 UCSC.fasta", minSize = 2330000, maxSize = 3500000, minDupA = 560000, maxDupA = 1130000, minDupB = 3460000, maxDupB = 4020000)

EXAMPLE:
kmer = kmer_count("22q11.21 hg38 UCSC.fasta", 18)
kmer = nplicate(kmer, 2)
printCSV(kmer, "2")
    where "22q11.21 hg38 UCSC.fasta" contains a single DNA sequence from
    which 18-mers will be identified and counted and the resulted printed as a csv
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from collections import Counter
import csv


def read_fasta(file_name):
    # Read a fasta file (can contain only one sequence) and return Bio.SeqRecord object
    return SeqIO.read(file_name, "fasta", alphabet = IUPAC.ambiguous_dna)


def kmer_count(file_name, kmer, minRange = 0, maxRange = -1):
    # Take in a fasta file with one DNA sequence and a value for length of k-mer
    
    # Read the sequence along with its reverse complement
    DNAseq = read_fasta(file_name).seq
    if maxRange == -1:
        DNAseq = DNAseq[minRange:]
    else:
        DNAseq = DNAseq[minRange:maxRange]

    DNAseqRC = DNAseq.reverse_complement()

    # Count kmers
    store = {}
    for k in range(len(DNAseq) - kmer + 1):
        # Subset the sequence for kmer and keep count as value in store
        oligo1 = str(DNAseq[k : k + kmer])
        store[oligo1] = store.get(oligo1, 0) + 1

        oligo2 = str(DNAseqRC[k : k + kmer])
        store[oligo2] = store.get(oligo2, 0) + 1

    # Return a list of tuples ('sequence', count)
    return store.items()
    

def nplicate(kmer_list, Nplicate_filter = 2):
    # Takes in return value from kmer_count method and filters for duplicates by default
    # Optional argument: Nplicate_filter takes a positive integer to filter the kmer by frequency of occurrence
    ret_list = []
    
    for kmer in kmer_list:
        if kmer[1] == Nplicate_filter:
            ret_list.append(kmer)

    # Returns a list of tuples('sequence', count)
    return ret_list
    

def printCSV(kmer_list, Nplicate = "N", dscrpt = ""):
    # Takes in return value from kmer_count method
    # Nplicate: to specify whether nplicate method was applied
    # dscrpt: for a more descriptive title
    
    kmer = len(kmer_list[0][0])
    file_name = str(kmer) + "mer_" + str(Nplicate) + "plicate" + dscrpt + ".csv"

    # Write CSV file
    csvwrite = open(file_name, "wb")
    csv_writer = csv.writer(csvwrite, delimiter = ',')
    #csv_writer.writerow(["Sequence", "Count"])

    for row in kmer_list:
        if not "N" in row[0]:
            # Conditional removes any N bases from the output
            csv_writer.writerow(row)

    csvwrite.close()

    print(file_name + " PRINTED.")
       

def scoreGuideDesign(file_name,
                     minSize = 1970000, maxSize = 4000000,
                     minDupA = 580000, maxDupA = 1640000,
                     minDupB = 3610000, maxDupB = 4170000,
                     guideLength = 20):
    # CRISPR design with 20-bp target with a 5'-G and ends with 5'-NGG-3' (pX330 or pX459 design)
    # for generating microdeletions using SCORE approach (Tai et al. Nature neuroscience 2016)
    # file_name: fasta file with a single DNA sequence
    # minSize: Acceptance limit for minimum distance between targets (defaults for 22q11.21 hg38 UCSC.fasta)
    # maxSize: Acceptance limit for maximum distance between targets (defaults for 22q11.21 hg38 UCSC.fasta)
    # min/max-DupA/DupB = lower and upper ranges for duplication region A and B (defaults for 22q11.21DS)
    
    DNAseq = read_fasta(file_name).seq
    DNAseqRC = DNAseq.reverse_complement()
    DNAseq = str(DNAseq)
    DNAseqRC = str(DNAseqRC)


    def seqListOnly(kmer):
        # Converts output from kmer_count and nplicate into a single list of UPPERCASE kmers
        return [oligo[0].upper() for oligo in kmer]


    def fivePrimeGandTTTTfilter(oligo):
        # Filters DNA targets for 5'-G and eliminate guides containing "TTTT" for U6 promoter transcriptional compatibility
        # Takes in a list of sequences (processed by seqListOnly method) as argument
        ret_list = []

        for target in oligo:
            if target[0] == "G" and "TTTT" not in target:
                ret_list.append(target)

        # Returns a list of sequences with a leading 5'-G base
        return ret_list


    def findPosition(oligo):
        # Given a sequence, finds the position on DNAseq and DNAseqRC
        # + position = target on (+)-sense DNA
        # - position = target on (-)-sense DNA
        pos = DNAseq.find(oligo)
        neg = DNAseqRC.find(oligo)

        if pos == -1:
            return (-1 * neg, -1 * DNAseqRC.find(oligo, neg + 1))
        elif neg == -1:
            return (pos, DNAseq.find(oligo, pos + 1))
        else:
            return (pos, -1 * neg)


    def NGGfilter(oligo):
        # Filters DNA targets for NGG PAM downstream
        # Takes in a list of "sequence" as argument
        ret_list = []
        init, end = guideLength + 1, guideLength + 3

        for seq in oligo:
            pos1, pos2 = findPosition(seq)
            
            if pos1 > 0:
                if pos2 > 0:
                    if DNAseq[pos1 + init:pos1 + end] == "GG" and DNAseq[pos2 + init:pos2 + end] == "GG":
                        ret_list.append(seq)
                else:
                    if DNAseq[pos1 + init:pos1 + end] == "GG" and DNAseqRC[-pos2 + init:-pos2 + end] == "GG":
                        ret_list.append(seq)
            else:
                if pos2 > 0:
                    if DNAseqRC[-pos1 + init:-pos1 + end] == "GG" and DNAseq[pos2 + init:pos2 + end] == "GG":
                        ret_list.append(seq)
                else:
                    if DNAseqRC[-pos1 + init:-pos1 + end] == "GG" and DNAseqRC[-pos2 + init:-pos2 + end] == "GG":
                        ret_list.append(seq)

        # Returns a list of sequences filtered for the presence of NGG PAM
        return ret_list


    def distanceBetweenTargets(oligo):
        # Lays out all the positions of the targets and their distance apart
        ret_list = []

        for seq in oligo:
            pos1, pos2 = findPosition(seq)
            if pos1 > 0:
                if pos2 > 0:
                    ret_list.append((seq, pos1, pos2, abs(pos1 - pos2)))
                else:
                    ret_list.append((seq, pos1, pos2, abs(pos1 - (len(DNAseq) + pos2))))
            else:
                if pos2 > 0:
                    ret_list.append((seq, pos1, pos2, abs(pos2 - (len(DNAseq) + pos1))))
                else:
                    ret_list.append((seq, pos1, pos2, abs(pos1 - pos2)))

        # Returns a list of tuple ("sequence", target position1, target position2, distance between targets)
        return ret_list


    def delSizeFilter(dist_data):
        # Checks that size of the deletions is within set limits (minSize/maxSize)
        ret_list = []

        for data in dist_data:
            if minSize <= data[3] <= maxSize:
                ret_list.append(data)

        return ret_list


    def segDupFilter(dist_data):
        # Checks that the targets are within the provided segmental duplication regions (min/maxDupA/B)
        ret_list = []

        for data in dist_data:
            if data[1] > 0:
                pos1 = data[1]
            else:
                pos1 = len(DNAseq) + data[1]
            if data[2] > 0:
                pos2 = data[2]
            else:
                pos2 = len(DNAseq) + data[2]
                
            if (minDupA <= pos1 <= maxDupA and minDupB <= pos2 <= maxDupB) or (minDupA <= pos2 <= maxDupA and minDupB <= pos1 <= maxDupB):
                ret_list.append(data)
                
        # Returns a list of tuple ("sequence", target position1, target position2, distance between targets)
        return ret_list

    


    kmer = kmer_count(file_name, guideLength)
    print("Step 1/8 - kmer counting COMPLETE.")

    twoPlicate = nplicate(kmer, 2)
    print("Step 2/8 - filter for duplicates COMPLETE.")

    seq = seqListOnly(twoPlicate)
    print("Step 3/8 - clean data COMPLETE.")

    Gseq = fivePrimeGandTTTTfilter(seq)
    print("Step 4/8 - 5'-G filter and TTTT filter COMPLETE.")

    GseqNGG = NGGfilter(Gseq)
    print("Step 5/8 - NGG filter COMPLETE.")

    targetPosition = distanceBetweenTargets(GseqNGG)
    print("Step 6/8 - compute deletion size COMPLETE.")

    deletionSize = delSizeFilter(targetPosition)
    print("Step 7/8 - deletion size filter COMPLETE.")

    fullyFiltered = segDupFilter(deletionSize)
    print("Step 8/8 - LCR filter COMPLETE.")
	
    # Returns a CSV containing (candidate gRNA, target position 1, target position 2, distance between targets) for each row
    printCSV(fullyFiltered, 2, "_SCORE")


def CRISPRdesign(file_name, guideLength = 20, minRange = 0, maxRange = -1,
                 fiveGfilter = True):
    # CRISPR design with 20-bp target with a 5'-G and ends with 5'-NGG-3' (pX330 or pX459 design)
    # for generating a single double-stranded cut
    # file_name: fasta file with a single DNA sequence
    # minRange/maxRange = lower and upper boundaries for the the sequence to be searched
    # If maxRange = -1, it means to extract from minRange to end of sequence
    # fiveGfilter = True means filter for 5'-G in the sequence of (pX330 or pX459 design)
    DNAseq = read_fasta(file_name).seq
    if maxRange == -1:
        DNAseq = DNAseq[minRange:]
    else:
        DNAseq = DNAseq[minRange:maxRange]

    DNAseqRC = DNAseq.reverse_complement()
    DNAseq = str(DNAseq)
    DNAseqRC = str(DNAseqRC)


    def seqListOnly(kmer):
        # Converts output from kmer_count and nplicate into a single list of kmers
        ret_list = []
        for oligo in kmer:
            ret_list.append(oligo[0])

        return ret_list


    def fivePrimeGfilter(oligo):
        # Filters DNA targets for 5'-G
        # Takes in a list of sequences (processed by seqListOnly method) as argument
        ret_list = []

        for target in oligo:
            if target[0] == "G":
                ret_list.append(target)

        # Returns a list of sequences with a leading 5'-G base
        return ret_list


    def findPosition(oligo):
        # Given a sequence, finds the position on DNAseq and DNAseqRC
        # + position = target on (+)-sense DNA
        # - position = target on (-)-sense DNA
        pos = DNAseq.find(oligo)
        neg = DNAseqRC.find(oligo)

        if pos == -1:
            return -1 * neg
        else:
            return pos


    def NGGfilter(oligo):
        # Filters DNA targets for NGG PAM downstream
        # Takes in a list of "sequence" as argument
        ret_list = []
        init, end = guideLength + 1, guideLength + 3

        for seq in oligo:
            pos1 = findPosition(seq)
            if pos1 > 0:
                if DNAseq[pos1 + init : pos1 + end] == "GG":
                    ret_list.append((seq, pos1))
            else:
                if DNAseqRC[-pos1 + init : -pos1 + end] == "GG":
                    ret_list.append((seq, pos1))
                
        # Returns a list of tuple (sequence, position) for the presence of NGG PAM in targets
        return ret_list


    kmer = kmer_count(file_name, guideLength, minRange = minRange, maxRange = maxRange)
    print("Step 1/5 - kmer counting COMPLETE.")

    onePlicate = nplicate(kmer, 1)
    print("Step 2/5 - filter for unique targets COMPLETE.")

    seq = seqListOnly(onePlicate)
    print("Step 3/5 - clean data COMPLETE.")

    if fiveGfilter:
        Gseq = fivePrimeGfilter(seq)
        print("Step 4/5 - 5'-G filter COMPLETE.")
    else:
        Gseq = seq
        print("Step 4/5 - No 5'-G filter applied COMPLETE.")

    GseqNGG = NGGfilter(Gseq)
    print("Step 5/5 - NGG filter COMPLETE.")

    # Returns a CSV containing (candidate gRNA, target position) for each row
    printCSV(GseqNGG, 1, "_min" + str(minRange) + "_max" + str(maxRange) + "_CRISPR")    


def findOverlap(file_name1, file_name2):
    # Searches for targets found in both files and returns them
    # Useful for combinding scoreGuideDesign or CRISPRdesign results from two different assemblies (e.g. hg19 and hg38)
    # Takes the scoreGuideDesign output as arguments 
    csv1 = open(file_name1, 'r')
    csv2 = open(file_name2, 'r')

    csv1 = csv1.read()
    csv2 = csv2.read()

    csv1 = csv1.split('\n')
    csv2 = csv2.split('\n')

    combined = csv1 + csv2

    seq_list = []
    
    for row in combined:
        if row != "":
            seq_list.append(row.split(',')[0])

    seq_list = Counter(seq_list)

    # Write CSV file
    csvwrite = open("combined.csv", "wb")
    csv_writer = csv.writer(csvwrite, delimiter = ',')
    for seq in seq_list:
        if seq_list[seq] > 1:
            csv_writer.writerow((seq, ))

    csvwrite.close()

    print("combined.csv PRINTED.")
            


                            
""" THE CODE BELOW DOES COUNT UP ALL KMERS CORRECTLY! IT FILTERS FOR CRISPR
    TARGETS AND OVERLOOK IDENTICAL KMERS THAT ARE NOT VALID CRISPR TARGETS!!!
def CRISPRkmer_count(file_name, k):
    # CRISPR designe with k-bp target with a 5'-G and ends with 5'-NGG-3'
    # Take in a fasta file with one DNA sequence and a value for length of k-mer
    
    # Read the sequence along with its reverse complement
    DNAseq = read_fasta(file_name).seq
    # (-)-sense DNA
    DNAseqRC = DNAseq.reverse_complement()

    # Count kmers
    store = {}
    for i in range(len(DNAseq) - k -2):
        # Subset the sequence for kmer and keep count as value in store
        PAM_init = i + k + 1
        PAM_end = i + k + 3
        
        if DNAseq[i] == "G" and DNAseq[PAM_init: PAM_end] == "GG":
            oligo1 = str(DNAseq[i : i + k])
            store[oligo1] = store.get(oligo1, 0) + 1

        if DNAseqRC[i] == "G" and DNAseqRC[PAM_init: PAM_end] == "GG":
            oligo2 = str(DNAseqRC[i : i + k])
            store[oligo2] = store.get(oligo2, 0) + 1

    # Return a list of tuples ('sequence', count)        
    return store.items()
"""    
