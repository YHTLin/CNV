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
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW
from collections import Counter
from collections import defaultdict
import csv
import sys   # For terminating functions after exception is thrown
      

def scoreGuideDesign(file_name,
                     minSize = 1970000, maxSize = 4000000,
                     minDupA = 580000, maxDupA = 1640000,
                     minDupB = 3610000, maxDupB = 4170000,
                     guideLength = 20, GandTTTTfilter = True):
    # CRISPR design with 20-bp target with a 5'-G and ends with 5'-NGG-3' (pX330 or pX459 design)
    #   for generating microdeletions using SCORE approach (Tai et al. Nature neuroscience 2016)
    # file_name: fasta file with a single DNA sequence
    # minSize: Acceptance limit for minimum distance between targets (defaults for 22q11.21 hg38 UCSC.fasta)
    # maxSize: Acceptance limit for maximum distance between targets (defaults for 22q11.21 hg38 UCSC.fasta)
    # min/max-DupA/DupB = lower and upper ranges for duplication region A and B (defaults for 22q11.21DS)
    # GandTTTTfilter = True to run 5'-G and TTTT filter and False otherwise
    

    def read_fasta(file_name):
        # Read a fasta file (can contain only one sequence) and return Bio.SeqRecord object
        return SeqIO.read(file_name, "fasta", alphabet = IUPAC.ambiguous_dna)


    def kmer_count(file_name, kmer):
        # Take in a fasta file with one DNA sequence and a value for length of k-mer

        # Read the sequence along with its reverse complement
        DNAseq = read_fasta(file_name).seq
        DNAseqRC = DNAseq.reverse_complement()

        # Count kmers
        store = defaultdict(list)
        for k in range(len(DNAseq) - kmer + 1):
            # Subset the sequence for kmer and keep count as value in store
            oligo1 = str(DNAseq[k : k + kmer])
            store[oligo1].append(k + 1)

            oligo2 = str(DNAseqRC[k : k + kmer])
            store[oligo2].append(-(k + 1))

        # Return a list of tuples ('sequence', [coordinate#1, coordinate#2...])
        # coordinates are 1-based indices, negative = on the reverse complement strand
        return store.items()
        

    def nplicate(kmer_list, Nplicate_filter = 2):
        # Takes in return value from kmer_count method and filters for duplicates by default
        # Optional argument: Nplicate_filter takes a positive integer to filter the kmer by frequency of occurrence
        ret_list = []
        
        for kmer in kmer_list:
            if len(kmer[1]) == Nplicate_filter:
                ret_list.append(kmer)

        # Returns a list of tuples('sequence', [coordinate#1, coordinate#2...])
        return ret_list


    def fivePrimeGandTTTTfilter(kmer_list):
        # Filters DNA targets for 5'-G and eliminate guides containing "TTTT" for U6 promoter transcriptional compatibility
        # Takes in a list of sequences (processed by seqListOnly method) as argument
        ret_list = []

        for kmer in kmer_list:
            if kmer[0][0].upper() == "G" and "TTTT" not in kmer[0].upper():
                ret_list.append(kmer)

        # Returns a list of tuples('sequence', [coordinate#1, coordinate#2...]) with a leading 5'-G base and no "TTTT" in sequence
        return ret_list


    def distanceBetweenTargets(kmer_list):
        # Lays out all the positions of the targets and their distance apart
        ret_list = []
        init, end = guideLength - 1, guideLength + 2

        for kmer in kmer_list:
            pos1, pos2 = kmer[1]
            if pos1 > 0:
                if pos2 > 0:
                    ret_list.append((kmer[0], DNAseq[pos1 + init:pos1 + end], DNAseq[pos2 + init:pos2 + end],
                                     pos1, pos2, abs(pos1 - pos2)))
                else:
                    ret_list.append((kmer[0], DNAseq[pos1 + init:pos1 + end], DNAseqRC[-pos2 + init:-pos2 + end],
                                     pos1, pos2, abs(pos1 - (len(DNAseq) + pos2 + 1))))
            else:
                if pos2 > 0:
                    ret_list.append((kmer[0], DNAseqRC[-pos1 + init:-pos1 + end], DNAseq[pos2 + init:pos2 + end],
                                     pos1, pos2, abs(pos2 - (len(DNAseq) + pos1 + 1))))
                else:
                    ret_list.append((kmer[0], DNAseqRC[-pos1 + init:-pos1 + end], DNAseqRC[-pos2 + init:-pos2 + end],
                                     pos1, pos2, abs(pos1 - pos2)))

        # Returns a list of tuple ("sequence", NGG1, NGG2, target position1, target position2, distance between targets)
        return ret_list


    def delSizeFilter(dist_data):
        # Checks that size of the deletions is within set limits (minSize/maxSize)
        ret_list = []

        for data in dist_data:
            if minSize <= data[5] <= maxSize:
                ret_list.append(data)

        return ret_list


    def segDupFilter(dist_data):
        # Checks that the targets are within the provided segmental duplication regions (min/maxDupA/B)
        ret_list = []

        for data in dist_data:
            if data[3] > 0:
                pos1 = data[3] - 1   # Adjust to 1-indexing
            else:
                pos1 = len(DNAseq) + data[3]   # Adjust to 1-indexing
            if data[4] > 0:
                pos2 = data[4] - 1   # Adjust to 1-indexing
            else:
                pos2 = len(DNAseq) + data[4]   # Adjust to 1-indexing
                
            if (minDupA <= pos1 <= maxDupA and minDupB <= pos2 <= maxDupB) or (minDupA <= pos2 <= maxDupA and minDupB <= pos1 <= maxDupB):
                ret_list.append(data)
                
        # Returns a list of tuple ("sequence", target position1, target position2, distance between targets)
        return ret_list


    def NGGfilter(filtered_data):
        # Filters DNA targets for NGG PAM downstream
        # Takes in a list of "sequence" as argument
        ret_list = []
        init, end = guideLength, guideLength + 2

        for seq in filtered_data:
            pos1, pos2 = seq[3], seq[4]
            
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


    def printCSV(kmer_list, Nplicate = "N", dscrpt = ""):
        # Takes in return value from kmer_count method
        # Nplicate: to specify whether nplicate method was applied
        # dscrpt: for a more descriptive title

        try:
            kmer = len(kmer_list[0][0])
        except IndexError:
            print("ATTENTION: No SCORE guides identified.")
            sys.exit(1)
            
        file_name = str(kmer) + "mer_" + str(Nplicate) + "plicate" + dscrpt + ".csv"

        # Write CSV file
        csvwrite = open(file_name, "wb")
        csv_writer = csv.writer(csvwrite, delimiter = ',')
        csv_writer.writerow(["Sequence", "NGG 1", "NGG 2", "Target 1", "Target 2", "Distance between targets"])

        for row in kmer_list:
            if not "N" in row[0]:
                # Conditional removes any N bases from the output
                csv_writer.writerow(row)

        csvwrite.close()

        print(file_name + " PRINTED.")

    
    DNAseq = read_fasta(file_name).seq
    DNAseqRC = DNAseq.reverse_complement()
    DNAseq = str(DNAseq)
    DNAseqRC = str(DNAseqRC)

    print("Step 1/7 - kmer counting...")
    kmer = kmer_count(file_name, guideLength)
    print("Step 1/7 - COMPLETE.")
    
    print("Step 2/7 - filter for duplicates...")
    twoPlicate = nplicate(kmer, 2)
    print("Step 2/7 - COMPLETE.")

    if GandTTTTfilter:
        print("Step 3/7 - 5'-G filter and TTTT filter...")
        twoPlicate = fivePrimeGandTTTTfilter(twoPlicate)
        print("Step 3/7 - COMPLETE.")
    else:
        print("Step 3/7 - 5'-G filter and TTTT filter skipped.")

    print("Step 4/7 - compute deletion size...")
    targetPosition = distanceBetweenTargets(twoPlicate)
    print("Step 4/7 - COMPLETE.")

    print("Step 5/7 - deletion size filter...")
    deletionSize = delSizeFilter(targetPosition)
    print("Step 5/7 - COMPLETE.")

    print("Step 6/7 - LCR filter...")
    LCRfiltered = segDupFilter(deletionSize)
    print("Step 6/7 - COMPLETE.")

    print("Step 7/7 - NGG filter...")
    fullyFiltered = NGGfilter(LCRfiltered)
    print("Step 7/7 - COMPLETE.")
	
    # Returns a CSV containing (candidate gRNA, target position 1, target position 2, distance between targets) for each row
    printCSV(fullyFiltered, 2, "_SCORE")



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
            
