CNV_engineering
------
Copy number variations project. Engineering 22q11.21 deletion clones using CRISPR.

**GOAL:** Identify all possible guides that meet the conditions for CNV engineering using CRISPR/Cas9 system

**REQUIRED FILES:**
- FASTA sequence file(s) containing target of interest
- kmerDNA.py for computation (OBSOLETE)
- kmerDNA2.py for computation

**USAGE:** Run kmerDNA.py to find all possible CRISPR guides for a given sequence, with the option of applying SCORE method (Tai et al. *Nature*. 2016) to generate microdeletions/duplications by targeting segmental duplications

*NOTE:* Method for dual guide design is outdated (no fuctionality to eliminate sequences with "TTTT" which interferes with transcription, and to screen for cloning incompatibility - BbsI cut sites)

DETAILED INSTRUCTIONS:
------
1. Download the most current version of Python 2.7 at python.org/download.
2. Extract a single DNA sequence spanning the entire CNV region and save as a FASTA file. I do this in UCSC Genome Browser.
3. Download and save the kmerDNA.py file in the same directory as FASTA file.
4. Open the python development environment "IDLE (Python 2.7 GUI)"
5. Open the kmerDNA.py file in IDLE ("File" > "Open..." and select kmerDNA.py)
6. In the .py file, go to "Run" > "Run Module" or press F5.
7. Call the scoreGuideDesign function along with the necessary parameters. For example: scoreGuideDesign("22q11.21 hg19 UCSC.fasta", minSize = 2330000, maxSize = 3500000, minDupA = 560000, maxDupA = 1130000, minDupB = 3460000, maxDupB = 4020000, guideLength = 20).
8. min/maxSize = min/max distance between the two target sites for each guide
9. min/maxDupA = indices denoting the endpoints of LCR 1
10. min/maxDupB = indices denoting the endpoints of LCR 2
11. The output CSV file should contain the guide sequences, indices of guide target 1 (positive = 5' to 3' sense; negative = 5' to 3' antisense strand), indices of guide target 2, and distances between the targets.

Highly recommend extracting the hg19 and hg38 sequences for the CNV region of interest and running the scoreGuideDesign function on both. You can then use the findOverlap function to identify guides present in both files, ultimately ones that fit both assemblies. Command: findOverlap("CSVfile1", "CSVfile2").
