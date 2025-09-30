The code is for personal use only - not for profit

**MAIN-CODE_FILE = y-dna-analyzer.py**

This code is still a work in progress, and I expect to update it from time to time.

**A Windows-ready exe file is now available here: https://github.com/sgbmzm/y-dna-analyzer/releases/download/y-dna-analyzer/y-dna-analyzer_30_Sep_25.exe**
The file is not signed, because signing costs money, but there really are no viruses in it.

or run the python code (It is recommended through the software https://thonny.org/) 


Use the latest Reference file from here
https://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz
This file is updated every day with new variants

You also need the Reference file "Msnps_hg19.vcf.gz" that I made from the file https://ybrowse.org/gbrowse2/gff/snps_hg19.csv (because there is a problem of mixing genomic locations between the two references in the file https://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz)

Both files should be in the folder where the code is saved.

This code analyzes which branch of the Y chromosome tree the sample is on.
It extracts Y chromosome data from raw files from FTDNA, MyHeritage, Ancestry, 23&Me, and from whole genome (or only-y) VCF files (including ftdna-big-y-vcf).
**There is no need to extract the raw files. The current software can process them both when compressed and when extracted.**


Some files contain only a small amount of Y chromosome information, so the result is a high branch that is not up-to-date.
Usually, VCF files will yield a up-to-date branch

This code was born only because of https://pypi.org/project/yclade/
Thank them for that.

This code was created in collaboration with GPT-Chat. It stopped a lot for me.

<img width="802" height="608" alt="image" src="https://github.com/user-attachments/assets/3c646631-ac8e-4c47-9c9a-af65f5ab7cbd" />

