import csv
import gzip
import os

# נתיב קובץ ה-CSV
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
csv_file = os.path.join(desktop, "snps_hg19.csv")
vcf_file = os.path.join(desktop, "Msnps_hg19.vcf.gz")

# כותרות סטנדרטיות ל-VCF
vcf_header = [
    "##fileformat=VCFv4.2",
    "##fileDate=20171217",
    "##source=Edited from https://ybrowse.org/gbrowse2/gff/snps_hg19.csv",
    "##reference=hg19",
    "##INFO=<ID=.,Number=.,Type=String,Description=\"No info\">",
    "##FORMAT=<ID=.,Number=.,Type=String,Description=\"No format\">",
]

# שורה עם שמות העמודות
columns = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]

with gzip.open(vcf_file, 'wt') as vcf_out:
    # כתיבת הכותרות
    for line in vcf_header:
        vcf_out.write(line + "\n")
    
    # כתיבת שורת הכותרות
    vcf_out.write("\t".join(columns) + "\n")
    
    # קריאת ה-CSV והמרה לשורות VCF
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            chrom = row['CHROM'].replace("hg19ChrY", "ChrY")
            pos = row['POS']
            snp_id = row['ID']
            ref = row['REF']
            alt = row['ALT']
            
            # השאר העמודות ריקות
            qual = "."
            filter_col = "."
            info = "."
            fmt = "."
            
            vcf_line = [chrom, pos, snp_id, ref, alt, qual, filter_col, info, fmt]
            vcf_out.write("\t".join(vcf_line) + "\n")

print(f"VCF.GZ file saved as: {vcf_file}")
