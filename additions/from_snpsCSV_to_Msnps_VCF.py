# באמצעות הקוד הזה ממירים קובץ רפרנס מ CSV ל VCF.GZ
# https://ybrowse.org/gbrowse2/gff/snps_hg38.csv
#https://ybrowse.org/gbrowse2/gff/snps_hg19.csv

# הקוד הזה צריך להיות שמור באותה תיקייה שבה שומרים את קובץ הרפרנס שאותו רוצים להמיר.
# את קובץ HG19 אין מה להמיר שוב ושוב כי הוא לא מתעדכן כבר שנים רבות.
# את קובץ HG38 כדאי לעדכן מידי פעם כי הוא מתעדכן על בסיס יום יומי בווריאנטים חדשים


import csv
import gzip
import os
from datetime import datetime

def snps_to_Msnps(csv_path):
    """
    המרת CSV ל-VCF לפי אינדקסים של העמודות הרצויות במילון col_indices
    """
    
    # באיזה עמודות בקובץ CSV נמצאים הטורים שאני צריך
    # הנתונים הבאים מתאימים לקבצים של תומס קראן
    col_indices = {
        'CHROM': 0,
        'POS': 3,
        'ID': 8,
        'REF': 10,
        'ALT': 11
    }


    # === תאריך נוכחי בפורמט יפה ===
    today = datetime.now().strftime("%Y-%m-%d")

    # === זיהוי reference לפי שם הקובץ ===
    csv_name = os.path.basename(csv_path).lower()
    if "hg19" in csv_name:
        reference = "hg19"
    elif "hg38" in csv_name:
        reference = "hg38"
    else:
        reference = "unknown"

    # === יצירת שם אוטומטי לקובץ הפלט ===
    vcf_path = f"Msnps_{reference}.vcf.gz"

    # === בניית כותרת VCF ===
    vcf_header = [
        "##fileformat=VCFv4.2",
        f"##fileDate={today}",
        "##source=Converted from CSV",
        f"##reference={reference}",
        "##INFO=<ID=.,Number=.,Type=String,Description=\"No info\">",
        "##FORMAT=<ID=.,Number=.,Type=String,Description=\"No format\">",
    ]

    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    # === כתיבת קובץ VCF.GZ ===
    with gzip.open(vcf_path, 'wt') as vcf_out:
        for line in vcf_header:
            vcf_out.write(line + "\n")
        vcf_out.write("\t".join(columns) + "\n")

        with open(csv_path, newline='') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)  # דילוג על שורת כותרת
            for row in reader:
                chrom = row[col_indices['CHROM']].replace("hg19ChrY", "ChrY")
                pos = row[col_indices['POS']]
                snp_id = row[col_indices['ID']]
                ref = row[col_indices['REF']]
                alt = row[col_indices['ALT']]

                qual = "."
                filter_col = "."
                info = "."
                fmt = "."

                vcf_line = [chrom, pos, snp_id, ref, alt, qual, filter_col, info, fmt]
                vcf_out.write("\t".join(vcf_line) + "\n")

    print(f"VCF.GZ file saved as: {vcf_path}")


# === דוגמה לשימוש ===

# אם הקובץ שמור באותה תיקייה שבה שמור הסקריפט הזה
snps38_file = "snps_hg38.csv"
snps19_file = "snps_hg19.csv"

# קריאה לפונקצייה
snps_to_Msnps(snps38_file)
#snps_to_Msnps(snps19_file)



