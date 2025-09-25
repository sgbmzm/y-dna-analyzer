import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import csv
import os
import re
import gzip
import zipfile
import io
import webbrowser
import yclade
from yclade import tree, snps, find
import networkx as nx

# טעינת עץ ווייפול
yfull_tree_data = tree.get_yfull_tree_data(version=None, data_dir=None)


# ------------------------
# משתנים גלובליים
# ------------------------
last_clades = []
last_positive_snp_string = ""
last_ref_type = ""
reference_snps = {}      # pos -> ref_snp_name + ref_allele
reference_names = {}      # id_snp_name -> pos
user_snaps = {}          # pos -> ref_snp_name + ref_allele
last_dna_file = ""
last_reference_file = ""
reference_loaded = False
user_loaded = False
last_dna_file_type = ""
last_yfull_link = ""
last_ftdna_link = ""
last_ab_data = None
last_dna_file_info = ""


# פונקצייה לקבל מיקום מוחלט של קובץ שאמור להיות כלול בתוכנה
# ראו: https://stackoverflow.com/questions/7674790/bundling-data-files-with-pyinstaller-onefile/44352931#44352931
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# ------------------------
# איפוס משתמש
# ------------------------
def reset_user():
    global last_clades, last_positive_snp_string, last_dna_file, user_loaded 
    
    last_clades = []
    last_positive_snp_string = ""
    last_dna_file = ""
    user_loaded = False

    result_var.set("")
    btn_yfull.grid_forget()
    btn_ftdna.grid_forget()
    dna_loading_label.config(text="No DNA_file loaded", fg="red")
    yclade_label.config(text="load user DNA file", fg="red")
    user_result_var.set("")
    user_result_label.config(text="", bg="SystemButtonFace")
    btn_save_results.grid_forget()
    btn_unload_dna.grid_forget()
    last_dna_file_type = ""
    last_dna_file_info = ""
    last_yfull_link = ""
    last_ftdna_link = ""

    
    
'''
# זה לא עבד טוב וזה לשמירה בלבד
הפוקצייה הבאה ברוך השם עובדת
# פונקצייה מאוד חשובה שמחזירה את שמות כל ענפי הצאצאים שתחת ענף מסויים בוויפול
def get_branch_and_descendants(snp_name: str, include_descendants: bool = True):
    """מחזיר רשימת ענפים לפי SNP, כולל צאצאים אם include_descendants=True.
       יוצר את tree_data בתוך הפונקציה.
    """
    # טען את העץ
    tree_data = tree.get_yfull_tree_data(version=None, data_dir=None)
    
    snp_name = snps.SnpResults(positive={snp_name.rstrip("+-")}, negative=set())
    
    snp_name = snps.normalize_snp_results(snp_results=snp_name,snp_aliases=tree_data.snp_aliases)

    # מציאת הצומת לפי ה-SNP
    node = next(
        (n for n, snps_set in tree_data.clade_snps.items() if snp_name.rstrip("+-") in snps_set),
        None
    )
    if node is None:
        raise ValueError(f"לא נמצא צומת מתאים ל-{snp_name}")

    # קביעת הצמתים לכלול
    nodes = [node] + list(nx.descendants(tree_data.graph, node)) if include_descendants else [node]

    # מיון מהשורש לעלים
    nodes_sorted = sorted(nodes, key=lambda n: len(list(nx.ancestors(tree_data.graph, n))))
    
    # זה לשימוש עתידי אם רוצים שיאסוף גם את רשימת כל הווריאנטים שתחת ענף מסויים
    #snps_list = [tree_data.clade_snps.get(n, set()) for n in nodes_sorted]
    
    return nodes_sorted
'''

# פונקצייה מאוד חשובה שמחזירה את שמות כל ענפי הצאצאים שתחת ענף מסויים בוויפול
def get_clade_and_descendants_lists(tree_data, snp_name: str, include_descendants: bool = True, merge_snps: bool = False):
    
    """מחזיר שתי רשימות מקבילות:
       - כל שמות הקליידים (או רק הענף עצמו)
       - כל ה-SNPים של אותם קליידים (או רשימה אחת מאוחדת אם merge_snps=True)
       בסדר מהשורש לעלים.
       
       פרמטרים:
           include_descendants: True - מחזיר את הענף וכל הצאצאים, False - רק הענף
           merge_snps: True - מחזיר את כל ה-SNPים של הענף והצאצאים כרשימה אחת מאוחדת
    """
    
    # נירמול ה-SNP לפי אליאסים
    snp_results = snps.SnpResults(positive={snp_name.rstrip("+-")}, negative=set())
    snp_results_normalized = snps.normalize_snp_results(
        snp_results=snp_results,
        snp_aliases=tree_data.snp_aliases
    )
    if not snp_results_normalized.positive:
        raise ValueError(f"SNP {snp_name} לא נמצא בעץ לאחר נירמול")
    canonical_snp = next(iter(snp_results_normalized.positive))

    # מציאת הצומת הקנוני בעץ
    node = None
    for n, snps_set in tree_data.clade_snps.items():
        if canonical_snp in snps_set:
            node = n
            break
    if node is None:
        raise ValueError(f"לא נמצא צומת מתאים ל-{snp_name}")

    # קביעת הצמתים לכלול
    if include_descendants:
        descendants = nx.descendants(tree_data.graph, node)
        #ancestors = nx.ancestors(tree_data.graph, node) # אם רוצים להחזיר את האבות ולא את הצאצאים
        nodes = [node] + list(descendants)
    else:
        nodes = [node]

    # מיון לפי סדר מהשורש לעץ
    def distance_from_root(n):
        return len(list(nx.ancestors(tree_data.graph, n)))

    nodes_sorted = sorted(nodes, key=distance_from_root)

    # יצירת הרשימות
    branch_names = nodes_sorted
    snps_list = [tree_data.clade_snps.get(n, set()) for n in nodes_sorted]

    if merge_snps:
        snps_list = [set().union(*snps_list)]  # מאחד את כל ה-SNPים לרשימה אחת

    return branch_names#, snps_list



    

# פונקצייה לפתיחת קבצים שונים כל קובץ לפי סוג הפתיחה הדרוש
def universal_opener(file_path, only_dna_if_zip=False):
    # אם זה קובץ דחוס בשיטת gz
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt", encoding="utf-8")
    
    # אם זה קובץ דחוס בשיטת zip
    elif file_path.endswith(".zip"):
        zf = zipfile.ZipFile(file_path, "r")
        # מניחים שמעוניינים דווקא בקובץ הראשון שיש בזיפ
        first_file_in_zip = zf.namelist()[0]    
        ####################################################
        # בדיקה אם רוצים דווקא קובץ מייהירטייג
        if only_dna_if_zip:
            lower_file_name = first_file_in_zip.lower()
            # תנאי: שיהיה גם במחרוזת השם "myheritage" וגם שהסיומת תתאים לקובץ של מייהירטייג
            # וגם שהסיומת לא תהיה vcf כי אם כן זה קובץ של ftdna של תוצאות big-y
            if ("dna" not in lower_file_name and not lower_file_name.endswith(".csv") and not lower_file_name.endswith(".vcf")):
                messagebox.showerror("Error", f"Not Fund 'dna' in first filename in zip: {first_file_in_zip}")
                zf.close()
                return None
        ####################################################
        # אם לא רוצים לוודא מייהירטייג אז מחזירים את הקובץ הראשון בזיפ
        inner_file = zf.open(first_file_in_zip)
        return io.TextIOWrapper(inner_file, encoding="utf-8")
    
    # אם זה קובץ רגיל
    else:
        return open(file_path, "rt", encoding="utf-8")
    
    
# פונקצייה שמחזירה את השם של הקובץ הראשון בתוך זיפ שבדרך כלל הוא העיקרי
def get_first_file_name_in_zip(zip_path):
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            return zf.namelist()[0]  # מחזיר את השם של הקובץ הראשון
    except:
        return None


# זה מיוחד לקבצי y של ftdna למי שעושה בדיקת פמילי פינדר שלא כתוב בהם איזה רפרנס הם אבל הם hg38
# פונקצייה שבודקת האם מסתיים רק ב gz ולא לדוגמא vcf.gz
def is_gz_only(filename):
    parts = filename.lower().split(".")
    return parts[-1] == "gz" and len(parts) == 2  # רק סיומת אחת לפני ה-gz

# פונקצייה שמנסה למצוא שם של הקובץ הפנימי שדחוס ב gz
def get_gz_internal_filename(filepath):
    try:
        with open(filepath, 'rb') as f:
            header = f.read(10)          # קופצים על ההדר
            if header[2] != 8 or not header[3] & 0x08:  # בדיקה קצרה: קומפרסיה = deflate ו-FNAME
                return None
            name = b''
            while (b := f.read(1)) not in (b'', b'\x00'):
                name += b
            return name.decode('latin-1')
    except:
        return None

##fileformat=MyHeritage
##format=MHv1.0
##chip=GSA
##timestamp=2022-03-03 06:27:05 UTC
##reference=build37
# DOWNLOADED DATA WILL NO LONGER BE PROTECTED BY OUR SECURITY MEASURES.
#RSID CHROMOSOME POSITION RESULT # אין סולמית במקור
#rs547237130 1 72526 AA # אין סולמית במקור
# כרומוזום Y נקרא בקובץ:   Y


#AncestryDNA raw data download
#This file was generated by AncestryDNA at: 01/23/2021 17:05:38 UTC
#Data was collected using AncestryDNA array version: V2.0
#Data is formatted using AncestryDNA converter version: V1.0
#Below is a text version of your DNA file from Ancestry.com DNA, LLC.  THIS
#Genetic data is provided below as five TAB delimited columns.  Each line 
#corresponds to a SNP.  Column one provides the SNP identifier (rsID where 
#possible).  Columns two and three contain the chromosome and basepair position 
#of the SNP using human reference build 37.1 coordinates.  Columns four and five
#rsid chromosome position allele1 allele2 # אין סולמית במקור
#rs3131972 1 752721 G G # אין סולמית במקור
# כרומוזום Y נקרא בקובץ:   24


# This data file generated by 23andMe at: Fri Jul 31 00:41:13 2020
# This file contains raw genotype data, including data that is not used in 23andMe reports.
# This data has undergone a general quality review however only a subset of markers have been 
# individually validated for accuracy. As such, this data is suitable only for research, 
# educational, and informational use and not for medical or other use.
# Below is a text version of your data.  Fields are TAB-separated
# Each line corresponds to a single SNP.  For each SNP, we provide its identifier 
# (an rsid or an internal id), its location on the reference human genome, and the 
# genotype call oriented with respect to the plus strand on the human reference sequence.
# We are using reference human assembly build 37 (also known as Annotation Release 104).
# Note that it is possible that data downloaded at different times may be different due to ongoing 
# improvements in our ability to call genotypes. More information about these changes can be found at:
# https://you.23andme.com/p/24128f3ee36ceb87/tools/data/download/
# More information on reference human assembly builds:
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/
# rsid chromosome position genotype # יש סולמית במקור
#rs548049170 1 69869 TT # אין סולמית במקור
# כרומוזום Y נקרא בקובץ:   Y

# מייהירטייג: כתוב שם החברה, כתוב רפרנס 


# פונקצייה לבדיקה האם כתוב בקובץ באיזה רפרנס הוא משתמש באיזה תאריך נוצר ומי היוצר
def detect_headlines(file_path):
    
    ref = "unknown"
    creator = "unknown"
    creation_date = "unknown"
       
    #######################################################################
    if "haplocaller" in file_path:
        ref = "hg38"
        creator = "ftdna"
        creation_date = "unknown"
    
    if is_gz_only(file_path): # אם זה קובץ שמסתיים רק ב gz ולא vcf.gz וכדומה
        internal_name = get_gz_internal_filename(file_path) # מנסה למצוא האם יש שם לקובץ הפנימי שבתוך הדחיסה
        # אם יש שם ובתוך השם יש את המילה haplocaller סימן שזה קובץ y של ftdna ואני יודע שהוא hg38
        if internal_name and "haplocaller" in internal_name:
            ref = "hg38"
            creator = "ftdna"
            creation_date = "unknown"
            #return {"ref": ref, "creator": creator, "creation_date": creation_date}
   ############################################################################
    
    ref_map = {
        "hg19": ["build 37", "grch37", "hg19"],
        "hg38": ["build 38", "grch38", "hg38", "hg 38"],
    }
    
    with universal_opener(file_path) as f:
        for line in f:
            if not line.startswith("#"):
                break  # יציאה – כבר עברנו את ה-header
            
            # מקטין את השורה לאותיות קטנות
            lower_line = line.lower()
            
            # בודק האם כתוב מה התאריך
            if "date=" in lower_line:
                creation_date = lower_line.split("date=", 1)[1].strip()
            
            # בודק האם מופיע שם החברה בכותרות הקובץ ואם לא מופיע מחזיר "unknown"
            for a_creator in ["myheritage", "23andme", "ancestry"]:
                if a_creator in lower_line:
                    creator = a_creator
            
            # בודק מה הרפרנס שכתוב בקובץ שבו עשו שימוש ליצירת הקובץ 
            for a_ref, keys in ref_map.items():
                if any(k in lower_line for k in keys):
                    ref = a_ref
                
    return {"ref": ref, "creator": creator, "creation_date": creation_date}


# ------------------------
# טעינת רפרנס
# ------------------------
def load_reference(ref_path="ask"):
    global reference_loaded, reference_snps, reference_names, last_reference_file, last_ref_type
    
    if ref_path != "ask":
        ref_path = resource_path(ref_path)
    
    if ref_path == "ask":
        
        messagebox.showwarning("select Reference file", f"Now Need to select Reference file. file can be downloaded from: https://ybrowse.org/gbrowse2/gff/")

        ref_path = filedialog.askopenfilename(
            title="Select Reference File",
            filetypes=[("VCF files", "*.vcf *.vcf.gz *.json"), ("All files", "*.*")]
        )
        
    if not ref_path:
        return
    
    reference_loading_label.config(text=f"Loading reference {os.path.basename(ref_path)} ...")
    root.update()
    
    last_ref_type = detect_headlines(ref_path)["ref"]

    reference_snps = {}
    reference_names = {}
    
    opener = gzip.open if ref_path.endswith(".gz") else open
    try:
        with opener(ref_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                try:
                    chrom = fields[0]
                    pos = int(fields[1])
                    snp_names = fields[2] # שים לב - יכול להכיל כמה שמות מופרדים בפסיק
                    ref = fields[3].upper()
                    alt = fields[4].upper()
                except Exception:
                    continue
                
                # מדלג על כל השורות שלא שייכות לכרומוזום Y
                chrom_norm = chrom.lower().replace("chr", "")
                if chrom_norm != "y":
                    continue
                
                # מוסיף את השורה למילון שמאונדקס לפי המיקום הגנומי
                reference_snps[pos] = {"pos": pos, "name": snp_names, "ref": ref, "alt": alt}
                
                # לפצל את השמות לפי פסיק ולשמור כל שם בנפרד כמפתח במילון נוסף שמאונדקס לפי שמות
                if snp_names != ".":
                    for snp_name in snp_names.split(","):
                        #reference_names[snp_name] = {"pos": pos, "name": snp_name, "ref": ref, "alt": alt}
                        reference_names[snp_name] = pos
                        
        reference_loaded = True
        reference_loading_label.config(fg="blue", text=f"Reference loaded: \nname: {os.path.basename(ref_path)}  \n Y-SNPs: {len(reference_snps)}  \ntype: {last_ref_type} \ndate: {detect_headlines(ref_path)['creation_date']}")
        last_reference_file = resource_path(ref_path)
        btn_unload_ref.grid(row=1, column=0, padx=5, pady=5)
        
        
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load reference: {e}")
        #reference_loading_label.config(text="")
        
        

# ------------------------
# ניתוח CSV
# ------------------------

        
def load_dna_file():
    
    # כל פעם כשבוחרים קובץ דנא חדש הכל מתאפס ומחושב מהתחלה
    reset_user()
    
    # בחירת קובץ הנא של המשתמש
    file_path = filedialog.askopenfilename(
        title="Select DNA file",
        filetypes=[("DNA files", "*.txt *.csv *.gz *.zip *.vcf"), ("All files", "*.*")]
    )
    if not file_path:
        return
    
    dna_file_info = detect_headlines(file_path)
    
    # בדיקת או בחירת הרפרנס המתאים לקובץ הדנא של המשתמש
    ref_auto_detect = dna_file_info["ref"]
    
    if ref_var.get() == "Aautodetect":
        ref_path = "Msnps_hg19.vcf.gz" if ref_auto_detect == "hg19" else "snps_hg38.vcf.gz" if ref_auto_detect == "hg38" else "ask"
    
        if not ref_auto_detect:
            choice = messagebox.askyesnocancel("Reference not Autodetected", "Autodetected Reference faild\nChoose hg19, hg38, or select file manually.\n\nYes = hg19, No = hg38, Cancel = choose manually")
            ref_path = "Msnps_hg19.vcf.gz" if choice else "snps_hg38.vcf.gz" if (choice == False) else "ask"
        
        # במקרה שקובץ הרפרנס לא נמצא בתיקיית הסקריפט הנוכחי יש לבחור קובץ באופן ידני
        if ref_path != "ask" and not resource_path(ref_path):
            messagebox.showwarning("Reference file missing", f"Reference file:     {ref_path}     missing \nplease select file manually")
            ref_path = "ask"       
    else:
        ref_path = "ask"
    
    global last_reference_file, reference_loaded
    
    if ref_path != "ask":
        ref_path = resource_path(ref_path)
        
    if not last_reference_file == ref_path:
        reference_loaded = False
    
    if not reference_loaded:    
        # טעינת הרפרנס הנבחר    
        load_reference(ref_path)
    
    # הצהרה על משתנים גלובליים הדרושים להלן
    global last_clades, last_positive_snp_string, last_dna_file, user_snps, user_loaded, last_dna_file_type, last_dna_file_info
    
    # אם לא בחרו קובץ רפרנס מאפסים הכל ולא ממשיכים
    if not reference_loaded:
        reset_user()
        return

    last_clades = []
    last_positive_snp_string = ""
    last_dna_file = file_path
    last_dna_file_type = ref_auto_detect
    last_dna_file_info = dna_file_info

    dna_loading_label.config(text=f"Loading {os.path.basename(file_path)} ...")
    root.update()

    positive_snps = []
    user_snps = {}
    
    ########################################################################################################
      
    # זה מטפל במקרה מיוחד של ביג Y של פטדנא שהוא בזיפ רגיל אבל הוא VCF
    is_ftdna_big_y_vcf = file_path.endswith(".zip") and get_first_file_name_in_zip(file_path) and get_first_file_name_in_zip(file_path).endswith(".vcf")
        
    ##########################################################################################################
    is_vcf_file = file_path.endswith(".vcf") or file_path.endswith(".vcf.gz") or is_ftdna_big_y_vcf
    
    '''
    # כבר טופל לעיל במשתנה is_ftdna_big_y_vcf ונשאר כאן רק לזיכרון
    if not is_vcf_file and file_path.endswith(".zip"):
        first_file_nane_in_zip = get_first_file_name_in_zip(file_path)
        if first_file_nane_in_zip and first_file_nane_in_zip.endswith(".vcf"):
            is_vcf_file = True        
    '''
    
    try:
        with universal_opener(file_path, only_dna_if_zip=True) as f:
            header_found = False
            for line in f:
                line = line.strip()
                if not line or line.startswith("##") or line.startswith("#"):
                    continue
                
                if is_vcf_file:
                    header_found = True
                
                # השורה הראשונה בלי # היא ה-header
                if not header_found:
                    # מזהים אוטומטית פסיקים או טאבים
                    delimiter = "\t" if "\t" in line else ","
                    header = [h.strip().strip('"') for h in line.split(delimiter)]
                    header_found = True
                    continue

                # מזהים את הדלימיטר בשורה הנוכחית
                delimiter = "\t" if "\t" in line else ","
                parts = [p.strip().strip('"') for p in line.split(delimiter)]

                if len(parts) < 4:
                    continue  # דילוג על שורות לא תקינות
                
                if is_vcf_file:
                    # לקובץ ויסיאף יש סדר מסויים ומוקפד לעמודות
                    chrom, pos_str, rsid, ref, alt = parts[:5]
                    format_fields = parts[8].split(":")              # לדוגמה ["GT","AD","DP","GQ","PL"]
                    sample_index = 0 # לפעמים יש כמה נבדקים שונים וכל אחד מהם בטור אחד אחרי השני כאן מדובר באחד בלבד
                    sample_info = parts[9 + sample_index].split(":") # זה מכיל את כל המידע על נתוני הנבדק הנוכחי וכדלהלן
                    sample_dict = dict(zip(format_fields, sample_info)) # הופך למילון: {"GT":"0/0", "AD":"13,0", ...}
                    # עכשיו אפשר לגשת בצורה בטוחה:
                    gt_str = sample_dict.get("GT") # לוקח רק את החלק של ה־GT (למשל "0/1" או "1/1"
                    ad_str = sample_dict.get("AD") # לוקח כמה קריאות לכל ווריאנט לדוגמא: 17,0 אומר 17 קריאות כמו הרפרנס ואפס קריאות כמו ה-אלט 
                    dp_str = sample_dict.get("DP") # כמה קריאות היו בסך הכל
                    pl_str = sample_dict.get("PL") # הסתברויות לכל קריאה
                    
                    '''
                    # הכנה לבדיקת איכות הקריאה לפי עומק הקריאה של החיובי לעומת השלילי
                    if ad_str:
                        ref_count, alt_count = map(int, ad_str.split(","))
                        # תנאים:
                        # 1. ALT לפחות פי 3 מ-REF
                        # 2. ALT לפחות 10 קריאות
                        if not alt_count >= 3 * ref_count and alt_count >= 5:
                            continue
                    '''
                    
                    alleles = [int(a) for a in re.split("[/|]", gt_str) if a.isdigit()]
                    if any(a > 0 for a in alleles): # כלומר אם יש שם משהו שהוא 1 אז יש משהו חיובי וזה אומר שהוא כמו ה alt
                        is_positive = True
                        alleles_str = alt
                    else:
                        is_positive = False
                        alleles_str = ref
                else: # זה לא ויסיאף אלא קובץ שרובו אוטוזומלי מהחברות הקטנות כמו מייהירטייג אנססטרי וכו ומכיל מעט מידע על כרומזום Y
                    # השדות הראשונים תמיד הם: rsid, chrom, pos, alleles
                    rsid, chrom, pos_str, alleles_str = parts[:4]
                    ad_str = ""
                
                # עושה אלל אחד גדול גם היה קטן וגם אם היו שניים יחד
                allele_str = alleles_str.upper().strip()[0] if alleles_str else ""
            
                # נרמול שם הכרומוזום ל־Y
                chrom = chrom.replace("chr", "").upper()
                if chrom not in ("Y", "24"):
                    continue
                
                # לא לכלול את המיקום הזה בקובץ של מייהירטייג כי הוא תמיד חיובי לכולם
                if dna_file_info["creator"] == "myheritage" and rsid == "rs570569843":
                    allele_str += "???" # או פשוט לדלג על השורה באמצעות: continue
                 
                # הוספת השורה למילון המשתמש לאחר שוודאנו שמובר בשורה של Y
                user_snps[int(pos_str)] = {"chrom": chrom, "pos_str": pos_str, "ref": ref_auto_detect, "snp_name": "?", "allele": allele_str, "is_positive": "?", "ad-R/A": ad_str}
                
                # בדיקה מול הרפרנס
                if not pos_str.isdigit():
                    continue
                ref_info = reference_snps.get(int(pos_str))
                if not ref_info:
                    continue

                if allele_str == ref_info["alt"]:
                    s = f"{ref_info['name']}+"
                    positive_snps.append(s)
                    user_snps[int(pos_str)]["is_positive"] = "Yes"   # או "כן" / "+" או כל מה שאתה רוצה
                    user_snps[int(pos_str)]["snp_name"] = ref_info['name']
                elif allele_str == ref_info["ref"] or allele_str == ".":
                    user_snps[int(pos_str)]["is_positive"] = "No"
                    user_snps[int(pos_str)]["snp_name"] = ref_info['name']

        last_positive_snp_string = ", ".join(positive_snps)
        #print("Final positive SNPs:", last_positive_snp_string)

        dna_loading_label.config(fg="blue", text=f"DNA file loaded: \nname: {os.path.basename(file_path)} \n{len(user_snps)} total Y-rows  \n{len(positive_snps)} Positive Y-SNPs in DNA_file  \ntype: {last_dna_file_type}")
        user_loaded = True
        
        if len(positive_snps) >= 100:
            run_calculate_clade()
        else:
            yclade_label.config(text="female / incorrect reference \n(Too little Y positive variants)", fg="red")
        
        # הצבת כפתור ביטול הטעינה    
        btn_unload_dna.grid(row=1, column=4, padx=5, pady=5)

    except Exception as e:
        messagebox.showerror("Error", f"Failed reading file: {e}")
        dna_loading_label.config(text="")


# פונקציה שטוענת מקובץ את הנתונים על קבוצות אבותינו וענפי הצאצאים שלהם ויכולה גם לשמור לקובץ
def get_ab_data():
    
    global last_ab_data, yfull_tree_data
    
    if last_ab_data: #  אם כבר הנתונים נטענו אין צורך לטעון אותם שוב מהקובץ
        return
    
    ab_groups_path = resource_path("ab_groups.csv")  # כאן את שם הקובץ שלך
    ab_data = []
    # טעינת כל שורה בקובץ לתוך מילון
    with open(ab_groups_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')  # אם הקובץ מופרד בפסיקים  
        for row in reader:
            # כל שורה היא מילון עם שמות העמודות כמפתחות
            ab_data.append(row)
    # בדיקה האם קיימים נתונים על ענפי הצאצים של כל שורה ואם לא קיימים אז להוסיף אותם               
    for row in ab_data:
        # אם אין מפתח 'sub_clades', ניצור אותו עם הערך "052"
        if 'sub_clades' not in row and row['Final SNP'] not in (None, '', 'None'):
            row['sub_clades'] = get_clade_and_descendants_lists(yfull_tree_data, row['Final SNP'])
    
    last_ab_data = ab_data
    
    write = False
    '''
    # זה כרגע לא טוב כי זה מפריע לקריאה של השדה row['sub_clades'] ###= ast.literal_eval(row['sub_clades']) 
    if write:
        # כתיבה חזרה (דורכת על הקובץ המקורי)
        fieldnames = list(ab_data[0].keys())  # שמות העמודות אחרי ההוספה
        with open(ab_groups_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')
            writer.writeheader()
            writer.writerows(ab_data)
    '''


# פונקצייה שמקבלת שם של ענף לדוגמא J2 או J-L243 ומחזירה את השם של הווריאנט שמגדיר אותו או את אחד מהם
def get_snp_from_clade(clade: str):
    if "-" not in clade: # לדוגמא J2 שלא בנוי J2-M172 ולכן אין לנו שם ווריאנט
        clade_snp = list(yfull_tree_data.clade_snps[clade])[0] # מחזיר אחד מהווריאנטים שמגדירים את הענף הנוכחי
    else:
        # הוצאת שם הווריאנט מתוך שם הענף דוגמא: J-L243 נהיה L243
        clade_snp = clade.split("-")[-1]
    return clade_snp
        
                      
# פונקצייה מאוד חשובה שמחזירה האם ענף מסויים נמצא בתוך קבוצת אבותינו או מעל אחת או יותר מקבוצות אבותינו
def get_ab_from_clade(clade: str, from_snp = False):    # הצהרה על משתנים גלובליים
    global last_ab_data, yfull_tree_data
    # אם רוצים מסניפ אז צריך קודם לדעת אל איזה ענף הסניפ הזה יושב בעץ ואז ממשיכים עם הענף המתאים
    if from_snp:
        clade = yclade.find_clade(clade_snp)
    # מוציאם משם הקלייד את הווריאנט שמגדיר אותו באמצעות פונקצייה
    clade_snp = get_snp_from_clade(clade)
    # חישוב כל הענפים שתחת הענף מחזיר רשימה
    clade_sub_clades = get_clade_and_descendants_lists(yfull_tree_data, clade_snp)
    
    # משתנים לצורך הספירה ומיכלים לאיסוף תוצאות
    clade_found_ab = False  
    sub_clades_found_ab = False  
    clade_found_ab_rows = []
    sub_clades_found_ab_rows = []

    # עוברים על כל שורה של קבוצות אבותינו כלומר על כל קבוצה ובודקים עבור כל ענף וענף בקבוצה את הדברים הבאים
    for row in last_ab_data:
        ab_clades_list = row.get('sub_clades', [])
        ab_clades_clean = [s.replace("*", "") for s in ab_clades_list]  # מסירים כוכביות

        # בדיקה האם הענף המבוקש נמצא בתוך אחת מקבוצות אבותינו
        if clade in ab_clades_clean:
            clade_found_ab = True
            clade_found_ab_rows.append(row)

        # בדיקה האם אחד ממתי הענפים של הענף הנבוקש נמצא בתוך אחת מקבוצות אבותינו
        if any(a_clade.replace("*", "") in ab_clades_clean for a_clade in clade_sub_clades):
            sub_clades_found_ab = True
            sub_clades_found_ab_rows.append(row)
    
    # A אומרת שיש בקבוצה הזו אשכנזים. זה אוסף לרשימה של כל קבוצות אבותינו שנמצאו
    # בתוך ab אמור לצאת רק אחד ולא יותר אבל לפעמים יוצא יותר בגלל שיש ענפי ab שהענף שלהם מוגדר רחב ולא ספציפי ברשימה שבקובץ
    in_ab = ", ".join([f"{r['AB-Group']}{'(A)' if 'A' in r['Communities'] else ''}" for r in clade_found_ab_rows])
    above_ab = ", ".join([f"{r['AB-Group']}{'(A)' if 'A' in r['Communities'] else ''}" for r in sub_clades_found_ab_rows])
    ab_string = f'in AB: \n[{in_ab}]' if in_ab else f'above ABs: \n[{above_ab}]' if above_ab else ""
    # מחזירים את רשימת ענפי אבותינו שבתוך או שמעל    
    return ab_string

    
               
def run_calculate_clade(Final_clade_index = 0):
    global last_positive_snp_string, last_clades, last_reference_file, ref_user_file, last_dna_file_type, last_ref_type, last_ab_data, last_dna_file_info
    result_var.set("") # תמיד לאפס קודם ולרוקן את הכיתוב הישן
    if not last_positive_snp_string:
        result_var.set("No matching SNPs were found in the file")
        yclade_label.config(text="Analysis finished - no matching SNPs found.")
        return

    try:
        yclade_label.config(text="Running clade calculation...")
        root.update()

        # קריאה ל־yclade עם מחרוזת אחת
        try:
            
            
            # clades = yclade.find_clade(last_positive_snp_string) זו הדרך הרגילה אבל היא ארוכה כי צריך כל פעם לטעון את העץ מחדש
            # לכן עושים את זה כך ישר על yfull_tree_data שכבר טעננו בתחילת הקוד פעם אחת
            snp_results = snps.parse_snp_results(last_positive_snp_string)
            snp_results = snps.normalize_snp_results(snp_results=snp_results, snp_aliases=yfull_tree_data.snp_aliases)
            clades = find.get_ordered_clade_details(tree=yfull_tree_data, snps=snp_results)
            
                    
        except Exception as e:
            messagebox.showerror("Error", f"yclade.find_clade failed: {e}")
            yclade_label.config(text="")
            return

        last_clades = clades

        if not clades:
            result_var.set("No clades returned by yclade.")
            yclade_label.config(text="Analysis finished - no clades returned.")
            return
        
        # האחרון הוא בעל הציון הכי גבוה ולכן בדרך כלל הכי נכון
        Final_clade = clades[Final_clade_index]
        
        ############################################################
        # Z2125 במיקום גנומי 19: 6892233 בקבצים הגולמיים של מייהירטייג חיובי לכולם ולכן גורם לבעיה בעיקר אצל AB-067 וצריך לדלג עליו לגמרי
        # השורה הלא נכונה היא: rs570569843 Y 6892233 TT
        # Known SNPs at this position: Z2125 C  to  T
        if Final_clade.name == "R-Z2125" and user_loaded and last_dna_file_info["creator"] == "myheritage":
            Final_clade = clades[Final_clade_index+1]    
        ########################################################
        
        # קבלת מידע על קבוצת אבותינו שהענף נמצא בה או שהיא נמצאת תחת הענף באמצעות פונקציה
        ab_string = get_ab_from_clade(Final_clade.name)
        
        # מוצאים את ה-score הגבוה ביותר
        #max_score = max(getattr(c, "score", 0) for c in clades)

        # לוקחים את כל הקלאדים עם ה-score הזה
        #top_clades = [c for c in clades if getattr(c, "score", 0) == max_score]

        # בטיחות: קח שדות רק אם קיימים 
        name = getattr(Final_clade, "name", "Unknown")
        age_info = getattr(Final_clade, "age_info", None)
        tmrca = getattr(age_info, "most_recent_common_ancestor", "Unknown") if age_info else "Unknown"
        formed = getattr(age_info, "formed", "Unknown") if age_info else "Unknown"

        global last_yfull_link, last_ftdna_link
        last_yfull_link = f"https://www.yfull.com/tree/{name}/" if name != "Unknown" else False
        last_ftdna_link = f"https://discover.familytreedna.com/y-dna/{name}/tree/" if name != "Unknown" else False
        
        more_result_warning = "\n\nwarning!: There is more result with the same 'score'\nOne of them may be wrong \nsave results and check it"
            
        #ref_result_label.config(fg="green", bg="SystemButtonFace")
        
        def split_ab_string(ab_string: str, per_line: int = 5) -> str:
            """
            מקבל מחרוזת ארוכה עם ערכים מופרדים בפסיקים
            ומחזיר מחרוזת מפוצלת לשורות, כך שכל שורה מכילה per_line ערכים.
            """
            items = [x.strip() for x in ab_string.split(",") if x.strip()]
            lines = []
            for i in range(0, len(items), per_line):
                lines.append(", ".join(items[i:i+per_line]))
            return "\n".join(lines)

        # פיצול רשימת ה AB לפי 5 לשורה באמצעות הפונקצייה הנל    
        ab_string_multiline = split_ab_string(ab_string, per_line=5)
        
        # הצגת התוצאה במסך
        result_var.set(
            f"run for:     {'User DNA-file' if user_loaded else 'check Y-SNP button'}\n"
            f"YFull predicted clade (used ref {last_ref_type}):\n"
            f" Name:    {name}\n"
            f" TMRCA:   {tmrca} ybp\n"
            f" FORMED:  {formed} ybp\n"
            f"{ab_string_multiline}\n"
            f"yfull tree version. {yfull_tree_data.version}"
        )

        btn_yfull.grid(row=3, column=2, padx=5, pady=5)
        btn_ftdna.grid(row=4, column=2, padx=5, pady=5)    
        # מניחים את הלחצן רק אחרי שיש תוצאות
        btn_save_results.grid(row=5, column=2, padx=5, pady=5)
         
        yclade_label.config(text="Analysis finished successfully.", fg="blue")
        
        # במקרה שיש יותר מענף אחד בציון הכי גבוה נותנת אופצייה להתקדם למיקום הבא ברשימת התוצאות
        if len(clades) >= 2 and clades[Final_clade_index].score == clades[(Final_clade_index+1)].score:
            yclade_label.config(text=more_result_warning, fg="red")
            if messagebox.askyesno("Warning", "There is more than one top-scoring clade.\nDo you want to see the next one?"):
                return run_calculate_clade(Final_clade_index = (Final_clade_index+1) % len(clades))
        '''
        # זו צורה נוספת שבודקת רק את המיקום האחרון והבא אחריו, וההודעה קובץ כל פעם מהמיקום האחרון לזה שאחריו וחוזר חלילה    
        if len(clades) >= 2 and clades[0].score == clades[(1)].score:
            yclade_label.config(text=more_result_warning, fg="red")
            if messagebox.askyesno("Warning", "There is more than one top-scoring clade.\nDo you want to see the next one?"):
                return run_calculate_clade(Final_clade_index = 1)
        '''
        '''
        # זו עוד דרך לדעת שיש יותר מענף אחד שקיבל אותו סקור. ראו להלן
        problem = False
        if len(clades) >= 2:
            next_Final_clade = clades[Final_clade_index + 1] # or clades[1]
            clade_snp = get_snp_from_clade(next_Final_clade.name)
            next_Final_clade_descendants = get_clade_and_descendants_lists(yfull_tree_data, clade_snp)
            if Final_clade.name not in next_Final_clade_descendants:
                #messagebox.showerror("next_Final_clade_descendants_Error", "next_Final_clade_descendants_Error")
                next_Final_clade_descendants_Error = "next_Final_clade_descendants_Error"
                yclade_label.config(text=next_Final_clade_descendants_Error, fg="red")
        ''' 
        
    except Exception as e:
        messagebox.showerror("Error", str(e))
        




# ------------------------
# שמירת תוצאות
# ------------------------
def save_clades_to_file():
    if not last_clades:
        messagebox.showerror("Error", "No Clades to save.")
        return
    file_path = filedialog.asksaveasfilename(
        initialfile=f"{os.path.basename(last_dna_file)}.yclade_results.txt",
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    if not file_path:
        return
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(f"Source DNA_file: {os.path.basename(last_dna_file)}\n")
            f.write(f"Reference SNPs loaded: {os.path.basename(last_reference_file)}\n\n")
            f.write(f"SNPs used for YClade analysis:\n{last_positive_snp_string}\n\n")
            f.write(f"The above SNPs can used in https://predict.yseq.net/clade-finder/index.php\n")
            
            f.write("\n\nClades returned by yclade:\n")
            f.write("\nNote: clade with highest score expected to be accurate prediction")
            f.write("\nif ultiple clades share the highest score Some of them may be incorrect.\n\n")
            
            for clade in last_clades:
                name = getattr(clade, 'name', '')
                score = getattr(clade, 'score', '')
                age = getattr(clade, 'age_info', None)
                formed = getattr(age, 'formed', '') if age else ''
                tmrca = getattr(age, 'most_recent_common_ancestor', '') if age else ''
                last_yfull_link = f"https://www.yfull.com/tree/{name}" if name != "Unknown" else ""
                # קבלת מידע על קבוצת אבותינו שהענף נמצא בה או שהיא נמצאת תחת הענף באמצעות פונקציה
                #ab_string = get_ab_from_clade(clade.name)
                line = f"Name: {name},     Score: {score}, Formed: {formed}, TMRCA: {tmrca},   Link: {last_yfull_link}\n\n"
                f.write(line)

        messagebox.showinfo("Success", f"Results saved to {file_path}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to save file: {e}")

# ------------------------
# בדיקת שם סניפ או מיקום גנומי בנתוני המשתמש וברפרנס
# ------------------------
       
def check_search_input(ref_search = True):
    
    global reference_names, reference_snaps, user_snps, user_loaded, reference_loaded, last_positive_snp_string
    
    if not reference_loaded:
        choice = messagebox.askyesnocancel("Reference not Autodetected", "Autodetected Reference faild\nChoose hg19, hg38, or select file manually.\n\nYes = hg19, No = hg38, Cancel = choose manually")
        ref_path = "Msnps_hg19.vcf.gz" if choice else "snps_hg38.vcf.gz" if (choice == False) else "ask"       
        load_reference(ref_path)
        
    if not reference_loaded:
        return
    
    search_input = entry_search.get().strip().upper()[:10] # מסירים רווחים הופכים לאותיות גדולות וחותכים את כל מה שמעבר ל 10 אותיות או ספרות
    
    if search_input.isdigit():
        pos = int(search_input)
        fields_reference = reference_snps.get(pos)
        fields_user = user_snps.get(pos) if user_loaded else False
    else:  
        pos = reference_names.get(search_input)
        fields_reference = reference_snps.get(pos)
        fields_user = user_snps.get(pos) if user_loaded else False
            
    if fields_reference:
        ref_result_var.set(fields_reference)
        ref_result_label.config(fg="green", bg="SystemButtonFace")
    else:
        ref_result_var.set(f"{search_input} not found in reference file")
        ref_result_label.config(fg="blue", bg="yellow")
        
    if fields_user:
        user_result_var.set(fields_user)
        fg_for_user_result_label = "green" if fields_user["is_positive"] == "Yes" else "red"
        user_result_label.config(fg=fg_for_user_result_label, bg="SystemButtonFace")
        
    else:
        msg = f"{search_input} not found in user DNA_file" if user_loaded else "user DNA_file_not_loaded"
        user_result_var.set(msg)
        user_result_label.config(fg="blue", bg="yellow")
        
    # במקרה שאין דנא של נבדק אז מריצים את חישוב המיקום על עץ ווייפול עבור הווריאנט המבוקש כאילו שהוא חיובי
    if not user_loaded and fields_reference:
        last_positive_snp_string = f"{fields_reference['name']}+"
        run_calculate_clade()
        
        
# ------------------------
# הדבקה מהלוח
# ------------------------
def paste_from_clipboard():
    clipboard_text = root.clipboard_get().strip()[:15]
    entry_search.delete(0, tk.END)
    entry_search.insert(0, clipboard_text)
    
    
def unload_ref():
    global reference_snps, reference_names, last_reference_file, reference_loaded
    reference_snps = {}      # pos -> ref_snp_name + ref_allele
    reference_names = {}      # id_snp_name -> pos
    last_reference_file = ""
    reference_loaded = False
    reference_loading_label.config(text="No reference_file loaded", fg="red")
    btn_unload_ref.grid_forget()
    ref_result_var.set("")
    ref_result_label.config(text="", bg="SystemButtonFace")
    reset_user()


# ------------------------
# GUI
# ------------------------
root = tk.Tk()
# קביעת גודל התחלתי (רוחב x גובה)
#root.geometry("650x500")

# קביעת מינימום גודל
#root.minsize(500, 500)

root.title("Y_Tree SNP Analyzer | by Dr. simcha-gershon Bohrer (Phd.) | versin: 26 sep 2025")

# מפריד אנכי
ttk.Separator(root, orient="vertical").grid(row=0, column=1, sticky="ns", padx=5, rowspan=20)
ttk.Separator(root, orient="vertical").grid(row=0, column=3, sticky="ns", padx=5, rowspan=20)


tk.Label(root, text="Reference file", font="david 14 bold").grid(row=0, column=0, padx=75, pady=10)
tk.Label(root, text="Yclade", font="david 14 bold").grid(row=0, column=2, padx=75, pady=10)
tk.Label(root, text="User DNA file", font="david 14 bold").grid(row=0, column=4, padx=75, pady=10)


# יצירת כפתורי רדיו
# משתנה בחירה (מחזיק את הערך של הכפתור הנבחר)
ref_var = tk.StringVar(value="Aautodetect")
tk.Radiobutton(root, text="Aautodetect Reference File     ", variable=ref_var, value="Aautodetect").grid(row=2, column=0)
tk.Radiobutton(root, text="Ask Reference File (vcf/vcf.gz)", variable=ref_var, value="ask").grid(row=3, column=0)

# כפתור לביטול טעינת קובץ הרפרנס
btn_unload_ref = tk.Button(root, text="unload_ref_file", command=unload_ref)
# המיקום גריד שלו מתבצע בפונקציית טעינת הרפרנס

# כפתור לביטול טעינת קובץ דנא של המשתמש
btn_unload_dna = tk.Button(root, text="unload_dna_file", command=reset_user)
# המיקום גריד שלו מתבצע בפונקציית טעינת קובץ דנא של המשתמש

btn_csv = tk.Button(root, text="Choose \nUser RAW-DNA File \nvcf/vcf.gz/txt/csv/gz/zip", command=load_dna_file)
btn_csv.grid(row=1, column=4, rowspan=2)

reference_loading_label = tk.Label(root, text="No reference-file loaded", fg="red")
reference_loading_label.grid(row=4, column=0, padx=5, pady=5, rowspan=3)

dna_loading_label = tk.Label(root, text="No DNA-file loaded", fg="red")
dna_loading_label.grid(row=4, column=4, padx=5, pady=5, rowspan=3)

yclade_label = tk.Label(root, text="Check SNP or load DNA-file", anchor="w", fg="red")
yclade_label.grid(row=1, column=2, padx=5, pady=5)

result_var = tk.StringVar()
result_label = tk.Label(root, textvariable=result_var, fg="green")
result_label.grid(row=2, column=2, padx=5, pady=10)

def open_yfull():
    webbrowser.open_new(last_yfull_link)
def open_ftdna():
    webbrowser.open_new(last_ftdna_link)

btn_yfull = tk.Button(root, text="open clade in Yfull Tree", command=open_yfull)
btn_ftdna = tk.Button(root, text="open clade in FTDNA Tree", command=open_ftdna)
btn_save_results = tk.Button(root, text="Save Clades to TXT", command=save_clades_to_file)
# הגריד שלהם נמצא בפונקציית חישוב הקלייד

# בדיקת מיקום גנומי עם הכיתוב hg38
ttk.Separator(root, orient="horizontal").grid(row=10, column=2, sticky="ew", padx=5, pady=30)
tk.Label(root, text="Check Y-SNP", font="david 14 bold").grid(row=11, column=2, padx=5, pady=5)

# כפתור הדבקה משמאל עם בדיקה
btn_paste = tk.Button(root, text="Paste", command=paste_from_clipboard)
btn_paste.grid(row=12, column=2, padx=5, pady=5)

# בדיקת מיקום גנומי (עם תווית שמציינת hg38/hg19)
entry_search = tk.Entry(root, width=15)
entry_search.grid(row=13, column=2, padx=5, pady=5)

# כפתור בדיקה בנתוני המשתמש
btn_check = tk.Button(root, text="Check: SNP name / Genomic position", command=check_search_input) 
btn_check.grid(row=14, column=2, padx=5, pady=5)
    
ref_result_var = tk.StringVar()
ref_result_label = tk.Label(root, textvariable=ref_result_var, fg="green")
ref_result_label.grid(row=14, column=0)

user_result_var = tk.StringVar()
user_result_label = tk.Label(root, textvariable=user_result_var, fg="green")
user_result_label.grid(row=14, column=4)

tk.Label(root, text="NOTE: Each reference has different positions").grid(row=15, column=2, padx=5, pady=5)

get_ab_data()


root.mainloop()


