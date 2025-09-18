import tkinter as tk
from tkinter import filedialog, messagebox
import csv
import os
import re
import gzip
import yclade
import webbrowser

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
    link_label_yfull.config(text="", fg="blue", cursor="")
    link_label_ftdna.config(text="", fg="blue", cursor="")
    dna_loading_label.config(text="No DNA_file loaded", fg="red")
    info_label.config(text="")
    user_result_var.set("")
    user_result_label.config(text="", bg="SystemButtonFace")
    btn_save_results.grid_forget()
    btn_unload_ref.grid_forget()
    last_dna_file_type = ""
    
    
# פונקצייה לבדיקה האם כתוב בקובץ באיזה רפרנס הוא משתמש
def detect_reference(file_path):
    
    ref_map = {
        "hg19": ["build 37", "grch37", "hg19"],
        "hg38": ["build 38", "grch38", "hg38"],
    }
    
    opener = gzip.open if file_path.endswith(".gz") else open
    with opener(file_path, "rt", encoding="utf-8") as f:
        for line in f:
            if not line.startswith("#"):
                break  # יציאה – כבר עברנו את ה-header
            lower_line = line.lower()
            for ref, keys in ref_map.items():
                if any(k in lower_line for k in keys):
                    return ref
    return None  # לא זוהה


# ------------------------
# טעינת רפרנס
# ------------------------
def load_reference(ref_path="ask"):
    global reference_loaded, reference_snps, reference_names, last_reference_file, last_ref_type
    
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
    
    last_ref_type = detect_reference(ref_path) 

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
        reference_loading_label.config(fg="blue", text=f"Reference loaded: {os.path.basename(ref_path)}  |  {len(reference_snps)} Y-SNPs  |  type: {last_ref_type}")
        last_reference_file = os.path.basename(ref_path)
        btn_unload_ref.grid(row=2, column=0, padx=5, pady=5)
        #update_buttons_state()
        
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
        filetypes=[("DNA files", "*.txt *.csv *.vcf *.vcf.gz"), ("All files", "*.*")]
    )
    if not file_path:
        return
    
    # בדיקת או בחירת הרפרנס המתאים לקובץ הדנא של המשתמש
    if ref_var.get() == "Aautodetect":
        ref_auto_detect = detect_reference(file_path)
        ref_path = "Msnps_hg19.vcf.gz" if ref_auto_detect == "hg19" else "snps_hg38.vcf.gz" if ref_auto_detect == "hg38" else "ask"
    
        if not ref_auto_detect:
            choice = messagebox.askyesnocancel("Reference not Autodetected", "Autodetected Reference faild\nChoose hg19, hg38, or select file manually.\n\nYes = hg19, No = hg38, Cancel = choose manually")
            ref_path = "Msnps_hg19.vcf.gz" if choice else "snps_hg38.vcf.gz" if (choice == False) else "ask"
        
        # במקרה שקובץ הרפרנס לא נמצא בתיקיית הסקריפט הנוכחי יש לבחור קובץ באופן ידני
        if ref_path != "ask" and not os.path.exists(ref_path):
            messagebox.showwarning("Reference file missing", f"Reference file:     {ref_path}     missing \nplease select file manually")
            ref_path = "ask"       
    else:
        ref_path = "ask"
    
    global last_reference_file, reference_loaded
        
    if not last_reference_file == ref_path:
        reference_loaded = False
    
    if not reference_loaded:    
        # טעינת הרפרנס הנבחר    
        load_reference(ref_path)
    
    # הצהרה על משתנים גלובליים הדרושים להלן
    global last_clades, last_positive_snp_string, last_dna_file, user_snps, user_loaded, last_dna_file_type
    
    # אם לא בחרו קובץ רפרנס מאפסים הכל ולא ממשיכים
    if not reference_loaded:
        reset_user()
        return

    last_clades = []
    last_positive_snp_string = ""
    last_dna_file = file_path
    last_dna_file_type = ref_auto_detect

    dna_loading_label.config(text=f"Loading {os.path.basename(file_path)} ...")
    root.update()

    positive_snps = []
    user_snps = {}
    
    is_vcf_file = file_path.endswith(".vcf") or file_path.endswith(".vcf.gz")
    
    opener = gzip.open if file_path.endswith(".gz") else open
    try:  
        with opener(file_path, "rt", encoding="utf-8") as f:
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

        dna_loading_label.config(fg="blue", text=f"DNA file loaded: {os.path.basename(file_path)} {len(user_snps)} total Y-rows  |  {len(positive_snps)} Positive Y-SNPs in DNA_file  |  type: {last_dna_file_type}")
        user_loaded = True
        
        run_calculate_clade()

    except Exception as e:
        messagebox.showerror("Error", f"Failed reading file: {e}")
        dna_loading_label.config(text="")

                  
def run_calculate_clade():
    global last_clades, last_reference_file, ref_user_file, last_dna_file_type, last_ref_type

    if not last_positive_snp_string:
        result_var.set("No matching SNPs were found in the file")
        info_label.config(text="Analysis finished - no matching SNPs found.")
        return

    try:
        info_label.config(text="Running clade calculation...")
        root.update()

        # קריאה ל־yclade עם מחרוזת אחת
        try:
            clades = yclade.find_clade(last_positive_snp_string)
        except Exception as e:
            messagebox.showerror("Error", f"yclade.find_clade failed: {e}")
            info_label.config(text="")
            return

        last_clades = clades

        if not clades:
            result_var.set("No clades returned by yclade.")
            info_label.config(text="Analysis finished - no clades returned.")
            return

        Final_clade = clades[0]
            
        # מוצאים את ה-score הגבוה ביותר
        max_score = max(getattr(c, "score", 0) for c in clades)

        # לוקחים את כל הקלאדים עם ה-score הזה
        top_clades = [c for c in clades if getattr(c, "score", 0) == max_score]

        # בטיחות: קח שדות רק אם קיימים 
        name = getattr(Final_clade, "name", "Unknown")
        age_info = getattr(Final_clade, "age_info", None)
        tmrca = getattr(age_info, "most_recent_common_ancestor", "Unknown") if age_info else "Unknown"
        formed = getattr(age_info, "formed", "Unknown") if age_info else "Unknown"

        YFull_Final_clade_Url = f"https://www.yfull.com/tree/{name}/" if name != "Unknown" else ""
        FTDNA_Final_clade_Url = f"https://discover.familytreedna.com/y-dna/{name}/tree/" if name != "Unknown" else ""
        
        warning = "\n\nwarning!: There is more result with the same 'score'\nOne of them may be wrong \nsave results and check it"

        # הצגת התוצאה במסך
        result_var.set(
            f"ref_in_dna_file:     {last_dna_file_type}\n"
            f"YFull predicted clade (used ref {last_ref_type}):\n"
            f" Name:    {name}\n"
            f" TMRCA:   {tmrca} ybp\n"
            f" FORMED:  {formed} ybp"
            f"{warning if clades[0].score == clades[1].score else ''}"
        )

        if YFull_Final_clade_Url:
            link_label_yfull.config(text=f' LINK:     {YFull_Final_clade_Url}', fg="blue", cursor="hand2")
            link_label_yfull.bind("<Button-1>", lambda e: webbrowser.open_new(YFull_Final_clade_Url))
        
        if FTDNA_Final_clade_Url:
            link_label_ftdna.config(text=f' LINK:     {FTDNA_Final_clade_Url}', fg="blue", cursor="hand2")
            link_label_ftdna.bind("<Button-1>", lambda e: webbrowser.open_new(FTDNA_Final_clade_Url))
                
        info_label.config(text="Analysis finished successfully.")
        
        # מניחים את הלחצן רק אחרי שיש תוצאות
        btn_save_results.grid(row=9, column=0, columnspan=3, pady=5)

    except Exception as e:
        messagebox.showerror("Error", str(e))
        #info_label.config(text="")




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
                YFull_Final_clade_Url = f"https://www.yfull.com/tree/{name}" if name != "Unknown" else ""
                line = f"Name: {name},     Score: {score}, Formed: {formed}, TMRCA: {tmrca},     Link: {YFull_Final_clade_Url}\n\n"
                f.write(line)

        messagebox.showinfo("Success", f"Results saved to {file_path}")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to save file: {e}")

# ------------------------
# בדיקת שם סניפ או מיקום גנומי בנתוני המשתמש וברפרנס
# ------------------------
       
def check_search_input(ref_search = True):
    
    global reference_names, reference_snaps, user_snps, user_loaded, reference_loaded
    
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
        ref_result_label.config(fg="blue", bg="red")
        
    if fields_user:
        user_result_var.set(fields_user)
        user_result_label.config(fg="green", bg="SystemButtonFace")
        
    else:
        msg = f"{search_input} not found in user DNA_file" if user_loaded else "user DNA_file_not_loaded"
        user_result_var.set(msg)
        user_result_label.config(fg="blue", bg="red")
    
        
        
# ------------------------
# הדבקה מהלוח
# ------------------------
def paste_from_clipboard():
    clipboard_text = root.clipboard_get().strip()[:15]
    entry_search.delete(0, tk.END)
    entry_search.insert(0, clipboard_text)

# ------------------------
# GUI
# ------------------------
root = tk.Tk()
root.title("Y_Tree SNP Analyzer | by Dr. simcha-gershon Bohrer (Phd.) | versin: 17 sep 2025")

# יצירת כפתורי רדיו
# משתנה בחירה (מחזיק את הערך של הכפתור הנבחר)
ref_var = tk.StringVar(value="Aautodetect")
tk.Radiobutton(root, text="Aautodetect Reference File", variable=ref_var, value="Aautodetect").grid(row=0, column=0, padx=5, pady=5)
tk.Radiobutton(root, text="Manual reference File (VCF/VCF.GZ)", variable=ref_var, value="ask").grid(row=0, column=1, padx=5, pady=5)

def unload_ref():
    global reference_snps, reference_names, last_reference_file, reference_loaded
    reference_snps = {}      # pos -> ref_snp_name + ref_allele
    reference_names = {}      # id_snp_name -> pos
    last_reference_file = ""
    reference_loaded = False
    reference_loading_label.config(text="No reference_file loaded", fg="red")
    btn_unload_ref.grid_forget()
    reset_user()

# כפתור הדבקה משמאל עם בדיקה
btn_unload_ref = tk.Button(root, text="unload", command=unload_ref)
# המיקום גריד שלו מתבצע בפונקציית 

btn_csv = tk.Button(root, text=" Click To Choose RAW_DNA_File & Analyze", command=load_dna_file)
btn_csv.grid(row=1, column=1, padx=5, pady=5)

reference_loading_label = tk.Label(root, text="No reference_file loaded", fg="red")
reference_loading_label.grid(row=2, column=1, padx=5, pady=5)

dna_loading_label = tk.Label(root, text="No DNA_file loaded", fg="red")
dna_loading_label.grid(row=3, column=0, columnspan=3, padx=5, pady=5)

info_label = tk.Label(root, text="", anchor="w", justify="left")
info_label.grid(row=4, column=0, columnspan=3, padx=5, pady=5)

result_var = tk.StringVar()
result_label = tk.Label(root, textvariable=result_var, justify="left", fg="green")
result_label.grid(row=6, column=0, columnspan=3, padx=5, pady=10)

link_label_yfull = tk.Label(root, text="", fg="blue", cursor="")
link_label_yfull.grid(row=7, column=0, columnspan=3, padx=5, pady=5)

link_label_ftdna = tk.Label(root, text="", fg="blue", cursor="")
link_label_ftdna.grid(row=8, column=0, columnspan=3, padx=5, pady=5)


btn_save_results = tk.Button(root, text="Save Clades to TXT", command=save_clades_to_file)
# המיקום גריד שלו מתבצע בפונקציית כלכולייט קלייד

# בדיקת מיקום גנומי עם הכיתוב hg38
tk.Label(root, text="_____________________________").grid(row=10, column=1, sticky="e")

# כפתור הדבקה משמאל עם בדיקה
btn_paste = tk.Button(root, text="Paste", command=paste_from_clipboard)
btn_paste.grid(row=12, column=2, padx=5, pady=5)

# בדיקת מיקום גנומי עם הכיתוב hg38
tk.Label(root, text="Reference:").grid(row=13, column=0, sticky="we")

# בדיקת מיקום גנומי עם הכיתוב hg38
tk.Label(root, text="User:").grid(row=13, column=3, sticky="we")

# בדיקת מיקום גנומי (עם תווית שמציינת hg38/hg19)
entry_search = tk.Entry(root, width=15)
entry_search.grid(row=13, column=1, padx=5, pady=5)

# כפתור בדיקה בנתוני המשתמש
btn_check = tk.Button(root, text="Check Y-SNP", command=check_search_input) 
btn_check.grid(row=14, column=1, padx=5, pady=5)
    
ref_result_var = tk.StringVar()
ref_result_label = tk.Label(root, textvariable=ref_result_var, justify="left", fg="green")
ref_result_label.grid(row=14, column=0, padx=5, pady=10)

user_result_var = tk.StringVar()
user_result_label = tk.Label(root, textvariable=user_result_var, justify="left", fg="green")
user_result_label.grid(row=14, column=2, padx=5, pady=10)

root.mainloop()


