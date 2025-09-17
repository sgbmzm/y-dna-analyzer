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
last_positive_snps_list = []
last_positive_snp_string = ""
reference_snps = {}      # pos -> ref_snp_name + ref_allele
reference_names = {}      # id_snp_name -> pos
user_snaps = {}          # pos -> ref_snp_name + ref_allele
current_dna_file = ""
current_reference_file = ""
reference_loaded = False
user_loaded = False

# ------------------------
# GUI כפתורים
# ------------------------
def set_buttons_state(state):
    btn_calculate_clade.config(state=state)
    btn_check_user.config(state=state)
    btn_check_ref.config(state=state)
    btn_paste.config(state=state)
    #if last_clades:
    btn_save_results.config(state=state)

def update_buttons_state():
    if reference_loaded and user_loaded:
        set_buttons_state("normal")
    else:
        set_buttons_state("disabled")

# ------------------------
# איפוס
# ------------------------
def reset_all():
    global last_clades, last_positive_snp_string, last_positive_snps_list
    global reference_snps, reference_names, current_dna_file, reference_loaded, user_loaded
    
    last_clades = []
    last_positive_snps_list = []
    last_positive_snp_string = ""
    reference_snps = {}
    reference_names={}
    current_dna_file = ""
    current_reference_file = ""
    reference_loaded = False
    user_loaded = False
    #ref_var.set(value="Aautodetect")

    result_var.set("")
    link_label_yfull.config(text="", fg="blue", cursor="")
    link_label_ftdna.config(text="", fg="blue", cursor="")
    dna_loading_label.config(text="No DNA_file loaded", fg="red")
    reference_loading_label.config(text="No reference_file loaded", fg="red")
    info_label.config(text="")
    update_buttons_state()
    
    

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
    global reference_loaded, reference_snps, reference_names, current_reference_file
    
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
        reference_loading_label.config(fg="blue", text=f"Reference loaded: {os.path.basename(ref_path)}     {len(reference_snps)} Y-SNPs")
        current_reference_file = os.path.basename(ref_path)
        update_buttons_state()
        
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load reference: {e}")
        #reference_loading_label.config(text="")
        
        

# ------------------------
# ניתוח CSV
# ------------------------

        
def load_dna_file():
    
    # כל פעם כשבוחרים קובץ נא חדש הכל מתאפס ומחושב מהתחלה
    reset_all()
    
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
    
    # טעינת הרפרנס הנבחר    
    load_reference(ref_path)
    
    # הצהרה על משתנים גלובליים הדרושים להלן
    global last_clades, last_positive_snp_string, last_positive_snps_list, current_dna_file, user_snps, user_loaded, reference_loaded
    
    # אם לא בחרו קובץ רפרנס מאפסים הכל ולא ממשיכים
    if not reference_loaded:
        reset_all()
        return

    last_clades = []
    last_positive_snps_list = []
    last_positive_snp_string = ""
    current_dna_file = file_path

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
                    chrom, pos_str, rsid, ref, alt, qual = parts[:6]
                    sample_index = 0 # לפעמים יש כמה נבדקים שונים וכל אחד מהם בטור אחד אחרי השני כאן מדובר באחד בלבד
                    gt_str = parts[9 + sample_index].split(":")[0] # לוקח רק את החלק של ה־GT (למשל "0/1" או "1/1" מתוך העמודה של הנבדק הראשון
                    alleles = [int(a) for a in re.split("[/|]", gt_str) if a.isdigit()]
                    if any(a > 0 for a in alleles): # כלומר אם יש שם משהו שהוא 1 אז יש משהו חיובי וזה אומר שהוא כמו ה alt
                        is_positive = True
                        alleles_str = alt
                    else:
                        is_positive = False
                        alleles_str = ref
                else:
                    # נניח שהשדות הראשונים תמיד הם: rsid, chrom, pos, alleles
                    rsid, chrom, pos_str, alleles_str = parts[:4]
                    ref = "?"
                
                # עושה אלל אחד גדול גם היה קטן וגם אם היו שניים יחד
                allele_str = alleles_str.upper().strip()[0] if alleles_str else ""
            
                # נרמול שם הכרומוזום ל־Y
                chrom = chrom.replace("chr", "").upper()
                if chrom not in ("Y", "24"):
                    continue
                
                # הוספת השורה למילון המשתמש לאחר שוודאנו שמובר בשורה של Y
                user_snps[int(pos_str)] = {"chrom": chrom, "pos_str": pos_str, "snp_name": "?", "ref": "?", "alt": allele_str, "is_positive": "?"}
                
                # בדיקה מול הרפרנס
                if not pos_str.isdigit():
                    continue
                ref_info = reference_snps.get(int(pos_str))
                if not ref_info:
                    continue

                if allele_str == ref_info["alt"]:
                    s = f"{ref_info['name']}+"
                    positive_snps.append(s)
                    last_positive_snps_list.append(s)
                    user_snps[int(pos_str)]["is_positive"] = "Yes"   # או "כן" / "+" או כל מה שאתה רוצה
                    user_snps[int(pos_str)]["snp_name"] = ref_info['name']
                elif allele_str == ref_info["ref"] or allele_str == ".":
                    user_snps[int(pos_str)]["is_positive"] = "No"
                    user_snps[int(pos_str)]["snp_name"] = ref_info['name']

        last_positive_snp_string = ", ".join(positive_snps)
        #print("Final positive SNPs:", last_positive_snps_list)

        dna_loading_label.config(fg="blue", text=f"DNA file loaded: {os.path.basename(file_path)}     {len(user_snps)} total Y-rows,     {len(positive_snps)} Positive Y-SNPs,     in DNA_file")
        user_loaded = True
        update_buttons_state()

    except Exception as e:
        messagebox.showerror("Error", f"Failed reading file: {e}")
        loading_label.config(text="")

                  
def run_calculate_clade():
    global last_clades

    if not last_positive_snp_string:
        result_var.set("No matching SNPs were found in the file")
        loading_label.config(text="Analysis finished - no matching SNPs found.")
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
            loading_label.config(text="Analysis finished - no clades returned.")
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
            f"YFull Final clade:\n"
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
        
        btn_save_results.config(state="normal")

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
        initialfile=f"{os.path.basename(current_dna_file)}.yclade_results.txt",
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    if not file_path:
        return
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(f"Source DNA_file: {os.path.basename(current_dna_file)}\n")
            f.write(f"Reference SNPs loaded: {os.path.basename(current_reference_file)}\n\n")
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
# בדיקת מיקום גנומי מהירה (עם שימוש במיפוי העמודות)
# ------------------------
def check_position():
    if not user_loaded:
        return
    
    pos_input = entry_search.get().strip()
    
    if pos_input.isdigit():
        pos = int(pos_input)
        fields = user_snps.get(pos)
        
    if fields:
        user_result_var.set(fields)
        
    else:
        user_result_var.set(f"Position {pos} not found in user DNA_file")
        user_result_label.config(fg="blue", bg="red")
        
        
def check_search_input(ref_search = True):
    if not reference_loaded:
        return
    
    global reference_names, reference_snaps, user_snps
    
    search_input = entry_search.get().strip().upper()
    
    if search_input.isdigit():
        pos = int(search_input)
        fields_reference = reference_snps.get(pos)
        fields_user = user_snps.get(pos)
    
    else:  
        pos = reference_names.get(search_input)
        fields_user = user_snps.get(pos)
        fields_reference = reference_snps.get(pos)
            
    if fields_reference:
        ref_result_var.set(fields_reference)
    else:
        ref_result_var.set(f"{search_input} not found in reference file")
        ref_result_label.config(fg="blue", bg="red")
        
    if fields_user:
        user_result_var.set(fields_user)
    else:
        user_result_var.set(f"{search_input} not found in user DNA_file")
        user_result_label.config(fg="blue", bg="red")
    
        
        
# ------------------------
# הדבקה מהלוח
# ------------------------
def paste_from_clipboard():
    clipboard_text = root.clipboard_get().strip()
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

btn_csv = tk.Button(root, text=" Click To Choose RAW_DNA_File & Analyze", command=load_dna_file)
btn_csv.grid(row=1, column=1, padx=5, pady=5)

reference_loading_label = tk.Label(root, text="No reference_file loaded", fg="red")
reference_loading_label.grid(row=2, column=0, columnspan=3, padx=5, pady=5)

dna_loading_label = tk.Label(root, text="No DNA_file loaded", fg="red")
dna_loading_label.grid(row=3, column=0, columnspan=3, padx=5, pady=5)

info_label = tk.Label(root, text="", anchor="w", justify="left")
info_label.grid(row=4, column=0, columnspan=3, sticky="we", padx=5, pady=5)

btn_calculate_clade = tk.Button(root, text="Calculate Clade", command=run_calculate_clade, state="disabled")
btn_calculate_clade.grid(row=5, column=0, columnspan=3, pady=5)

result_var = tk.StringVar()
result_label = tk.Label(root, textvariable=result_var, justify="left", fg="green")
result_label.grid(row=6, column=0, columnspan=3, padx=5, pady=10)

link_label_yfull = tk.Label(root, text="", fg="blue", cursor="")
link_label_yfull.grid(row=7, column=0, columnspan=3, padx=5, pady=5)

link_label_ftdna = tk.Label(root, text="", fg="blue", cursor="")
link_label_ftdna.grid(row=8, column=0, columnspan=3, padx=5, pady=5)


btn_save_results = tk.Button(root, text="Save Clades to TXT", command=save_clades_to_file, state="disabled")
btn_save_results.grid(row=9, column=0, columnspan=3, pady=5)


# בדיקת מיקום גנומי (עם תווית שמציינת hg38/hg19)
#tk.Label(root, text="Genomic position (numeric):").grid(row=10, column=0, sticky="e")
entry_search = tk.Entry(root, width=20)
entry_search.grid(row=10, column=1, padx=5, pady=5)

# כפתור הדבקה משמאל עם בדיקה
btn_paste = tk.Button(root, text="Paste", command=paste_from_clipboard, state="disabled")
btn_paste.grid(row=10, column=2, padx=5, pady=5)

# כפתור בדיקה בנתוני המשתמש
btn_check_user = tk.Button(root, text="Check Y-SNP in user DNA_file", command=check_search_input, state="disabled")
btn_check_user.grid(row=10, column=3, padx=5, pady=5)

user_result_var = tk.StringVar()
user_result_label = tk.Label(root, textvariable=user_result_var, justify="left", fg="green")
user_result_label.grid(row=11, column=3, columnspan=2, padx=5, pady=10)
    
# כפתור בדיקה בנתוני הרפרנס
btn_check_ref = tk.Button(root, text="Check Y-SNP in reference file", command=check_search_input, state="disabled")
btn_check_ref.grid(row=10, column=0, padx=5, pady=5)

ref_result_var = tk.StringVar()
ref_result_label = tk.Label(root, textvariable=ref_result_var, justify="left", fg="green")
ref_result_label.grid(row=11, column=0, columnspan=2, padx=5, pady=10)

root.mainloop()

