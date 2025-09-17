import tkinter as tk
from tkinter import filedialog, messagebox
import re
import yclade
import webbrowser
import os
import gzip

# ------------------------
# משתנים גלובליים
# ------------------------
last_clades = []
last_positive_snp_string = ""
vcf_data = []
vcf_dict = {}  # key = position (int), value = fields
vcf_loaded = False
reference_snps = {}  # key = position (int), value = SNP name

# ------------------------
# טעינת קובץ רפרנס SNP (snps_hg38.vcf)
# ------------------------
def load_reference_snps(ref_path="snps_hg38.vcf.gz"):
    global reference_snps
    reference_snps = {}
    try:
        with gzip.open(ref_path, "rt" ,encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                chrom = fields[0]
                if chrom != "chrY":
                    continue
                pos = int(fields[1])
                id_field = fields[2]
                if not id_field or id_field == ".":
                    continue
                reference_snps[pos] = id_field
        print(f"Reference SNPs loaded: {len(reference_snps)} entries")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load reference SNPs file: {e}")

# ------------------------
# פונקציה לעדכון מצב כפתורים
# ------------------------
def set_buttons_state(state):
    btn_calculate_clade.config(state=state)
    btn_save_results.config(state=state)
    btn_paste.config(state=state)
    btn_check_position.config(state=state)

# ------------------------
# קריאה חד-פעמית של VCF (רגיל או gz) וטעינת מילון
# ------------------------
def load_vcf_once():
    global vcf_data, vcf_loaded, vcf_dict
    if not vcf_loaded:
        vcf_path = entry_file.get()
        if not vcf_path:
            messagebox.showerror("Error", "No VCF file selected")
            return False
        # השבתת כפתורים בזמן טעינה
        set_buttons_state("disabled")
        loading_label.config(text=f"Loading {os.path.basename(vcf_path)}...\nPlease wait until analysis is finished. Buttons will be enabled afterwards.")
        root.update()
        try:
            if vcf_path.endswith(".gz"):
                with gzip.open(vcf_path, "rt", encoding="utf-8") as f:
                    vcf_data = f.readlines()
            else:
                with open(vcf_path, encoding="utf-8") as f:
                    vcf_data = f.readlines()
            # בונה מילון מיקום -> fields
            vcf_dict = {}
            for line in vcf_data:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                chrom = fields[0]
                if chrom != "chrY":
                    continue
                pos_field = int(fields[1])
                vcf_dict[pos_field] = fields
            vcf_loaded = True
            
            # טעינת קובץ הרפרנס עם השמות
            load_reference_snps("snps_hg38.vcf.gz")
        
            # הפעלת כפתורים לאחר סיום טעינת המילון
            set_buttons_state("normal")
            loading_label.config(text="VCF loaded successfully. You may now use all buttons.")
            return True
        except Exception as e:
            messagebox.showerror("Error", str(e))
            set_buttons_state("normal")
            loading_label.config(text="")
            return False
    return True


# ------------------------
# בחירת קובץ VCF
# ------------------------
def choose_file():
    file_path = filedialog.askopenfilename(
        title="VCF file select",
        filetypes=[("VCF files", "*.vcf;*.vcf.gz"), ("All files", "*.*")]
    )
    if file_path:
        entry_file.delete(0, tk.END)
        entry_file.insert(0, file_path)
        
        # טעינת הקובץ ומילון
        load_vcf_once()

# ------------------------
# ניתוח קובץ VCF
# ------------------------
def analyze_vcf():
    
    global last_clades, last_positive_snps, last_positive_snp_string, last_snp_string
    if not load_vcf_once():
        return
    try:
         
        positive_snps = []
        sample_index=0
        for pos, fields in vcf_dict.items():
            # לוקחים שם SNP מהרפרנס לפי פוזיציה
            id_field = reference_snps.get(pos, None)
            if not id_field:
                continue
            gt_str = fields[9 + sample_index].split(":")[0]
            alleles = [int(a) for a in re.split("[/|]", gt_str) if a.isdigit()]
            
            if any(a > 0 for a in alleles):
                s = f"{id_field}+ "
                positive_snps.append(s.strip())
        
        print("positive_snps_len", len(positive_snps))        
        print(positive_snps)
            
        positive_snp_string = ", ".join(positive_snps)
        last_positive_snp_string = positive_snp_string
        
        if positive_snps:
            clades = yclade.find_clade(positive_snp_string) # לא לעשות snp_string כי זה תוקע המון
            last_clades = clades
            Final_clade = clades[0]
            YFull_Final_clade_Url = f"https://www.yfull.com/tree/{Final_clade.name}/"
            result_var.set(
                f"YFull Final clade:\n"
                f" Name:    {Final_clade.name}\n"
                f" TMRCA:   {Final_clade.age_info.most_recent_common_ancestor} ybp\n"
                f" FORMED:  {Final_clade.age_info.formed} ybp"
            )
            link_label.config(text=f' LINK:     {YFull_Final_clade_Url}', fg="blue", cursor="hand2")
            link_label.bind("<Button-1>", lambda e: webbrowser.open_new(YFull_Final_clade_Url))
        else:
            result_var.set("No matching SNPs were found in the file")
    except Exception as e:
        messagebox.showerror("Error", str(e))

# ------------------------
# שמירת תוצאות לקובץ TXT כולל SNPs ל-yclade
# ------------------------
def save_results_to_file(clades, last_positive_snp_string, default_name):
    file_path = filedialog.asksaveasfilename(
        initialfile=default_name,
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
    )
    if not file_path:
        return
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("SNPs used for YClade analysis:\n")
        f.write(f"{last_positive_snp_string}\n\n")
        for clade in clades:
            f.write(f"Name: {clade.name}\n")
            f.write(f"Score: {clade.score}\n")
            f.write(f"Formed: {clade.age_info.formed}\n")
            f.write(f"Formed CI: {clade.age_info.formed_confidence_interval}\n")
            f.write(f"TMRCA: {clade.age_info.most_recent_common_ancestor}\n")
            f.write(f"TMRCA CI: {clade.age_info.most_recent_common_ancestor_confidence_interval}\n")
            f.write("---\n")
    messagebox.showinfo("Success", f"Results saved to {file_path}")

def save_button_action():
    if not last_clades:
        messagebox.showerror("Error", "No results to save. Run analysis first.")
        return
    vcf_name = os.path.basename(entry_file.get())
    default_save_name = os.path.splitext(vcf_name)[0].replace(".gz","") + "_yclade.txt"
    save_results_to_file(last_clades, last_positive_snp_string, default_save_name)

# ------------------------
# בדיקת מיקום גנומי מהירה
# ------------------------
def check_position():
    if not load_vcf_once():
        return
    pos_input = entry_pos.get().strip()
    if not pos_input.isdigit():
        messagebox.showerror("Error", "Please enter a valid numeric position")
        return
    pos = int(pos_input)
    fields = vcf_dict.get(pos)
    if fields:
        gt_str = fields[9].split(":")[0]
        alleles = [int(a) for a in re.split("[/|]", gt_str) if a.isdigit()]
        allele_str = fields[4] if any(a > 0 for a in alleles) else fields[3]
        # שם SNP מהרפרנס
        snp_name = reference_snps.get(pos, "Unknown")
        if any(a > 0 for a in alleles):
            pos_result_var.set(f"Position {pos} ({snp_name}): Positive (+) ({allele_str})")
            pos_result_label.config(fg="green", bg=root.cget("bg"))
        else:
            pos_result_var.set(f"Position {pos} ({snp_name}): Negative (-) ({allele_str})")
            pos_result_label.config(fg="red", bg=root.cget("bg"))
    else:
        pos_result_var.set(f"Position {pos} not found in VCF")
        pos_result_label.config(fg="blue", bg="red")

# ------------------------
# GUI ראשי
# ------------------------
root = tk.Tk()
root.title("Y-Clade Finder - YFull")

# שורת קובץ
tk.Label(root, text="VCF File:").grid(row=0, column=0, sticky="e")
entry_file = tk.Entry(root, width=50)
entry_file.grid(row=0, column=1, padx=5, pady=5)
btn_select_file = tk.Button(root, text="Select VCF file", command=choose_file)
btn_select_file.grid(row=0, column=2, padx=5, pady=5)

# כפתור חישוב Clade
btn_calculate_clade = tk.Button(root, text="Calculate Clade", command=analyze_vcf)
btn_calculate_clade.grid(row=1, column=0, columnspan=3, pady=10)

# תוצאה
result_var = tk.StringVar()
result_label = tk.Label(root, textvariable=result_var, justify="left", fg="red")
result_label.grid(row=2, column=0, columnspan=3, padx=5, pady=10)

# לינק לחיץ
link_label = tk.Label(root, text="", fg="blue", cursor="")
link_label.grid(row=3, column=0, columnspan=3, padx=5, pady=5)

# כפתור שמירה
btn_save_results = tk.Button(root, text="Save Results to TXT", command=save_button_action)
btn_save_results.grid(row=4, column=0, columnspan=3, pady=10)

# בדיקת מיקום גנומי עם הכיתוב hg38
tk.Label(root, text="Genomic position (hg38):").grid(row=5, column=0, sticky="e")

entry_pos = tk.Entry(root, width=20)
entry_pos.grid(row=5, column=2, padx=5, pady=5)

# כפתור הדבקה משמאל עם בדיקה
def paste_from_clipboard():
    try:
        clipboard_text = root.clipboard_get().strip()
        if not clipboard_text.isdigit():
            messagebox.showerror("Error", "Clipboard must contain a numeric value")
            return
        entry_pos.delete(0, tk.END)
        entry_pos.insert(0, clipboard_text)
    except tk.TclError:
        messagebox.showerror("Error", "No data in clipboard")

btn_paste = tk.Button(root, text="Paste", command=paste_from_clipboard)
btn_paste.grid(row=5, column=1, padx=5, pady=5)

# כפתור בדיקה
btn_check_position = tk.Button(root, text="Check Position", command=check_position)
btn_check_position.grid(row=5, column=3, padx=5, pady=5)

pos_result_var = tk.StringVar()
pos_result_label = tk.Label(root, textvariable=pos_result_var, justify="left", fg="green")
pos_result_label.grid(row=6, column=0, columnspan=4, padx=5, pady=10)

# תווית טעינה
loading_label = tk.Label(root, text="", fg="blue", justify="left")
loading_label.grid(row=7, column=0, columnspan=4, padx=5, pady=5)


root.mainloop()

