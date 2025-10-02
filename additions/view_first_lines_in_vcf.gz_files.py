import tkinter as tk
from tkinter import filedialog, scrolledtext
import gzip
import zipfile

def read_vcf_file(file_path, y_only=False, max_lines=200):
    text_area.delete(1.0, tk.END)  # נקה את ה־Text לפני טעינה
    try:
        # קביעה איך לפתוח לפי סוג הקובץ
        if file_path.endswith(".gz"):
            f = gzip.open(file_path, 'rt', encoding="utf-8")
        elif file_path.endswith(".zip"):
            z = zipfile.ZipFile(file_path, 'r')
            # נניח שהקובץ הראשון בתוך ה־zip הוא ה־VCF
            inner_file = z.namelist()[0]
            f = z.open(inner_file, 'r')
            f = (line.decode("utf-8") for line in f)  # המרה ל־str
        else:
            f = open(file_path, 'r', encoding="utf-8")

        # קריאה בפועל
        count = 0
        for line in f:
            if y_only:
                if line.startswith('#'):  # שמור כותרות
                    text_area.insert(tk.END, line)
                    continue
                parts = line.strip().split('\t')
                if parts[0].lower().replace("chr", "") == "y":
                    text_area.insert(tk.END, line)
                    count += 1
            else:
                text_area.insert(tk.END, line)
                count += 1

            if count >= max_lines:
                break

        # סגירה (רק אם זה אובייקט קובץ רגיל)
        if hasattr(f, "close"):
            f.close()

    except Exception as e:
        text_area.insert(tk.END, f"Error reading file: {e}")

def open_vcf():
    file_path = filedialog.askopenfilename(
        title="Select VCF / VCF.GZ / ZIP file",
        filetypes=[("VCF files", "*.vcf *.vcf.gz *.zip")]
    )
    if file_path:
        read_vcf_file(file_path, y_only=False, max_lines=200)

def open_vcf_y_rows_only():
    file_path = filedialog.askopenfilename(
        title="Select VCF / VCF.GZ / ZIP file",
        filetypes=[("VCF files", "*.vcf *.vcf.gz *.zip")]
    )
    if file_path:
        read_vcf_file(file_path, y_only=True, max_lines=1000)


# GUI setup
root = tk.Tk()
root.title("VCF/VCF.GZ/ZIP Viewer")

btn1 = tk.Button(root, text="Open File (First 200 lines)", command=open_vcf)
btn1.pack(pady=10)

btn2 = tk.Button(root, text="Open File (Y-ROWS ONLY)", command=open_vcf_y_rows_only)
btn2.pack(pady=10)

text_area = scrolledtext.ScrolledText(root, width=120, height=30)
text_area.pack(padx=10, pady=10)

root.mainloop()
