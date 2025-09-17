import tkinter as tk
from tkinter import filedialog, scrolledtext
import gzip

def open_vcf():
    file_path = filedialog.askopenfilename(
        title="Select VCF.GZ file",
        filetypes=[("VCF GZ files", "*.vcf.gz")]
    )
    if not file_path:
        return

    text_area.delete(1.0, tk.END)  # נקה את ה־Text לפני טעינה

    try:
        with gzip.open(file_path, 'rt') as f:  # 'rt' = text mode
            for i, line in enumerate(f):
                text_area.insert(tk.END, line)
                if i >= 200:  # הדפס 200 שורות
                    break
    except Exception as e:
        text_area.insert(tk.END, f"Error reading file: {e}")

# GUI setup
root = tk.Tk()
root.title("VCF.GZ Viewer (First 200 lines)")

btn = tk.Button(root, text="Open VCF.GZ File", command=open_vcf)
btn.pack(pady=10)

text_area = scrolledtext.ScrolledText(root, width=120, height=30)
text_area.pack(padx=10, pady=10)

root.mainloop()

