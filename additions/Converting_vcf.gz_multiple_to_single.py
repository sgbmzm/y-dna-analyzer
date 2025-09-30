import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import gzip
import time
import threading

def load_vcf():
    global headers, sample_names, input_file
    input_file = filedialog.askopenfilename(
        title="Select multi-sample VCF.GZ file",
        filetypes=[("VCF GZ files", "*.vcf.gz")]
    )
    if not input_file:
        return

    try:
        with gzip.open(input_file, 'rt') as fin:
            for line in fin:
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    sample_names = headers[9:]
                    break
        if not sample_names:
            messagebox.showerror("Error", "No samples found in VCF file")
            return
    except Exception as e:
        messagebox.showerror("Error", f"Failed to read VCF: {e}")
        return

    sample_dropdown['values'] = sample_names
    sample_dropdown.current(0)
    extract_btn['state'] = 'normal'

def extract_sample_thread():
    thread = threading.Thread(target=extract_sample)
    thread.start()

def extract_sample():
    sample_name = sample_dropdown.get()
    if not sample_name:
        messagebox.showerror("Error", "Please select a sample")
        return

    output_file = filedialog.asksaveasfilename(
        title="Save single-sample VCF",
        defaultextension=".vcf.gz",
        filetypes=[("VCF GZ files", "*.vcf.gz")]
    )
    if not output_file:
        return

    try:
        sample_index = headers.index(sample_name)
        last_report = time.time()

        with gzip.open(input_file, 'rt') as fin, gzip.open(output_file, 'wt') as fout:
            for i, line in enumerate(fin, 1):
                if line.startswith('#CHROM'):
                    fout.write('\t'.join(headers[:9] + [sample_name]) + '\n')
                elif line.startswith('#'):
                    fout.write(line)
                else:
                    fields = line.strip().split('\t')
                    fout.write('\t'.join(fields[:9] + [fields[sample_index]]) + '\n')

                # דיווח כל חצי דקה
                if time.time() - last_report >= 30:
                    print(f"Processed {i} lines...")
                    last_report = time.time()

        print(f"Done! Single-sample VCF saved to: {output_file}")
        messagebox.showinfo("Done", f"Single-sample VCF saved to:\n{output_file}")

    except Exception as e:
        messagebox.showerror("Error", f"Failed to write VCF: {e}")

# GUI setup
root = tk.Tk()
root.title("Extract single-sample VCF")

tk.Button(root, text="Load VCF.GZ File", command=load_vcf).pack(pady=10)

sample_dropdown = ttk.Combobox(root, state="readonly")
sample_dropdown.pack(pady=5, padx=10, fill='x')

extract_btn = tk.Button(root, text="Extract Sample", command=extract_sample_thread, state='disabled')
extract_btn.pack(pady=10)

root.mainloop()

