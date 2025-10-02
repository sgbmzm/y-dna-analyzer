# יבוא הספריות הדרושות לתוכנה
import tkinter as tk
from tkinter import *
from tkinter import filedialog, messagebox, ttk
from platformdirs import user_data_dir
import traceback # ללכוד את כל השגיאות באופן מפורט 
import csv
import os
import re
import ast
import gzip
import zipfile
import io
import platform
import webbrowser
import subprocess
from urllib.request import urlopen
from pathlib import Path
import json
#import yclade
from yclade import tree, snps, find, const
import networkx as nx

# משתנה ששומר את תאריך הגרסה של התוכנה עבור הדפסה בשורת הכותרת של התוכנה
yda_version_date = "03 Oct 2025"

# משתנה מאוד חשוב שקובע האם מערכת ההפעלה הנוכחית היא ווינדוס כי אם היא לא אז אי אפשר לעשות חלק מהפעולות
is_windows = platform.system() == "Windows"

##########################################################################
# תוכן המידע על התוכנה
INFORMATION = '''
This project is hosted here:
https://github.com/sgbmzm/y-dna-analyzer/tree/main

This software analyzes raw DNA files.
Such as files from the testing companies Ancestry, 23andME, MyHeritage, FTDNA, as well as BIG-Y VCF files or whole genome WGS.
From the above files, only the lines dealing with the Y chromosome are read. Everything else is skipped.

This software has two modes of operation:
1. With raw DNA file: Used to find the appropriate branch on the yfull tree for the raw DNA file, and also used to manually check whether a particular variant is positive or negative in the subject using the search box, by searching by SNP name or by genomic position.
The raw DNA files can be loaded into the software as they come from the testing company, even without extraction.

2. Without raw DNA file: Used to search for information on a specific variant on the Y chromosome, using the search box, by searching by SNP name or by genomic position, and also to find the appropriate branch on the yfull tree for the specific variant typed in the search box. This mode is only possible when no user DNA file is loaded.
All modes of operation require loading an appropriate reference file, HG19 or HG38. The appropriate reference for a raw DNA file is usually identified automatically. If not automatically detected, and in the case of no raw DNA file, the user must select which reference to load.

Note: The genomic position is different between the two references for each variant.

The HG38 reference file is updated every day with new variant names.
You can re-download and update all files needed for the software, including the references, through the menu.

The Yclade column also contains information for the "Avotaynu" branch, i.e. the AB where the corresponding branch is located, or information about ABs that are below the corresponding branch in the case that the corresponding branch is an ancient branch.
If (A) appears next to AB, this means that there are Ashkenazim in this branch (but maybe others too).

The information about AB is not always accurate for two reasons:
1. Only for some of the ABs is it known which branch they are on in the tree.
2. Sometimes there is an error in the database for the corresponding branch. I try to check and update the information so that it is more accurate. The database can be updated - along with all the other files - through the menu.
All information about the ABs is here:
https://jewishdna.net/index.html

'''

# פונקצייה להצגת המידע על התוכנה
def show_information():
    # יצירת חלון חדש מעל root
    info_win = tk.Toplevel(root)
    info_win.title("Information")
    info_win.geometry("500x400")

    # יצירת Frame פנימי שיחזיק את הטקסט וה־Scrollbar
    frame = ttk.Frame(info_win)
    frame.pack(fill="both", expand=True, padx=10, pady=10)

    # יצירת Scrollbar אנכי
    scrollbar = ttk.Scrollbar(frame, orient="vertical")
    scrollbar.pack(side="right", fill="y")

    # תיבת טקסט לקריאה בלבד
    text_box = tk.Text(frame, wrap="word", height=20, width=60, yscrollcommand=scrollbar.set)
    text_box.pack(side="left", fill="both", expand=True)

    # חיבור ה־Scrollbar לתיבה
    scrollbar.config(command=text_box.yview)

    # הוספת המידע
    text_box.insert("1.0", INFORMATION)

    # נעילת התיבה לעריכה, עדיין מאפשרת העתקה
    text_box.config(state="disabled")


#####################################################################################

# נתיב לתיקיית הנתונים של קבצים הדרושים לתוכנה
yda_dir_path = user_data_dir(appname="y_dna_analyzer", appauthor=False, roaming=False)
# יצירת התיקייה אם לא קיימת ולא תהיה שגיאה אם כן קיימת
os.makedirs(yda_dir_path, exist_ok=True)

#######################################################################################################

# פונקצייה שמחזירה את מספר הגרסה העדכנית ביותר של עץ וויפול
def get_latest_yfull_tree_version() -> str:
    url = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_version.txt"
    try:
        with urlopen(url) as resp:
            return resp.read().decode("utf-8").strip()
    except Exception as e:
        tb_str = traceback.format_exc() # תופס את כל הפרטים אודות השגיאה ולא רק את נוסח השגיאה עצמה שזה e
        print(tb_str)
        print(f"get_latest_yfull_tree_version Eror: {tb_str}")
        return None  # ערך ברירת מחדל

######################################################################################################
                    # איזור הטיפול בקבצים הדרושים לתוכנה
######################################################################################################
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


# פונקצייה לעדכון הקבצים החשובים הדרושים לפעולת התוכנה
def update_required_files():
    # שאלה למשתמש האם מעוניין להוריד את הקבצים
    proceed = messagebox.askyesno(
        "Download Files",
        "Are you sure you want to download/update required files?"
    )
    if not proceed:
        messagebox.showinfo("Download canceled", "The download was canceled.")
        return
    
    # תיקיית שמירה
    global yda_dir_path
    save_dir = Path(yda_dir_path)
    save_dir.mkdir(parents=True, exist_ok=True)

    # רשימת קבצים להורדה
    files_to_download = [
        ("https://raw.githubusercontent.com/sgbmzm/y-dna-analyzer/main/pyinstaller.txt", "pyinstaller.txt"),
        ("https://raw.githubusercontent.com/sgbmzm/y-dna-analyzer/main/y-dna-analyzer.py", "y-dna-analyzer.py"),
        
        ("https://raw.githubusercontent.com/sgbmzm/y-dna-analyzer/main/required_files/ab_groups_snp.csv", "ab_groups_snp.csv"),
        ("https://raw.githubusercontent.com/sgbmzm/y-dna-analyzer/main/required_files/yda_yfull_tree_version.txt", "yda_yfull_tree_version.txt"),
        ("https://raw.githubusercontent.com/sgbmzm/y-dna-analyzer/main/required_files/yda_yfull_tree.json", "yda_yfull_tree.json"),
        ("https://raw.githubusercontent.com/sgbmzm/y-dna-analyzer/main/required_files/Msnps_hg19.vcf.gz", "Msnps_hg19.vcf.gz"),
        ("https://ybrowse.org/gbrowse2/gff/snps_hg38.vcf.gz", "snps_hg38.vcf.gz"),
    ]
    '''
    # גרסת עץ YFull הכי עדכנית והוספתה לרשימת ההורדות
    version_url = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_version.txt"
    with urlopen(version_url) as resp:
        version = resp.read().decode("utf-8").strip()
    yfull_tree_url = f"https://github.com/YFullTeam/YTree/raw/refs/heads/master/ytree/tree_{version}.zip"
    files_to_download.append((yfull_tree_url, f"tree_{version}.zip"))
    '''
    # ממשג גרפי להצגת התקדמות ההורדה
    global root
    progress_bar = ttk.Progressbar(root, length=200, maximum=len(files_to_download))
    progress_bar.grid(row=4, column=2, padx=5, pady=5)
    root.update()
    
    # הורדה ושמירה
    failed_files = []
    for i, (url, filename) in enumerate(files_to_download, start=1):
        save_path = save_dir / filename
        try:
            with urlopen(url) as response:
                data = response.read()
                # בודק שלא ייווצר קובץ אם הוא ריק ואין בו מידע כי לא הצליחו להוריד 
                if not data:
                    raise ValueError("Downloaded file is empty")
                with open(save_path, "wb") as out_file:
                    out_file.write(data)
        except Exception as e:
            tb_str = traceback.format_exc() # תופס את כל הפרטים אודות השגיאה ולא רק את נוסח השגיאה עצמה שזה e
            print(tb_str)
            failed_files.append(f"{filename} ({tb_str})")
        progress_bar['value'] = i
        root.update()

    # הסרת החלון לאחר סיום ההורדה
    progress_bar.destroy()

    # הודעה על הצלחה/כשל והחזרת טרו במקרה של ל הצלחה ופאלס במקרה של כשלון
    if not failed_files:
        messagebox.showinfo("Download required_files Complete", "All files were downloaded successfully.\nThe software will now close, and you need to restart it.")
        # הסרת החלון לאחר סיום ההורדה
        root.destroy()
        return True
    else:
        messagebox.showerror("Download required_files Failed", f"Failed to download:\n" + "\n".join(failed_files))
        return False

# פונקצייה חשובה מאוד שטוענת את כל הקבצים הדרושים לפעילות התוכנה ומחזירה את המשתנים הדרושים 
def get_required_files():
    # משתנה גלובלי מאוד חשוב ששומר את המידע האם הקבצים הדרושים לתוכנה קיימים
    is_required_files_exist = False
        
    # הגדרת הנתיבים לקבצים הדרושים לתוכנה בתיקייה הייעודית לקבצי התוכנה
    ab_groups_snp_path = os.path.join(yda_dir_path, "ab_groups_snp.csv")
    snps_hg38_path = os.path.join(yda_dir_path, "snps_hg38.vcf.gz")
    Msnps_hg19_path = os.path.join(yda_dir_path, "Msnps_hg19.vcf.gz")
    yda_tree_path = os.path.join(yda_dir_path, f"yda_yfull_tree.json")
    yda_tree_version = os.path.join(yda_dir_path, f"yda_yfull_tree_version.txt")

    # כל הקבצים הדרושים בתוך רשימה
    required_files = [ab_groups_snp_path, snps_hg38_path, Msnps_hg19_path, yda_tree_path, yda_tree_version]

    # אם כל הקבצים קיימים בתיקייה הייעודית אפשר להמשיך בתוכנה
    if all(os.path.exists(f) for f in required_files):
        is_required_files_exist = True
    # אם אפילו אחד מהם לא קיים מנסים במיקום השני שהוא בתיקייה שממנה רצה התוכנה
    else:
        ab_groups_snp_path = resource_path("ab_groups_snp.csv")  # כאן את שם הקובץ שלך
        snps_hg38_path = resource_path("snps_hg38.vcf.gz")
        Msnps_hg19_path = resource_path("Msnps_hg19.vcf.gz")
        yda_tree_path = resource_path("yda_yfull_tree.json")
        yda_tree_version = resource_path("yda_yfull_tree_version.txt")

        # כל הקבצים הדרושים בתוך רשימה
        required_files = [ab_groups_snp_path, snps_hg38_path, Msnps_hg19_path, yda_tree_path, yda_tree_version]

        # אם כולם נמצאים במיקום השני אז אפשר להמשיך בתוכנה
        if all(os.path.exists(f) for f in required_files):
            is_required_files_exist = True
        # אם אפילו אחד מהם לא נמצא אז מנסים לעדכן את הקבצים ולהוריד אותם מהאינטרנט
        # אם לא מצליחים אז יודעים שאין לתוכנה את הקבצים הדרושים
        else:
            # שואלים ואז קוראים לפונקציית הורדת הקבצים ואם היא מצליחה היא מחזירה טרו ואז הקבצים בתיקייה הייעודית להם
            if messagebox.askyesno("required files question", "The software is missing the required files. \nWould you like to download them now?") and update_required_files():
                ab_groups_snp_path = os.path.join(yda_dir_path, "ab_groups_snp.csv")
                snps_hg38_path = os.path.join(yda_dir_path, "snps_hg38.vcf.gz")
                Msnps_hg19_path = os.path.join(yda_dir_path, "Msnps_hg19.vcf.gz")
                yda_tree_path = os.path.join(yda_dir_path, f"yda_yfull_tree.json")
                yda_tree_version = os.path.join(yda_dir_path, f"yda_yfull_tree_version.txt")
                is_required_files_exist = True
            else:
                ab_groups_snp_path = None
                snps_hg38_path = None
                Msnps_hg19_path = None
                yda_tree_path = None
                yda_tree_version = None
                is_required_files_exist = False

    # הודעה למשתמש אם הקבצים הדרושים חסרים        
    if not is_required_files_exist:
        messagebox.showerror("Required files are missing", f"Required files are missing\nThe software is useless without these files\nPlease connect to the internet and download them from the menu")

    # אם יש את כל הקבצים הדרושים בודקים מה תאריך הגרסה של עץ ווייפול שבו משתמשים וטוענים אותו לשימוש בתוכנה זו
    if is_required_files_exist:
        with open(yda_tree_version, "r", encoding="utf-8") as f:
            string_version = f.readline().strip()  # קורא את השורה הראשונה ומסיר רווחים/סוף שורה
        ###########################################################    
        # זה בודק מה הגרסה האחרונה ואם אין חיבור לרשת או כל שגיאה אחרת זה מחזיר none
        latest_tree_version = get_latest_yfull_tree_version()
        # רק אם משתמשים בעץ שנמצא בתיקיית הנתונים (yda_dir_path in yda_tree_path) ויש גרסה עדכנית יותר באתר 
        if latest_tree_version and (yda_dir_path in yda_tree_path) and latest_tree_version != string_version:
            # מוסיפים סימני קריאה למספר הגרסה
            string_version += "( !!!)"
        ############################################################    
        # טוען את עץ וייפול לשימוש בתוכנה כולל מידע על הגרסה (file_path: Path, version: str | None = None) -> YTreeData
        yfull_tree_data = tree.yfull_tree_to_tree_data(Path(yda_tree_path), version=string_version)
    else:
        yfull_tree_data = None
            
    return is_required_files_exist, ab_groups_snp_path, snps_hg38_path, Msnps_hg19_path, yda_tree_path, yda_tree_version, yfull_tree_data

####################################################################################################
                # משתנים גלובליים עבור התוכנה
##########################################################################################################

last_positive_snp_string = "" # שומר את רשימת הווריאנטים החיוביים שבהם משתמשים לחישוב הקלייד ב run_calculate_clade
last_clades = [] # שומר את רשימת התוצאות שיוצאת מפונקציית run_calculate_clade לצורך שימוש עתידי
last_yfull_link = "" # מחזיק את כתובת הקישור לעץ yfull למיקום הקלייד הסופי
last_ftdna_link = "" # מחזיק את כתובת הקישור לעץ ftdna למיקום הקלייד הסופי
last_ab_link = None # מחזיק לינק לקבוצת אבותינו רק אם האדם נמצא בתוך קבוצת אבותינו מוגדר בפונקציית get_ab_from_clade

reference_positions_dict = {}  # מילון ששומר את נתוני הרפרנס לפי מפתח מיקום גנומי 
reference_names_dict = {}      # מילון רפרנס לפי מפתח שמות snp מחזיר מיקום גנומי מתאים

user_snps_dict = {} # מילון ששומר את כל הווריאנטים של המשתמש מתוך קובץ הדנ"א שלו כשהמפתח הוא המיקום הגנומי לפי הרפרנס המתאים או הנבחר

user_loaded = False # משתנה ששומר האם יש קובץ דנא של נבדק-משתמש שנטען לתוכנה ומוכן לשימוש
last_dna_file = "" # שומר את שם קובץ הדנ"א של המשתמש
last_dna_file_info = "" # שומר נתוני מידע שונים על קובץ הדנא של המשתמש שמתקבלים מפונקציית detect_headlines

reference_loaded = False # משתנה ששומר האם יש קובץ רפרנס שנטען לתוכנה ומוכן לשימוש
last_reference_file = "" # שומר את שם קובץ הרפרנס 
last_ref_file_info = "" # שומר נתוני מידע שונים על קובץ הרפרנס שמתקבלים מפונקציית detect_headlines

last_ab_data = None # מחזיק את הנתונים עבור כל קבוצות אבותינו ותתי הענפים שלהם

########################################################################################################

# פונקצייה לביטול טעינת קובץ הדנא של המשתמש וכל הנתונים שחושבו עליו
def reset_user():
    global last_clades, last_positive_snp_string, last_dna_file, last_dna_file_info, user_loaded    
    last_clades = []
    user_snps_dict = {}
    last_positive_snp_string = ""
    last_dna_file = ""
    user_loaded = False
    result_var.set("")
    btn_yfull.grid_forget()
    btn_ftdna.grid_forget()
    btn_ab.grid_forget()
    dna_loading_label.config(text="No DNA_file loaded", fg="red")
    yclade_label.config(text="Check SNP or load DNA-file", fg="red")
    user_result_var.set("")
    user_result_label.config(text="", bg="SystemButtonFace")
    btn_save_results.grid_forget()
    btn_unload_dna.grid_forget()
    last_dna_file_info = ""
    last_yfull_link = ""
    last_ftdna_link = ""
    last_ab_link = None
    

# פונקצייה לביטול טעינת קובץ הרפרנס    
def unload_ref():
    global reference_positions_dict, reference_names_dict, last_reference_file, last_ref_file_info, reference_loaded
    reference_positions_dict = {}      
    reference_names_dict = {}      
    last_reference_file = ""
    last_ref_file_info = ""
    reference_loaded = False
    reference_loading_label.config(text="No reference_file loaded", fg="red")
    btn_unload_ref.grid_forget()
    ref_result_var.set("")
    ref_result_label.config(text="", bg="SystemButtonFace")
    reset_user() # מוסיפים ביטול של טעינת קובץ המשתמש כי בלי רפרנס אי אפשר להציג נתוני משתמש
    
    
# clades = yclade.find_clade(last_positive_snp_string) זו הדרך הרגילה אבל היא ארוכה כי צריך כל פעם לטעון את העץ מחדש
# לכן עושים את זה כך ישר על yfull_tree_data שכבר טעננו בתחילת הקוד פעם אחת
# זה במקום yclade.find_clade 
def yda_find_clade(snp_string):
    snp_results = snps.parse_snp_results(snp_string)
    snp_results = snps.normalize_snp_results(snp_results=snp_results, snp_aliases=yfull_tree_data.snp_aliases)
    clades = find.get_ordered_clade_details(tree=yfull_tree_data, snps=snp_results)
    return clades



# פונקצייה מאוד חשובה שמחזירה את שמות כל ענפי הצאצאים שתחת ענף מסויים בוויפול
def get_clade_and_descendants_lists(tree_data, snp_name: str, include_descendants: bool = True, merge_snps: bool = False):
    
    """מחזיר שתי רשימות מקבילות:
       - כל שמות תתי הענפים (או רק הענף עצמו)
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
        # בדיקה אם רוצים דווקא קובץ זיפ של מייהירטייג ולא שיבחרו סתם זיפ אחר
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


# פונקצייה לבדיקה האם כתוב בקובץ באיזה רפרנס הוא משתמש באיזה תאריך נוצר ומי היוצר
def detect_headlines(file_path):
    
    ref = None
    creator = "unknown"
    creation_date = "unknown"
       
    #######################################################################
    if "haplocaller" in file_path:
        ref = "hg38"
        creator = "ftdna"
        creation_date = "unknown"
    
    elif is_gz_only(file_path): # אם זה קובץ שמסתיים רק ב gz ולא vcf.gz וכדומה
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
    # החזרת כל הנתונים הסופיים שנאספו כמילון
    return {"ref": ref, "creator": creator, "creation_date": creation_date}


# פונקצייה לטעינת קובץ הרפרנס
def load_reference(ref_path):
    global reference_loaded, reference_positions_dict, reference_names_dict, last_reference_file, last_ref_file_info
    
    if not ref_path:
        return
    
    reference_loading_label.config(text=f"Loading reference {os.path.basename(ref_path)} ...", fg="green")
    root.update()
    
    ref_file_info = detect_headlines(ref_path)
    
    # בדיקת או בחירת הרפרנס המתאים לקובץ הדנא של המשתמש
    ref_type = ref_file_info["ref"]
    
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
                
                # אם המיקום עדיין לא קיים במילון – צור רשימה חדשה שתכיל את כל הווריאנטים הקיימים עבור המיקום הגנומי הזה
                # צריך רשימה כי יש לפעמים שבמיקום גנומי אחד יש כמה ווריאנטים
                # אם כבר יש רשימה עבור המיקום הגומי הזה אז לא יוצרים רשימה אלא רק מוסיפים ערך לרשימה הקיימת
                if pos not in reference_positions_dict:
                    reference_positions_dict[pos] = []
                
                # הוסף את הווריאנט לרשימה של אותו מיקום
                reference_positions_dict[pos].append({
                    "pos": pos,
                    "name": snp_names,
                    "ref": ref,
                    "alt": alt
                })
                
                # לפצל את השמות לפי פסיק ולשמור כל שם בנפרד כמפתח במילון נוסף שמאונדקס לפי שמות
                # לשמות לא צריך רשימה כי כל שם מוביל בהכרח למיקום גנומי יחיד.
                if snp_names != ".":
                    for snp_name in snp_names.split(","):
                        #reference_names_dict[snp_name] = {"pos": pos, "name": snp_name, "ref": ref, "alt": alt}
                        reference_names_dict[snp_name] = pos
                        
        reference_loaded = True
        reference_loading_label.config(fg="blue", text=f"Reference loaded: \nname: {os.path.basename(ref_path)}  \n Y-SNPs: {len(reference_positions_dict)}  \ntype: {ref_type} \ndate: {ref_file_info['creation_date']}")
        last_reference_file = ref_path
        last_ref_file_info = ref_file_info
        btn_unload_ref.grid(row=1, column=0, padx=5, pady=5)
        
        
    except Exception as e:
        tb_str = traceback.format_exc() # תופס את כל הפרטים אודות השגיאה ולא רק את נוסח השגיאה עצמה שזה e
        print(tb_str)
        messagebox.showerror("load_reference Error", f"Failed to load reference: {tb_str}")
        #reference_loading_label.config(text="")

# המיקומים הללו בקובץ של מייהירטייג תמיד חיובי לכולם כולל נשים ולכן הם לא נכונים ויטופלו באופן שונה
myheritage_incorrect_variants_positions = [
    2679100, 3062994, 3105747, 3436239, 4179056, 4388629, 6739772, 6753387, 6892233, 7100362,
    9986197, 10023926, 10031793, 10032021, 13312756, 13450642, 13476461, 13484571, 14001289,
    14249991, 14491706, 14515040, 15425723, 17334694, 17361054, 17467526, 19179335, 28485084
    ]

# פונקצייה לטעינת קובץ הדנא של המשתמש מסוגים שונים ומחברות בדיקה שונות     
def load_user_dna_file():
    
    # כל פעם כשבוחרים קובץ דנא חדש הכל מתאפס ומחושב מהתחלה
    reset_user()
    
    # בחירת קובץ הנא של המשתמש
    file_path = filedialog.askopenfilename(
        title="Select DNA file",
        filetypes=[("DNA files", "*.txt *.csv *.gz *.zip *.vcf"), ("All files", "*.*")]
    )
    if not file_path:
        return
    
    # בדיקה האם בקובץ יש מידע על הרפרנס שבו צריכים להשתמש עבורו ומידע נוסף
    global last_dna_file_info
    last_dna_file_info = detect_headlines(file_path)
    
    # בדיקת או בחירת הרפרנס המתאים לקובץ הדנא של המשתמש
    ref_auto_detect = last_dna_file_info["ref"]
    
    ref_path = Msnps_hg19_path if ref_auto_detect == "hg19" else snps_hg38_path if ref_auto_detect == "hg38" else None #@@@@@@@@@
    
    if not ref_auto_detect:
        choice = messagebox.askyesnocancel("Reference not Autodetected", "Autodetected Reference faild\nChoose hg19, hg38, or None.\n\nYes = hg38, No = hg19, Cancel = None")
        ref_path = snps_hg38_path if choice else Msnps_hg19_path if (choice == False) else None #@@@@@@@@@
    
    global last_reference_file, reference_loaded
    
        
    if not last_reference_file == ref_path:
        reference_loaded = False
    
    if not reference_loaded:
        dna_loading_label.config(text="Waiting for Loading reference for this file ...", fg="green")
        yclade_label.config(text="Waiting for Loading reference & DNA-file ...", fg="green")
        # טעינת הרפרנס הנבחר    
        load_reference(ref_path)
    
    # הצהרה על משתנים גלובליים הדרושים להלן
    global last_clades, last_positive_snp_string, last_dna_file, user_snps_dict, user_loaded
    
    # אם לא בחרו קובץ רפרנס מאפסים הכל ולא ממשיכים
    if not reference_loaded:
        reset_user()
        return

    last_clades = []
    last_positive_snp_string = ""
    last_dna_file = file_path

    dna_loading_label.config(text=f"Loading {os.path.basename(file_path)} ...")
    root.update()

    positive_snps = []
    user_snps_dict = {}
          
    ##########################################################################################################
    # זה מטפל במקרה מיוחד של ביג Y של פטדנא שהוא בזיפ רגיל אבל הוא VCF
    is_ftdna_big_y_vcf = file_path.endswith(".zip") and get_first_file_name_in_zip(file_path) and get_first_file_name_in_zip(file_path).endswith(".vcf")
    # חייבים לדעת אם זה vcf או לא כי מבנה העמודות ב vcf שונה מהקבצים של החברות הרגילות
    is_vcf_file = file_path.endswith(".vcf") or file_path.endswith(".vcf.gz") or is_ftdna_big_y_vcf
    ##########################################################################################################
    
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
                        is_positive = "Yes"
                        alleles_str = alt
                    else:
                        is_positive = "No"
                        alleles_str = ref
                
                else: # זה לא ויסיאף אלא קובץ שרובו אוטוזומלי מהחברות הקטנות כמו מייהירטייג אנססטרי וכו ומכיל מעט מידע על כרומזום Y
                    # השדות הראשונים תמיד הם: rsid, chrom, pos, alleles
                    rsid, chrom, pos_str, alleles_str = parts[:4]
                    ad_str = ""
                    is_positive = "?"
                
                # עושה אלל אחד גדול גם היה קטן וגם אם היו שניים יחד
                allele_str = alleles_str.upper().strip()[0] if alleles_str else ""
            
                # נרמול שם הכרומוזום ל־Y
                chrom = chrom.replace("chr", "").upper()
                if chrom not in ("Y", "24"):
                    continue
                
                # המיקומים הללו בקובץ של מייהירטייג תמיד חיובי לכולם כולל נשים ולכן להוסיף אחריהם סימני שאלה וזה גם גורם שחישוב ענף וויפול לא יבוצע על פיהם
                if last_dna_file_info["creator"] == "myheritage" and int(pos_str) in myheritage_incorrect_variants_positions:
                    allele_str += "???" # או פשוט לדלג על השורה באמצעות: continue אבל אז מאבדים מידע
                 
                # הוספת השורה למילון המשתמש לאחר שוודאנו שמובר בשורה של Y
                user_snps_dict[int(pos_str)] = {"chrom": chrom, "pos_str": pos_str, "ref_type": ref_auto_detect, "snp_name": "?", "allele": allele_str, "is_positive": is_positive, "ad-R/A": ad_str}
                
                # בדיקה מול הרפרנס
                if not pos_str.isdigit():
                    continue
                
                # זה בודק האם יש מפתח מיקום גנומי כזה במילון הרפרנס ומחזיר רשימה של מילון אחד או יותר עבור מיקום גנומי זה
                # לפעמים יש למיקום גנומי אחד כמה ווריאנטים ואז מקבלים רשימה עם כמה איברים
                # אם יש רק ווריאנט אחד למיקום גנומי זה מקבלים רשימה עם איבר אחד
                ref_infos = reference_positions_dict.get(int(pos_str))
                
                # אם זה מחזיר none אז אין מידע ברפרנס על המיקום הגנומי הזה ולכן פשוט מדלגים עליו כי אי אפשר לבדוק האם חיובי או שלילי ומה שם הווריאנט
                if not ref_infos:
                    continue
                
                # עוברים על כל הווראינטים שיש ברפרנס למיקום הגנומי הנוכחי
                # רק ווריאנטים שעוברים את הביקורת החיובית כאן יכולים להיות מחושבים לחישוב הענף בווייפול
                # הגדרת Yes או No עבור is_positive חיונית בעיקר עבור קבצי דנ"א שאינם vcf, שאצלם בלי בדיקה זו אי אפשר לדעת האם חיוביים או שליליים
                for ref_info in ref_infos:
                    if allele_str == ref_info["alt"]:
                        s = f"{ref_info['name']}+"
                        positive_snps.append(s)
                        user_snps_dict[int(pos_str)]["is_positive"] = "Yes"   # או "כן" / "+" או כל מה שאתה רוצה
                        user_snps_dict[int(pos_str)]["snp_name"] = ref_info['name']
                    elif allele_str == ref_info["ref"] or allele_str == ".":
                        user_snps_dict[int(pos_str)]["is_positive"] = "No"
                        user_snps_dict[int(pos_str)]["snp_name"] = ref_info['name']

        last_positive_snp_string = ", ".join(positive_snps)
        #print("Final positive SNPs:", last_positive_snp_string)

        dna_loading_label.config(fg="blue", text=f"DNA file loaded: \nname: {os.path.basename(file_path)} \n{len(user_snps_dict)} total Y-rows  \n{len(positive_snps)} Positive Y-SNPs in DNA_file  \nref type: {ref_auto_detect}")
        user_loaded = True
        
        # אם יש פחות מחמישים ווריאנטים חיוביים זה אומר שיש בעיה ואין מה לחשב את ענף וייפול המתאים
        if len(positive_snps) >= 50:
            run_calculate_clade()
        else:
            yclade_label.config(text="female / incorrect reference \n(Too little Y positive variants)", fg="red")
        
        # הצבת כפתור ביטול הטעינה    
        btn_unload_dna.grid(row=1, column=4, padx=5, pady=5)

    except Exception as e:
        tb_str = traceback.format_exc() # תופס את כל הפרטים אודות השגיאה ולא רק את נוסח השגיאה עצמה שזה e
        print(tb_str)
        messagebox.showerror("load_user_dna_file Error", f"Failed reading file: {tb_str}")
        dna_loading_label.config(text="")
        yclade_label.config(text="Check SNP or load DNA-file", fg="red")


# פונקציה שטוענת מקובץ את הנתונים על קבוצות אבותינו וענפי הצאצאים שלהם ויכולה גם לשמור לקובץ
def get_ab_data():
    
    global last_ab_data, yfull_tree_data, ab_groups_snp_path
    
    if last_ab_data: #  אם כבר הנתונים נטענו אין צורך לטעון אותם שוב מהקובץ
        return
    
    ab_data = []
    
    # טעינת כל שורה בקובץ לתוך מילון
    with open(ab_groups_snp_path, newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')  # אם הקובץ מופרד בפסיקים  
        for row in reader:
            # כל שורה היא מילון עם שמות העמודות כמפתחות
            ab_data.append(row)
    
    # מעבר על כל השורות והוספת נתונים על תתי הענפים של כל ענף במידה וידוע ה- SNP המגדיר               
    # תמיד בודקים מחדש מהם תתי הענפים כי לפעמים יש עדכונים בעץ שמשפיעים על תתי הענפים
    for row in ab_data:     
        if row['Final SNP'] not in (None, '', 'None'): # בהמשך אולי להוסיף and row['snp_verified'] not in (None, '', 'None')
            row['sub_clades'] = get_clade_and_descendants_lists(yfull_tree_data, row['Final SNP'])
        
        '''
        # אם משתמשים בקובץ שיש בו מידע על תתי הענפים ורוצים להשתמש בהם במקום לחשב מחדש
        # לא טוב להשתמש בזה אז זה מבוטל ונמצא רק לשמירה
        # עדיף להשתמש באפשרות הקודמת ולבדוק כל פעם מחדש מהם תתי הענפים כי לפעמים יש עדכונים בעץ
        if 'sub_clades' in row and row['sub_clades'] not in (None, '', 'None'):
            row['sub_clades'] = ast.literal_eval(row['sub_clades'])
        '''
        
    # מיון לפי המספר שב־AB-Group
    ab_data.sort(key=lambda r: int(r["AB-Group"].split("-")[1]))
          
    # הגדרת המשתנה הגלובלי שיחזיק את כל המידע כדי שלא נטצרך לטעון כל פעם מחדש    
    last_ab_data = ab_data
    
    # האם לכתוב את הנתונים לקובץ. שימושי כדי לעשות קובץ עם תתי ענפים מעודכנים בקובץ
    write = False
    if write and messagebox.askyesno("Warning", "Do you want to write File?"):
        # כתיבה חזרה (דורכת על הקובץ המקורי)
        fieldnames = [k for k in ab_data[0].keys() if k != "sub_clades"] # שמות העמודות בלי עמודת תתי הקבוצות
        #fieldnames = [k for k in ab_data[0].keys()] # כל העמודות
        with open(ab_groups_snp_path, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')
            writer.writeheader()
            for row in ab_data:
                row_to_write = {k: v for k, v in row.items() if k in fieldnames} # מסנן רק את העמודות שהוגדרו בכותרות
                writer.writerow(row_to_write)


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
    global last_ab_data, yfull_tree_data, last_ab_link
    
    # איפוס המשתנש שזוכר את הלינק לקבוצת אבותינו הקודמת
    last_ab_link = None
    
    # אם רוצים מסניפ אז צריך קודם לדעת אל איזה ענף הסניפ הזה יושב בעץ ואז ממשיכים עם הענף המתאים
    if from_snp:
        #clade = yclade.find_clade(clade_snp)
        clade = yda_find_clade(snp_string)
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
        # התוספת and len(ab_clades_list) < 300 נועדה לפתור בעיה של ענפים שהוגדרו מידי למעלה בעץ ותוסר בהמשך לאחר שנוודא שכל הענפים מוגדרים נכון 
        if clade in ab_clades_clean and len(ab_clades_list) < 300:
            clade_found_ab = True
            clade_found_ab_rows.append(row)
            # רק אם זה בתוך ab יוצרים לינק. הבעיה היא שלפעמים זה בתוך כמה ab בגלל ענפים שמוגדרים לא נכון ואז יקבל את הלינק ל ab האחרון
            last_ab_link = f"https://jewishdna.net/{row['AB-Group']}.html"
        
        # בדיקה האם אחד מתתי הענפים של הענף הנבוקש נמצא בתוך אחת מקבוצות אבותינו
        elif any(a_clade.replace("*", "") in ab_clades_clean for a_clade in clade_sub_clades):
            sub_clades_found_ab = True
            sub_clades_found_ab_rows.append(row)
    
    # A אומרת שיש בקבוצה הזו אשכנזים. זה אוסף לרשימה של כל קבוצות אבותינו שנמצאו
    # בתוך ab אמור לצאת רק אחד ולא יותר אבל לפעמים יוצא יותר בגלל שיש ענפי ab שהענף שלהם מוגדר רחב ולא ספציפי ברשימה שבקובץ
    in_ab = ", ".join([f"{r['AB-Group']}{'(A)' if 'A' in r['Communities'] else ''}" for r in clade_found_ab_rows])
    above_ab = ", ".join([f"{r['AB-Group']}{'(A)' if 'A' in r['Communities'] else ''}" for r in sub_clades_found_ab_rows])
    ab_string = f'in AB: \n[{in_ab}]' if in_ab else f'above ABs: \n[{above_ab}]' if above_ab else ""
    # מחזירים את רשימת ענפי אבותינו שבתוך או שמעל    
    return ab_string

    
# פונקצייה שמחשבת על איזה ענף בעץ וויפול יושב הנבדק לפי הווריאנטים החיוביים שלו
# כברירת מחדל הענף הכי מדוייק הוא הענף הראשון במערך שמתקבל, שהוא בעל הציון סקור הכי גבוה
def run_calculate_clade(Final_clade_index = 0):
    
    # דבר ראשון קוראים לפונקצייה שמקימה את מערך המידע על קבוצות אבותינו
    get_ab_data()
    
    global last_positive_snp_string, last_clades, last_reference_file, last_ab_data, last_dna_file_info, last_ref_file_info
    result_var.set("") # תמיד לאפס קודם ולרוקן את הכיתוב הישן
    if not last_positive_snp_string:
        result_var.set("No matching SNPs were found in the file")
        yclade_label.config(text="Analysis finished - no matching SNPs found.")
        return

    try:
        yclade_label.config(text="Running clade calculation...", fg="green")
        root.update()

        # קריאה ל־yclade עם מחרוזת אחת
        try:                      
            clades = yda_find_clade(last_positive_snp_string)
                            
        except Exception as e:
            tb_str = traceback.format_exc() # תופס את כל הפרטים אודות השגיאה ולא רק את נוסח השגיאה עצמה שזה e
            print(tb_str)
            messagebox.showerror("yda_find_clade in run_calculate_clade Error", f"find_clade failed: {tb_str}")
            yclade_label.config(text="")
            return
        
        # שמירת הנתונים במשתנה הגלובלי
        last_clades = clades

        if not clades:
            result_var.set("No clades returned by yclade.")
            yclade_label.config(text="Analysis finished - no clades returned.")
            return
        
        # האחרון הוא בעל הציון הכי גבוה ולכן בדרך כלל הכי נכון
        Final_clade = clades[Final_clade_index]
        
        # קבלת מידע על קבוצת אבותינו שהענף נמצא בה או שהיא נמצאת תחת הענף באמצעות פונקציה
        ab_string = get_ab_from_clade(Final_clade.name)
        
        # בטיחות: קח שדות רק אם קיימים 
        name = getattr(Final_clade, "name", "Unknown")
        age_info = getattr(Final_clade, "age_info", None)
        tmrca = getattr(age_info, "most_recent_common_ancestor", "Unknown") if age_info else "Unknown"
        formed = getattr(age_info, "formed", "Unknown") if age_info else "Unknown"
        
        # הצהרה על משתנים גלובליים עבור הלינקים
        global last_yfull_link, last_ftdna_link, last_ab_link
        
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
            f"YFull predicted clade (used ref {last_ref_file_info['ref']}):\n"
            f"yfull tree version: {yfull_tree_data.version}\n"
            f" Name:    {name}\n"
            f" TMRCA:   {tmrca} ybp\n"
            f" FORMED:  {formed} ybp\n"
            f"{ab_string_multiline}"
        )
        
        #######################################################################################
        # הצבת כפתורים שרלוונטיים לתוצאות
        # לאחר שיש תוצאות הנחת כפתורי הלינקים לאילנות וייפול ופטדנא וכפתור שמירת התוצאות
        btn_yfull.grid(row=3, column=2, padx=5, pady=5)
        btn_ftdna.grid(row=4, column=2, padx=5, pady=5)
        if last_ab_link: # זה אומר שנמצא בתוך ab ולכן יש לינק. מוגדר בפונקציית get_ab_from_clade
            btn_ab.grid(row=5, column=2, padx=5, pady=5)
        else:
            btn_ab.grid_forget()
        btn_save_results.grid(row=6, column=2, padx=5, pady=5)
        #####################################################################################
        
        yclade_label.config(text="Analysis finished successfully.", fg="blue")
        
        
        # בדיקה בטוחה עם wrap-around
        next_index = (Final_clade_index + 1) % len(clades)
        if len(clades) >= 2 and clades[Final_clade_index].score == clades[next_index].score:
            yclade_label.config(text=more_result_warning, fg="red")
            if messagebox.askyesno("Yclade Warning", "There is more than one top-scoring clade.\nDo you want to see the next one?"):
                return run_calculate_clade(Final_clade_index=next_index)

        '''
        # במקרה שיש יותר מענף אחד בציון הכי גבוה נותנת אופצייה להתקדם למיקום הבא ברשימת התוצאות
        if len(clades) >= 2 and clades[Final_clade_index].score == clades[(Final_clade_index+1)].score:
            yclade_label.config(text=more_result_warning, fg="red")
            if messagebox.askyesno("Yclade Warning", "There is more than one top-scoring clade.\nDo you want to see the next one?"):
                return run_calculate_clade(Final_clade_index = (Final_clade_index+1) % len(clades))
        '''
        '''
        # זו צורה נוספת שבודקת רק את המיקום האחרון והבא אחריו, וההודעה קובץ כל פעם מהמיקום האחרון לזה שאחריו וחוזר חלילה    
        if len(clades) >= 2 and clades[0].score == clades[(1)].score:
            yclade_label.config(text=more_result_warning, fg="red")
            if messagebox.askyesno("Warning", "There is more than one top-scoring clade.\nDo you want to see the next one?"):
                return run_calculate_clade(Final_clade_index = 1)
        '''
        '''
        # זו עוד דרך לדעת שיש יותר מענף אחד שקיבל אותו סקור. לפי זה שהענף האחרון הוא לא צאצא של הענף שלפניו
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
        tb_str = traceback.format_exc() # תופס את כל הפרטים אודות השגיאה ולא רק את נוסח השגיאה עצמה שזה e
        print(tb_str)
        messagebox.showerror("run_calculate_clade Error", tb_str)
        
# פונקצייה לשמירת התוצאות המלאות של חישוב הענפים החיוביים מ run_calculate_clade
def save_clades_to_file():
    if not last_clades:
        messagebox.showerror("save_clades_to_file Error", "No Clades to save.")
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
        tb_str = traceback.format_exc() # תופס את כל הפרטים אודות השגיאה ולא רק את נוסח השגיאה עצמה שזה e
        print(tb_str)
        messagebox.showerror("save_clades_to_file Error", f"Failed to save file: {tb_str}")

# ------------------------
# פונקצייה לבדיקת שם סניפ או מיקום גנומי בנתוני המשתמש וברפרנס ומה הענף המתאים בוויפול ובאבותינו לווריאנט המבוקש     
def check_search_input(ref_search = True):
    
    # הצהרה על משתנים גלובליים.
    # צריך לזכור שreference_positions_dict מכיל רשימה (הכוללת מילון אחד או יותר) עבור כל מיקום גנומי כי לפעמים יש כמה ווריאנטים במיקום גנומי אחד
    # שאר המילונים מכילים ערך מילוני רגיל ולא רשימה
    global reference_positions_dict, reference_names_dict, user_snps_dict, user_loaded, reference_loaded, last_positive_snp_string
    
    # טוענים רפרנס לפי בחירת המשתמש אם עדיין לא נטען
    if not reference_loaded:
        choice = messagebox.askyesnocancel("Reference not Autodetected", "Autodetected Reference faild\nChoose hg19, hg38, or None.\n\nYes = hg38, No = hg19, Cancel = None")
        ref_path = snps_hg38_path if choice else Msnps_hg19_path if (choice == False) else None #@@@@@@@@@
        yclade_label.config(text="Waiting for Loading reference ...", fg="green")
        load_reference(ref_path) # טוען את קובץ הרפרנס שנבחר
        yclade_label.config(text="Check SNP or load DNA-file", fg="red")
        
    if not reference_loaded:
        return
    
    # לשכוח את הכפתורים שאולי הונחו עבור החיפוש הקודם. אבל אם כן טעון דנא של משתמש אז לא למחוק אותם כי הם עבור הדנא של המשתמש
    if not user_loaded:
        yclade_label.config(text="Check SNP or load DNA-file", fg="red")
        btn_yfull.grid_forget()
        btn_ftdna.grid_forget()
        btn_ab.grid_forget()
        btn_save_results.grid_forget()
    
    
    # לוקחים את מה שהמשתמש הקליד בתיבת החיפוש
    search_input = entry_search.get().strip().upper()[:10] # מסירים רווחים הופכים לאותיות גדולות וחותכים את כל מה שמעבר ל 10 אותיות או ספרות
    
    if not search_input:   # אם אין כלום בשדה החיפש
        return
    
    # אם הקליד מספר זה מיקום גנומי ובודקים ברפרנס ובנתוני המשתמש מה כתוב בשורה המתאימה למיקום גנומי זה
    elif search_input.isdigit():
        pos = int(search_input)
        fields_reference = reference_positions_dict.get(pos)
        fields_user = user_snps_dict.get(pos) if user_loaded else False
    # אם זה לא מספר אז מניחים שזה שם ווראינט ומחפשים אותו במילון של שמות הווראינטים שברפרנס
    else:  
        pos = reference_names_dict.get(search_input)
        fields_reference = reference_positions_dict.get(pos)
        fields_user = user_snps_dict.get(pos) if user_loaded else False
    
    # אם יש תוצאות מהרפרנס למיקום גנומי זה מציגים אותם בעמודת הרפרנס
    if fields_reference:
        # אם יש כמה ווריאנטים לאותו מיקום גנומי זה יציג כל אחד בשורה נפרדת
        fields_reference_for_result_var = "\n".join(str(item) for item in fields_reference)
        ref_result_var.set(fields_reference_for_result_var)
        ref_result_label.config(fg="green", bg="SystemButtonFace")
    else:
        ref_result_var.set(f"{search_input} not found in reference file")
        ref_result_label.config(fg="blue", bg="yellow")
    
    # אם יש תוצאות מהמשתמש למיקום גנומי זה מציגים בעמודת המשתמש
    if fields_user:
        user_result_var.set(fields_user)
        fg_for_user_result_label = "green" if fields_user["is_positive"] == "Yes" else "red"
        user_result_label.config(fg=fg_for_user_result_label, bg="SystemButtonFace")       
    else:
        msg = f"{search_input} not found in user DNA_file" if user_loaded else "user DNA_file_not_loaded"
        user_result_var.set(msg)
        user_result_label.config(fg="blue", bg="yellow")
        
    # רק במקרה שלא טענו קובץ דנא של נבדק אז מריצים את חישוב המיקום על עץ ווייפול עבור הווריאנט המבוקש כאילו שהוא חיובי כולל בדיקת קבוצת אבותינו המתאימה
    if not user_loaded and fields_reference:
        # במקרה שבחיפוש הוזן מספר אז מדובר במיקום גנומי וייתכן שיש כמה ווריאנטים במיקום גנומי זה
        # אבל אם הוזן משהו שאינו מספר אז מדובר בשם של ווריאנט וצריך לוודא שלוקחים רק את הווריאנט המבוקש
        positive_snps = fields_reference if search_input.isdigit() else [item for item in fields_reference if search_input in item['name']]
        # זה מטפל גם במקרים שעבור המיקום הגנומי שהוזן יש כמה ווריאנטים וזה יוצר סטרינג מתאים לפורמט של run_calculate_clade
        last_positive_snp_string = ", ".join(f"{item['name']}+" for item in positive_snps)
        # אחרי שהגדרנו את המשתנה הגלובלי last_positive_snp_string מריצים עליו את run_calculate_clade כדי לחשב את הענף על העץ
        run_calculate_clade()
        
        
# פונקצייה להדבקה מהלוח
def paste_from_clipboard():
    clipboard_text = root.clipboard_get().strip()[:15]
    entry_search.delete(0, tk.END)
    entry_search.insert(0, clipboard_text)
    
###############################################################################################################
                           # איזור הגדרת החלון הגראפי של התוכנה
###############################################################################################################

root = tk.Tk()
#root.attributes("-topmost", True) # מבטיח שהחלון יהיה מעל כל האחרים
# קביעת גודל התחלתי (רוחב x גובה)
#root.geometry("650x500")

# קביעת מינימום גודל
#root.minsize(500, 500)

# כותרת לחלון
root.title(f"Y-DNA-Analyzer | by Dr. simcha-gershon Bohrer (PhD.) | versin date: {yda_version_date}")

# מפרידים אנכיים בין העמודות
ttk.Separator(root, orient="vertical").grid(row=0, column=1, sticky="ns", padx=5, rowspan=20)
ttk.Separator(root, orient="vertical").grid(row=0, column=3, sticky="ns", padx=5, rowspan=20)

# כותרות לעמודות
tk.Label(root, text="Reference file", font="david 14 bold").grid(row=0, column=0, padx=75, pady=10)
tk.Label(root, text="Yclade", font="david 14 bold").grid(row=0, column=2, padx=75, pady=10)
tk.Label(root, text="User DNA file", font="david 14 bold").grid(row=0, column=4, padx=75, pady=10)

# כפתור לביטול טעינת קובץ הרפרנס # המיקום גריד שלו מתבצע בפונקציית טעינת הרפרנס
btn_unload_ref = tk.Button(root, text="unload ref file", command=unload_ref)

# כפתור לביטול טעינת קובץ דנא של המשתמש # המיקום גריד שלו מתבצע בפונקציית טעינת קובץ דנא של המשתמש
btn_unload_dna = tk.Button(root, text="unload dna file", command=reset_user)

# כפתור לבחירת קובץ הדנא של המשתמש
btn_csv = tk.Button(root, text="Choose \nUser RAW-DNA File \nvcf/vcf.gz/txt/csv/gz/zip", command=load_user_dna_file)
btn_csv.grid(row=2, column=4, rowspan=2)

# תווית מידע על קובץ הרפרנס
reference_loading_label = tk.Label(root, text="No reference-file loaded", fg="red")
reference_loading_label.grid(row=4, column=0, padx=5, pady=5, rowspan=3)

# תווית מידע על קובץ המשתמש
dna_loading_label = tk.Label(root, text="No DNA-file loaded", fg="red")
dna_loading_label.grid(row=4, column=4, padx=5, pady=5, rowspan=3)

####################################################################################################

# תווית מידע על הצלחת חישוב yclade
yclade_label = tk.Label(root, text="Check SNP or load DNA-file", anchor="w", fg="red")
yclade_label.grid(row=1, column=2, padx=5, pady=5)

# תווית מידע של תוצאת חישובי yclade
result_var = tk.StringVar()
result_label = tk.Label(root, textvariable=result_var, fg="green")
result_label.grid(row=2, column=2, padx=5, pady=10)

# כפתורי לינקיף לעצים של y-dna. ולאתר גואיש דנא, וכפתור לשמירת תוצאות מלאות של yclade. # הגריד שלהם נמצא בפונקציית חישוב הקלייד
btn_yfull = tk.Button(root, text="open clade in Yfull Tree", command= lambda: webbrowser.open_new(last_yfull_link))
btn_ftdna = tk.Button(root, text="open clade in FTDNA-Discover Tree", command= lambda: webbrowser.open_new(last_ftdna_link))
btn_ab = tk.Button(root, text="open AB in JewishDNA website", command= lambda: webbrowser.open_new(last_ab_link))
btn_save_results = tk.Button(root, text="Save Clades to TXT", command=save_clades_to_file)

##########################################################################################################3

# קו מפריד מאוזן וגם תווית כותרת עבור איזור בדיקת מיקום גנומי ידני
ttk.Separator(root, orient="horizontal").grid(row=10, column=2, sticky="ew", padx=5, pady=30)
tk.Label(root, text="Check Y-SNP", font="david 14 bold").grid(row=11, column=2, padx=5, pady=5)

# כפתור הדבקת תוצאות מהלוח אל תיבת החיפוש
btn_paste = tk.Button(root, text="Paste", command=paste_from_clipboard)
btn_paste.grid(row=12, column=2, padx=5, pady=5)

# תיבת חיפוש 
entry_search = tk.Entry(root, width=15)
entry_search.grid(row=13, column=2, padx=5, pady=5)

# כפתור בדיקה עבור הנתונים מתיבת החיפוש
btn_check = tk.Button(root, text="Check: SNP name / Genomic position", command=check_search_input) 
btn_check.grid(row=14, column=2, padx=5, pady=5)

# איזור לתוצאות החיפוש בנתוני הרפרנס
ref_result_var = tk.StringVar()
ref_result_label = tk.Label(root, textvariable=ref_result_var, fg="green")
ref_result_label.grid(row=14, column=0)

# איזור לתוצאות החיפוש בנתוני המשתמש
user_result_var = tk.StringVar()
user_result_label = tk.Label(root, textvariable=user_result_var, fg="green")
user_result_label.grid(row=14, column=4)

# תווית מידע
tk.Label(root, text="NOTE: Each reference has different positions").grid(row=15, column=2, padx=5, pady=5)

########################################################################################

# תפריט לחצנים נוספים לאפשרויות נוספות
mb=  Menubutton ( root, text= "Menu: Info & Options", relief=RAISED ,bg="gray87")
mb.grid(column=4, row=12)
mb.menu =  Menu ( mb, tearoff = 0 )
mb["menu"] =  mb.menu
mb.menu.add_command ( label= "Information", command= show_information)
mb.menu.add_command ( label= "open yda dir", command= lambda: subprocess.run(["explorer" if is_windows else "xdg-open", str(yda_dir_path)]))
mb.menu.add_command ( label= "Download/Update required files", command= update_required_files)

'''
# אם רוצים תפריט רגיל בחלק העליון של המסך
menubar = Menu(root)
root.config(menu=menubar)
options_menu = Menu(menubar, tearoff=0)
options_menu.add_command(label="open_yda_dir", command=lambda: subprocess.run(["explorer" if is_windows else "xdg-open", str(yda_dir_path)]))
options_menu.add_command(label="update_basic_files", command=update_required_files)
menubar.add_cascade(label="Help & Options", menu=options_menu)
'''
##########################################################################
# הבאת כל המשתנים של הקבצים הדרושים לתוכנה באמצעות פונקציית get_required_files() שהוגדרה למעלה
is_required_files_exist, ab_groups_snp_path, snps_hg38_path, Msnps_hg19_path, yda_tree_path, yda_tree_version, yfull_tree_data = get_required_files()

if not is_required_files_exist:
    # תווית מידע
    tk.Label(root, text="!!!!!!! Required files are missing !!!!!!!!\nThe software is useless without these files\nconnect to the internet and download them from the menu", fg="blue", bg="yellow").grid(row=4, column=2, padx=5, pady=5)

# זה מריץ כל הזמן את החלון הראשי שיהיה קיים תמיד
root.mainloop()


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





