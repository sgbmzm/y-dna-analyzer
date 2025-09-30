import requests
from bs4 import BeautifulSoup
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

main_url = "https://jewishdna.net/index.html"
base_group_url = "https://jewishdna.net/{}.html"

headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/139.0.0.0 Safari/537.36"
}

# הורדת העמוד הראשי
r = requests.get(main_url, headers=headers)
r.raise_for_status()
soup = BeautifulSoup(r.text, "html.parser")

# חילוץ קבוצות מהעמוד הראשי
table_rows = soup.find_all("tr")
groups = []
for tr in table_rows:
    tds = tr.find_all("td")
    if len(tds) < 5:
        continue
    group_id = tds[0].text.strip()
    if not group_id.startswith("AB-"):
        continue  # סינון רק קבוצות חוקיות
    group_name = tds[1].text.strip()
    communities = [td.text.strip() for td in tds[6:15] if td.text.strip()]
    communities_str = ",".join(communities)
    groups.append({
        "AB-Group": group_id,
        "Communities": communities_str,
        "Group": group_name
    })

total_groups = len(groups)
print(f"Found {total_groups} valid groups.")

# פונקציה להורדת דף קבוצה והוצאת Final SNP
def fetch_final_snp(group, idx=None):
    group_url = base_group_url.format(group["AB-Group"])
    try:
        r = requests.get(group_url, headers=headers)
        r.raise_for_status()
        soup = BeautifulSoup(r.text, "html.parser")

        final_snp = None

        # מחפשים את כל הקישורים שמובילים ל-YFull tree
        snp_links = [a for a in soup.find_all("a") if "yfull.com/tree/" in a.get("href", "")]
        if snp_links:
            final_snp = snp_links[-1].text.strip()  # האחרון ברשימה

        group["Final SNP"] = final_snp

    except Exception as e:
        print(f"Error fetching {group['AB-Group']}: {e}")
        group["Final SNP"] = None

    if idx is not None:
        print(f"[{idx+1}/{total_groups}] Fetched {group['AB-Group']}, Final SNP: {group['Final SNP']}")
    return group

# הורדה מקבילית עם ThreadPoolExecutor
results = []
with ThreadPoolExecutor(max_workers=10) as executor:
    futures = {executor.submit(fetch_final_snp, group, idx): group for idx, group in enumerate(groups)}
    for future in as_completed(futures):
        results.append(future.result())

# שמירה ל-CSV
df = pd.DataFrame(results)
df.to_csv("jewishdna_groups_final.csv", index=False,  na_rep="None")
print("Done! CSV saved as 'jewishdna_groups_final.csv'.")
