from yclade import tree, snps
import networkx as nx
import csv

def get_clade_and_descendants_lists(
    tree_data, snp_name: str, include_descendants: bool = True, merge_snps: bool = False
):
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

    return branch_names, snps_list


# שימוש
tree_data = tree.get_yfull_tree_data(version=None, data_dir=None)

# דוגמה: כל הענפים והצאצאים + מאחד את כל ה-SNPים לרשימה אחת
branches_all, variants_all = get_clade_and_descendants_lists(
    tree_data, "L243+", include_descendants=True, merge_snps=True
)

print("ענפים:", branches_all)
print("SNPים מאוחדים:", variants_all)


