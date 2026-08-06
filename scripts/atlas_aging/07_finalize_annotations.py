"""
07_finalize_annotations.py

依05(marker score)+06(差異表現)交叉比對後人工核對的結果,套用成最終標籤,
邏輯跟scripts/atlas_annotation/04_finalize_annotations.py完全一樣。

CORRECTIONS目前是空的——04當初的內容(cluster 17改成Skeletal muscle cells、
18/20改成T cells...等)是針對「3個月WT-only、25個cluster」那次分群結果
逐一人工核對出來的,cluster編號、甚至需要修正的cluster是否存在,都是那次
特定分群結果才成立的結論,不能直接套用到這裡(3個月+16個月合併、cluster
數量、每群的DE基因都會不一樣)。

跑完05+06之後,要重新打開cluster_marker_scores.xlsx(05的輸出)+
cluster_top_de_genes.xlsx(06的輸出),依同樣的方法人工核對:05判斷的
cell_type有沒有被06的top DE基因支持,對不上的cluster才需要在下面的
CORRECTIONS填入(最終標籤, 修改理由)。

Input:  results/atlas_aging/wt_aging_annotated.h5ad (05的輸出)
Output: results/atlas_aging/wt_aging_final_annotated.h5ad
        results/atlas_aging/cluster_annotation_corrections.xlsx (老師要看的對照表)
        results/atlas_aging/figures/umap_celltype_final.png
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 10,
})

# cluster(字串)-> (最終標籤, 修改理由)。跑完05+06、人工核對過後填入的結果
# (跟老師/AI一起核對過,05判斷vs06差異表現top20基因逐個cluster比對出來的)。
CORRECTIONS = {
    "9": ("Skeletal muscle cells",
          "top20基因是Ttn/Neb/Ryr1/Dmd/Myh4/Mybpc1/Myom1/Nrap/Tnnt3/Atp2a1等骨骼肌"
          "收縮蛋白基因的教科書等級組合,跟間皮細胞完全無關,05判Mesothelial cells證據不足"
          "(跟3個月版本cluster 17同一種錯誤模式:marker清單沒有專門對到骨骼肌的類型)"),
    "16": ("T cells",
           "top20基因出現Themis/Lef1/Prkcq,加上Skap1(T細胞受體訊號傳導蛋白),"
           "跟3個月版本cluster 18/20把ILCs改判T cells用的是同一組決定性基因,05判ILCs證據不足"),
    "18": ("T cells",
           "top基因出現Bcl11b(T細胞系別分化主控轉錄因子),加上Runx1/Tox/Ikzf2等"
           "T/淋巴球分化轉錄因子,05判ILCs證據不足"),
    "21": ("B cells",
           "top基因出現Pax5/Bank1/Bach2/Ebf1,皆為B細胞主控轉錄因子/高度特異marker,"
           "跟3個月版本cluster 19/22/23同一種錯誤模式,05判Endothelial cells錯誤"),
}


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    input_path = os.path.join(out_dir, "wt_aging_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Input: {adata.shape[0]} cells x {adata.shape[1]} genes, "
          f"{adata.obs['leiden'].nunique()} Leiden clusters")

    def apply_correction(row_leiden, row_celltype):
        leiden_str = str(row_leiden)
        if leiden_str in CORRECTIONS:
            return CORRECTIONS[leiden_str][0]
        return row_celltype

    adata.obs["cell_type_final"] = [
        apply_correction(l, c) for l, c in zip(adata.obs["leiden"], adata.obs["cell_type"])
    ]

    print("\nFinal cell_type distribution:")
    print(adata.obs["cell_type_final"].value_counts())

    cluster_summary = (
        adata.obs.drop_duplicates("leiden")[["leiden", "cell_type"]]
        .rename(columns={"cell_type": "cell_type_marker_based"})
        .set_index("leiden")
    )
    cluster_sizes = adata.obs["leiden"].value_counts()
    cluster_summary.insert(0, "n_cells", cluster_sizes.reindex(cluster_summary.index))
    cluster_summary["cell_type_final"] = [
        CORRECTIONS.get(str(l), (cluster_summary.loc[l, "cell_type_marker_based"], ""))[0]
        for l in cluster_summary.index
    ]
    cluster_summary["reason_for_change"] = [
        CORRECTIONS.get(str(l), (None, ""))[1] for l in cluster_summary.index
    ]
    cluster_summary = cluster_summary.sort_index(key=lambda idx: idx.astype(int))

    print("\nCorrection summary:")
    print(cluster_summary.to_string())

    excel_path = os.path.join(out_dir, "cluster_annotation_corrections.xlsx")
    cluster_summary.to_excel(excel_path)
    print(f"\nSaved: {excel_path}")

    fig = sc.pl.umap(
        adata, color="cell_type_final", show=False, return_fig=True,
        size=8, frameon=False, title="3mo+16mo WT atlas: cell type (final, manually reviewed)",
        legend_fontsize=9,
    )
    fig_path = os.path.join(figures_dir, "umap_celltype_final.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    output_path = os.path.join(out_dir, "wt_aging_final_annotated.h5ad")
    adata.write_h5ad(output_path)
    print(f"\nSaved final annotated atlas to: {output_path}")


if __name__ == "__main__":
    main()