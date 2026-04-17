import scanpy as sc
import os
import matplotlib.pyplot as plt
import numpy as np
# ------ Settings ---------

sc.settings.figdir = "0_explore/outputs/figures"
figdir = "0_explore/outputs/figures"
sc.set_figure_params(dpi=500,dpi_save=500)

# read  data

adata = sc.read_h5ad("data/atlas_light.h5ad")

# explore shape and metadata

adata.shape # (622529, 29923)
adata.obs_keys() 

# ['n_genes_by_counts',
#  'doublet_score',
#  'n_genes',
#  'condition',
#  'batch',
#  'sample_id',
#  'technology',
#  'study',
#  'leiden_0.40',
#  'cell_type_level0',
#  'cell_type_level1']


adata.obs["cell_type_level0"].value_counts()
# Goblet cells                       105146
# Epithelial cells-SI enterocytes    103506
# Epithelial cells-colon              83855
# Epithelial cells                    76828
# Fibroblasts                         49307
# T cells                             31796
# B cells                             30506
# Endothelial cells                   28256
# Prolif. cells                       28023
# Neuronal cells                      18433
# Smooth muscle cells                 15345
# Macrophages                         10762
# IgA plasma cells                    10463
# Paneth cells                         7923
# T cells-NK cells                     6460
# Tuft cells                           5251
# Enteroendocrine cells                4148
# ILC2                                 1826
# Mast cells                           1593
# Paneth-Goblet cells                  1253
# Neutrophils                          1075
# Epithelial cells?                     621
# Keratinocytes                         153



adata.obs["cell_type_level1"].value_counts()
# Epithelial cells-colon-1             50359
# Epithelial cells-SI enterocytes-1    49999
# Epithelial cells-SI enterocytes-2    46514
# Goblet cells-1                       34188
# Epithelial cells-colon-2             32622
# T cells                              31796
# Fibroblasts-1                        30923
# Prolif. cells                        28023
# Epithelial cells-2                   26047
# Goblet cells-2                       21919
# Goblet cells-3                       19422
# Epithelial cells-3                   17603
# B cells-1                            17503
# Epithelial cells-5                   16843
# Goblet cells-4                       16762
# Fibroblasts-2                        16722
# Neuronal cells-1                     14352
# Epithelial cells-4                   13309
# Endothelial cells-2                  12730
# Smooth muscle cells-1                11285
# Macrophages                          10762
# Endothelial cells-1                   9464
# IgA plasma cells-2                    9387
# Epithelial cells-SI enterocytes-3     6993
# Paneth cells-1                        6951
# B cells-3                             6942
# Goblet cells-5                        6635
# T cells-NK cells                      6460
# Endothelial cells-3                   6062
# B cells-2                             5428
# Tuft cells                            5251
# Goblet cells-6                        5059
# Enteroendocrine cells                 4148
# Neuronal cells-2                      4081
# Smooth muscle cells-2                 4060
# Epithelial cells-6                    2102
# ILC2                                  1826
# Mast cells                            1593
# Paneth-Goblet cells                   1253
# Goblet cells-7                        1161
# IgA plasma cells-1                    1076
# Neutrophils                           1075
# Paneth cells-2                         972
# Epithelial cells-1                     924
# Fibroblasts-4                          876
# Epithelial cells-colon-3               874
# B cells-5                              633
# Epithelial cells?                      621
# Fibroblasts-5                          395
# Fibroblasts-3                          391
# Keratinocytes                          153



adata.obs["study"].value_counts()
# Drokhlyansky    414307
# CRCDiet          84690
# Xu               55908
# Haber            26285
# Monarach         13790
# Parigi           10584
# Niec              9275
# Ayyaz             7690


adata.obs["sample_id"].value_counts()
# Drokhlyansky_SCP1038_mli                                  334723
# Drokhlyansky_SCP1038_msi                                   75242
# Xu_GSE124880_GSE124880_PP_LP_mm10                          55908
# ext_L24854_CD-AOM-DSS-Immune-single-cell                   19063
# ext_L24854_HFD-AOM-DSS-Immune-single-cell                  17130
# ext_L24854_LFD-AOM-DSS-Immune-single-cell                  15559
# ext_L24854_HFD-AOM-DSS-Epi_plus_DN-single-cell             14730
# Haber_GSE92332_GSE92332_SalmHelm_UMIcounts.txt              9610
# ext_L24854_LFD-AOM-DSS-Epi_plus_DN-single-cell              9578
# ext_L24854_CD-AOM-DSS-Epi_plus_DN-single-cell               8630
# Haber_GSE92332_GSE92332_atlas_UMIcounts.txt                 7057
# Parigi_GSE163638_GSM4983265_IEC-Stroma-control              6834
# Niec_GSE190037_GSM5712428_Colon_Lymphatics                  6307
# Haber_GSE92332_GSE92332_SalmonellaInfect_UMIcounts.txt      5016
# Haber_GSE92332_FAE_UMIcounts.txt                            4602
# Morarach_GSE149524_GSM4504450_P21_1                         4115
# Parigi_GSE163638_GSM4983266_IEC-Stroma-depleted             3750
# Morarach_GSE149524_GSM4504451_P21_2                         3424
# Morarach_GSE149524_GSM4504448_E15                           3384
# Niec_GSE190037_GSM5712427_SISC_Lymphatics                   2968
# Morarach_GSE149524_GSM4504449_E18                           2867
# Ayyaz_GSE123516_GSM3308718_C05                              2522
# Ayyaz_GSE123516_GSM3308720_C07                              2130
# Ayyaz_GSE123516_GSM3308717_C04                              1974
# Drokhlyansky_SCP1038_mli.neur                               1895
# Drokhlyansky_SCP1038_mli.glia                               1573
# Ayyaz_GSE123516_GSM3308719_C06                              1064
# Drokhlyansky_SCP1038_msi.neur                                467
# Drokhlyansky_SCP1038_msi.glia                                407
# plot annotation


adata.obsm_keys() # UMAP




# --------------- Plot UMAP with cell types --------------
sc.pl.umap(
    adata,
    color="cell_type_level0",
    save="_ct0.png")

sc.pl.umap(
adata,
color="cell_type_level1",
save="_ct1.png")


# --------- DGE ----------

# preprocess (not normalized?)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# --- ct0

sc.tl.rank_genes_groups(
    adata=adata,
    groupby="cell_type_level0",
    key_added="dge_ct0",
    n_genes=5
)

sc.pl.rank_genes_groups_dotplot(
    adata=adata,
    groupby="cell_type_level0",
    key="dge_ct0",
    n_genes=5,
    save="_ct0.png"
)




# --- ct1

sc.tl.rank_genes_groups(
    adata=adata,
    groupby="cell_type_level1",
    key_added="dge_ct1",
    n_genes=5
)

sc.pl.rank_genes_groups_dotplot(
    adata=adata,
    groupby="cell_type_level1",
    key="dge_ct1",
    n_genes=5,
    save="_ct1.png"
)

# ---- plot markers
marker_dict= { 
    "Macrophage"          : ["Adgre1", "Cd68", "Csf1r", "C1qa", "C1qb", "Cx3cr1","Spp1","Cd163","Trem2"],
    "Monocyte"            : ["Ly6c2", "Ccr2", "Cd14", "Cd16"],
    "Neutrophil"          : ["Ly6g", "Csf3r", "Cxcr2","S100a8", "S100a9"],
    "cDC1"                : ["Xcr1", "Clec9a", "Itgae", "Batf3"],
    "cDC2"                : ["Cd209a", "Clec10a", "Sirpa", "H2-Ab1"],
    "mRegDCs"             : ["Mreg","Ccr7","Cd274","Cd200"],
    "Mast"                : ["Kit", "Mcpt4", "Cpa3", "Fcer1a"],
    "T"                   : ["Cd8a", "Cd8b1", "Gzmb", "Gzma", "Gzmk", "Cd3d","Cd3e", "Ifng",
                                "Cd4", "Tbx21","Rorc","Foxp3","Gata3","Il10","Il17a",
                                "Tcrg","Tcrd","Trdc"],
    "Treg"                : ["Foxp3", "Il2ra", "Ctla4", "Ikzf2"],
    "ILC1_NK"             : ["Nkg7", "Klrb1c", "Ncr1", "Gzma", "Klrd1","Eomes","Il7r"],
    "ILC2"                : ["Gata3", "Il13", "Il5", "Areg"],
    "ILC3"                : ["Rorc", "Il22", "Il17a", "Kit"],
    "Stem"                : ["Lgr5", "Olfm4", "Ascl2", "Axin2"],
    "Colonocyte"          : ["Slc26a3", "Fabp1", "Car1", "Car4"],
    "Goblet"              : ["Muc2", "Tff3", "Agr2", "Spdef", "Clca1"],
    "Enteroendocrine"     : ["Chga", "Chgb", "Tac1", "Neurod1"],
    "Tuft"                : ["Dclk1", "Trpm5", "Il25"],
    "Proliferating"       : ["Mki67", "Top2a", "Pcna", "Stmn1"],
    "Fibroblast"          : ["Wnt2b","Col15a1","Ccn2","Pdgfra","Wnt5a","Bmp5","Bmp7","Adamdec1","Tagln","Gli1","Ncam1","Itga8","Pi16","Ackr4", "Col14a1", "Cd34", "Grem1","Jam2", "Pdgfra", "Col27a1", "Sfrp1", "Ngf", "Gdf10","C3","Thbs1","Timp1","Col1a2","Igfbp6"],
    "Smooth_muscle"       : ["Acta2", "Tagln", "Des", "Cnn1","Cspg2","Myh11"],
    "Pericyte"            : ["Rgs5", "Pdgfrb", "Notch3"],                                                    
    "Epithelial"          : ["Epcam", "Cdh1", "Krt8", "Krt18", "Krt19"],             
    "Goblet"              : ["Muc2", "Tff3", "Agr2", "Spdef", "Clca1"],                  
    "T_cell"              : ["Cd3e", "Cd3d", "Trac", "Cd8a", "Cd4"],
    "B_cell"              : ["Cd19", "Ms4a1", "Cd79a", "Cd79b", "Pax5"],                 
    "Plasma_cell"         : ["Jchain", "Mzb1", "Igkc", "Sdc1", "Ighg1","Igha"],            
    "Macrophage"          : ["Adgre1", "Cd68", "Csf1r", "C1qa", "C1qb", "Cx3cr1"],   
    "Dendritic_cell"      : ["Flt3", "Itgax", "Xcr1", "Siglech", "H2-Ab1", "Cd209a","Clec9a","Clec10a"],                                                                 
    "Neutrophil"          : ["S100a8", "S100a9", "Ly6g", "Csf3r", "Cxcr2"],          
    "NK"                  : ["Nkg7", "Klrb1c", "Ncr1", "Gzma", "Gzmb", "Klrd1"],             
    "ILC"                 : ["Il7r", "Kit", "Thy1", "Gata3", "Rorc", "Il22"],               
    "Adamdec1_fibroblast" : ["Adamdec1", "Cd34", "Pi16"],  
    "Pdgfra_fibroblast"   : ["Pdgfra", "Foxl1", "Bmp4", "Grem1", "Wnt5a"],    
    "Smooth_muscle"       : ["Acta2", "Myh11", "Tagln", "Des", "Cnn1"],           
} 

# save marker plots

os.makedirs(os.path.join(figdir,"marker_plots"),exist_ok=True)
sc.settings.figdir=os.path.join(figdir,"marker_plots")
for markers in marker_dict:

    sc.pl.umap(
        adata,
        color = [marker for marker in marker_dict[markers] if marker in adata.var_names],
        cmap = "Oranges_r",
        save=f"_{markers}.png"
    )

sc.settings.figdir=figdir