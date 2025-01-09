# %%
import pandas as pd

# %%

dset_map = {
    "Abaumannii": "a_baum",
    "AfischeriES114": "a_fischeri_ES114",
    "Akunkeei": "a_kunk",
    "Banthracis": "b_anth",
    "Bfragilis": "b_frag",
    "Bpseudomallei": "b_pseudo",
    "Bsubtilis168": "b_sub_168",
    "Bthetaiotaomicron": "b_theta",
    "Bxylanisolvens": "b_xyla",
    "Cbasilensis": "c_basil",
    "Cdifficile": "c_diff",
    "Cvibrioides": "c_vibrioides",
    "DvulgarisHilden": "d_vulg",
    "Ecoli": "e_coli",
    "Efaecalis": "e_fae",
    "Elimosum": "e_limo",
    "Kpneumoniae": "k_pneum",
    "MsmegmatisMC2155": "m_tuber",
    "MtuberculosisH37Rv": "m_smeg",
    "Ngonorrhoeae": "n_gonorr",
    "Paeruginosa": "p_aeru",
    "Pchlororaphisaureofaciens30-84": "p_cholor_aureo3084",
    "Pfluorescens": "p_fluor",
    "Psalmonis": "p_salmo",
    "Scoelicolor": "s_coelicolor",
    "Selongatus": "s_elon",
    "SentericaTyph": "s_enter_typh",
    "Sepidermidis": "s_epi",
    "Smeliloti": "s_meli",
    "Spcc6803": "s_spcc6803",
    "SsanguinisSK36": "s_sanguinis",
    "Xalbilineans": "x_albi",
    "Xoryzae": "x_ory",
}

# %%
df = pd.read_csv("inf_results.csv", dtype=str)
df["dataset"] = df["dataset"].map(dset_map).fillna(df["dataset"])

# %%
df.to_csv("inf_results.csv", index=False)
# %%
