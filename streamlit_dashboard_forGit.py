import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np
import os
import json
import scanpy as sc
from scipy.stats import spearmanr, ttest_ind, fisher_exact, chi2_contingency

# --------- CONFIGURATION ---------
RNA_DATA_DIR = "SEA_AD_analysis/results/aggregated_multiomic/subclass_split_rna/filtered"
ATAC_DATA_DIR = "SEA_AD_analysis/results/aggregated_multiomic/subclass_split_atac/filtered"

CELLTYPE_TO_FILE = {
    "L2/3 IT": "scored_subclass_L2_3_IT_rna.h5ad",
    "Astrocyte": "scored_subclass_Astrocyte_rna.h5ad",
    "Oligodendrocyte": "scored_subclass_Oligodendrocyte_rna.h5ad",
    "L6 CT": "scored_subclass_L6_CT_rna.h5ad",
    "Lamp5": "scored_subclass_Lamp5_rna.h5ad",
    "Vip": "scored_subclass_Vip_rna.h5ad",
    "L4 IT": "scored_subclass_L4_IT_rna.h5ad",
    "Pvalb": "scored_subclass_Pvalb_rna.h5ad",
    "L6 IT": "scored_subclass_L6_IT_rna.h5ad",
    "Endothelial": "scored_subclass_Endothelial_rna.h5ad",
    "L5 IT": "scored_subclass_L5_IT_rna.h5ad",
    "Sncg": "scored_subclass_Sncg_rna.h5ad",
    "Microglia-PVM": "scored_subclass_Microglia-PVM_rna.h5ad",
    "Lamp5 Lhx6": "scored_subclass_Lamp5_Lhx6_rna.h5ad",
    "L6 IT Car3": "scored_subclass_L6_IT_Car3_rna.h5ad",
    "Sst": "scored_subclass_Sst_rna.h5ad",
    "Pax6": "scored_subclass_Pax6_rna.h5ad",
    "L5/6 NP": "scored_subclass_L5_6_NP_rna.h5ad",
    "L6b": "scored_subclass_L6b_rna.h5ad",
    "Chandelier": "scored_subclass_Chandelier_rna.h5ad",
    "OPC": "scored_subclass_OPC_rna.h5ad",
    "VLMC": "scored_subclass_VLMC_rna.h5ad",
    "L5 ET": "scored_subclass_L5_ET_rna.h5ad",
    "Sst Chodl": "scored_subclass_Sst_Chodl_rna.h5ad",
}

ALLOWED_FEATURES = [
    "Donor ID", "sex", "Gender", "age_death", "Age at Death", "Race", "specify other race",
    "Hispanic/Latino", "APOE_gt", "Primary Study Name", "Secondary Study Name",
    "Highest level of education", "Years of education", "pmi",
    "Fresh Brain Weight", "Brain pH", "RIN",
    "AD_pathology", "Overall AD neuropathological Change", "Thal",
    "Braak", "CERAD score", "plaque_d", "Overall CAA Score", "Highest Lewy Body Disease",
    "Total Microinfarcts (not observed grossly)", "Total microinfarcts in screening sections",
    "Atherosclerosis", "Arteriolosclerosis", "LATE",
    "diagnosis", "Cognitive Status", "Last CASI Score", "Interval from last CASI in months",
    "Last MMSE Score", "Interval from last MMSE in months", "Last MOCA Score", "Interval from last MOCA in months",
    "Used in analysis", "NeuN positive fraction on FANS", "leiden", "Class confidence", "Class",
    "Subclass confidence", "Subclass", "Supertype confidence", "Supertype (non-expanded)", "Supertype",
    "Continuous Pseudo-progression Score", "Severely Affected Donor",
]

def atac_file_from_celltype(celltype):
    name = celltype.replace("/", "_").replace(" ", "_")
    return f"scored_subclass_{name}_atac.h5ad"

st.set_page_config(
    page_title="Multiomic AD Dashboard",
    layout="wide",
    initial_sidebar_state="expanded",
)

def build_mask_from_criteria(obs, criteria_list):
    mask = np.full(len(obs), True)
    for crit in criteria_list:
        if not crit or not crit.get("feature"):
            continue
        feature = crit["feature"]
        if crit["type"] == "continuous":
            op = crit["operator"]
            cutoff = crit["cutoff"]
            if op == ">":
                mask = mask & (obs[feature] > cutoff)
            elif op == "<":
                mask = mask & (obs[feature] < cutoff)
            else:
                mask = mask & (obs[feature] == cutoff)
        elif crit["type"] == "categorical":
            cats = crit["categories"]
            if cats:
                mask = mask & (obs[feature].astype(str).isin(cats))
    return mask

@st.cache_data(show_spinner="Loading AnnData...")
def load_anndata(path, ucell_suffix=None):
    adata = sc.read_h5ad(path)
    if ucell_suffix is not None:
        ucell_cols = [col for col in adata.obs.columns if "_UCell" in col and not col.endswith(ucell_suffix)]
        adata.obs.rename(columns={col: col+ucell_suffix for col in ucell_cols}, inplace=True)
    obs_df = adata.obs
    if "Continuous Pseudo-progression Score" in obs_df.columns:
        cps = obs_df["Continuous Pseudo-progression Score"].astype(float)
        adata.obs["twoGroupCPS"] = np.select(
            [cps > 0.50, cps <= 0.50],
            ["CPS_high", "CPS_low"],
            default="NA"
        )
        adata.obs["threeGroupCPS"] = np.select(
            [cps > 0.66, (cps <= 0.66) & (cps > 0.33), cps <= 0.33],
            ["CPS_high", "CPS_moderate", "CPS_low"],
            default="NA"
        )
        adata.obs["fourGroupCPS"] = np.select(
            [cps > 0.66,
             (cps <= 0.66) & (cps > 0.33),
             (cps <= 0.33) & (cps > 0),
             cps == 0],
            ["CPS_high", "CPS_moderate", "CPS_low", "CPS_none"],
            default="NA"
        )
        age_col = "Age at Death" if "Age at Death" in obs_df.columns else "age_death"
        if age_col in obs_df.columns:
            age = pd.to_numeric(obs_df[age_col], errors="coerce")
            cps_by_age = cps.replace(0, np.nan) / age
            cps_by_age_valid = cps_by_age.replace([np.inf, -np.inf], np.nan).dropna()
            min_abc, max_abc = cps_by_age_valid.min(), cps_by_age_valid.max()
            scaled_abc = (cps_by_age - min_abc) / (max_abc - min_abc)
            adata.obs["CPS by Age"] = scaled_abc
            abc = scaled_abc
            adata.obs["twoGroupCPSByAge"] = np.select(
                [abc > 0.50, abc <= 0.50],
                ["CPSByAge_high", "CPSByAge_low"],
                default="NA"
            )
            adata.obs["threeGroupCPSByAge"] = np.select(
                [abc > 0.66, (abc <= 0.66) & (abc > 0.33), abc <= 0.33],
                ["CPSByAge_high", "CPSByAge_moderate", "CPSByAge_low"],
                default="NA"
            )
            adata.obs["fourGroupCPSByAge"] = np.select(
                [abc > 0.66,
                 (abc <= 0.66) & (abc > 0.33),
                 (abc <= 0.33) & (abc > 0),
                 abc == 0],
                ["CPSByAge_high", "CPSByAge_moderate", "CPSByAge_low", "CPSByAge_none"],
                default="NA"
            )
    return adata

def run_umap_pipeline(adata, n_top_genes=2000):
    adata = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3")
    adata = adata[:, adata.var.highly_variable]
    sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(adata, use_rep='X_pca', n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    return adata

def parse_obs_features(obs, max_categories=10):
    continuous = []
    categorical = []
    for col in obs.columns:
        series = obs[col]
        n_unique = series.nunique(dropna=True)
        if pd.api.types.is_categorical_dtype(series) or pd.api.types.is_bool_dtype(series):
            categorical.append(col)
        elif pd.api.types.is_numeric_dtype(series) and n_unique > max_categories:
            continuous.append(col)
        elif pd.api.types.is_object_dtype(series) or n_unique <= max_categories:
            categorical.append(col)
        else:
            categorical.append(col)
    return continuous, categorical

def is_continuous(col, obs):
    return pd.api.types.is_numeric_dtype(obs[col]) and obs[col].nunique() > 10

def add_significance_annotation(fig, text, x=0.05, y=0.95):
    fig.add_annotation(
        xref="paper", yref="paper", x=x, y=y, showarrow=False,
        text=text, font=dict(size=14, color="red"), bgcolor="white", opacity=0.8
    )
    return fig

# --- STATE: TRACK DATA LOAD BUTTON, CELL TYPE, ETC. ---
if "cell_type" not in st.session_state:
    st.session_state.cell_type = None
if "data_loaded" not in st.session_state:
    st.session_state.data_loaded = False
if "survival_model_table" not in st.session_state:
    st.session_state.survival_model_table = pd.DataFrame()
cell_type_sidebar = st.selectbox("Choose a cell type", list(CELLTYPE_TO_FILE.keys()))
if cell_type_sidebar != st.session_state.cell_type:
    st.session_state.data_loaded = False
    st.session_state.cell_type = cell_type_sidebar
cell_type = st.session_state.cell_type

h5ad_path = os.path.join(RNA_DATA_DIR, CELLTYPE_TO_FILE[cell_type])
atac_file = os.path.join(ATAC_DATA_DIR, atac_file_from_celltype(cell_type))

# --- SIDEBAR: LOGOUT BUTTON ---
#if st.sidebar.button('Logout'):
#    authenticator.logout('main')
#    st.rerun()
# --- SIDEBAR: LOAD BUTTON ---
with st.sidebar:
    st.title("🧰 Controls")
    st.subheader("🧠 Cell Type")
    st.write(f"RNA file: `{h5ad_path}`")
    st.write(f"ATAC file: `{atac_file}`")
    load_data_clicked = st.button("📥 Load Data")

# --- ONLY LOAD DATA ON BUTTON PRESS ---
adata, atac_data = None, None
if load_data_clicked:
    st.session_state.data_loaded = False
    st.session_state.adata = None
    st.session_state.adata_path = None
    st.session_state.atac_data = None
    st.session_state.atac_data_path = None
    st.session_state.filters = []
    st.session_state.groupA = []
    st.session_state.plots = []
    st.session_state.table_dfs = []
    st.session_state.last_umap_run = {"key": None, "result": None}
    st.session_state.survival_model_table = pd.DataFrame()
    with st.spinner(f"Loading RNA/ATAC for {cell_type} ..."):
        try:
            adata = load_anndata(h5ad_path, ucell_suffix="_rna")
            atac_data = load_anndata(atac_file, ucell_suffix="_atac")
            common_idx = adata.obs.index.intersection(atac_data.obs.index)
            adata._inplace_subset_obs(common_idx)
            atac_data._inplace_subset_obs(common_idx)
            st.session_state.adata = adata
            st.session_state.atac_data = atac_data
            st.session_state.adata_path = h5ad_path
            st.session_state.atac_data_path = atac_file
            st.session_state.data_loaded = True
            st.success("Loaded and aligned RNA & ATAC data!")
        except Exception as e:
            st.session_state.adata = None
            st.session_state.atac_data = None
            st.session_state.adata_path = None
            st.session_state.atac_data_path = None
            st.session_state.data_loaded = False
            st.error(f"Failed to load AnnData: {e}")

elif (
    st.session_state.data_loaded
    and st.session_state.adata is not None
    and st.session_state.atac_data is not None
    and st.session_state.adata_path == h5ad_path
    and st.session_state.atac_data_path == atac_file
):
    adata = st.session_state.adata
    atac_data = st.session_state.atac_data
else:
    adata = None
    atac_data = None

# --- REST OF APP: ONLY SHOW IF DATA IS LOADED ---
if adata is not None and atac_data is not None:
    obs = adata.obs
    all_obs_continuous, all_obs_categorical = parse_obs_features(obs)
    continuous_features = [f for f in all_obs_continuous if f in ALLOWED_FEATURES or "_UCell" in f or "by CPS" in f or "GroupCPS" in f]
    categorical_features = [f for f in all_obs_categorical if f in ALLOWED_FEATURES or "_UCell" in f or "by CPS" in f or "GroupCPS" in f]
    atac_ucell_cols = [col for col in atac_data.obs.columns if col.endswith("_UCell_atac")]
    for col in ["twoGroupCPS", "threeGroupCPS", "fourGroupCPS", "twoGroupCPSByAge", "threeGroupCPSByAge", "fourGroupCPSByAge"]:
        if col in obs.columns and col not in categorical_features:
            categorical_features.append(col)
    for col in ["CPS by Age"]:
        if col in obs.columns and col not in continuous_features:
            continuous_features.append(col)
    feature_list = sorted(set(continuous_features + categorical_features)) + sorted(atac_ucell_cols)
    if "User_selected_group" not in feature_list:
        feature_list.append("User_selected_group")
        if "User_selected_group" not in categorical_features:
            categorical_features.append("User_selected_group")

    # --------- SIDEBAR: FILTER & GROUPS & MODULES ---------
    with st.sidebar:
        st.subheader("🔎 Filtering")
        if "filters" not in st.session_state or st.session_state.filters is None:
            st.session_state.filters = []
        filter_col, filter_op_col = st.columns([3, 1])
        with filter_col:
            selected_feature = st.selectbox("Feature to filter on", feature_list, key="filter_feat")
        filter_kwargs = {}
        if selected_feature in continuous_features or selected_feature in atac_ucell_cols:
            with filter_op_col:
                operator = st.selectbox("Operator", [">", "<", "="], key="filter_op")
            min_val = float(np.nanmin(obs[selected_feature]))
            max_val = float(np.nanmax(obs[selected_feature]))
            cutoff = st.slider("Cutoff", min_val, max_val, float((min_val + max_val)/2), key="filter_cut")
            filter_kwargs = {"feature": selected_feature, "type": "continuous", "operator": operator, "cutoff": cutoff}
        elif selected_feature in categorical_features:
            with filter_op_col:
                st.write("")
            options = sorted([str(x) for x in obs[selected_feature].dropna().unique()])
            selected_categories = st.multiselect("Select categories", options, key="filter_cat")
            filter_kwargs = {"feature": selected_feature, "type": "categorical", "categories": selected_categories}
        col_apply, col_reset = st.columns([1,1])
        with col_apply:
            if st.button("➕ Add Filter", key="add_filter"):
                st.session_state.filters.append(filter_kwargs)
        with col_reset:
            if st.button("♻️ Reset Filters", key="reset_filters"):
                st.session_state.filters = []
        if st.session_state.filters:
            st.markdown("**Active Filters:**")
            for i, f in enumerate(st.session_state.filters):
                if f["type"] == "continuous":
                    st.markdown(f"• `{f['feature']} {f['operator']} {f['cutoff']}`")
                else:
                    cats = ', '.join(map(str, f['categories']))
                    st.markdown(f"• `{f['feature']} in [{cats}]`")
        # GROUP STACK
        st.subheader("🧬 Define Group (A)")
        if "groupA" not in st.session_state or st.session_state.groupA is None:
            st.session_state.groupA = []
        groupA_col, groupA_op_col = st.columns([3, 1])
        with groupA_col:
            group_feat_a = st.selectbox("Feature", feature_list, key="grpA_feat")
        groupA_kwargs = {}
        if group_feat_a in continuous_features or group_feat_a in atac_ucell_cols:
            with groupA_op_col:
                op_a = st.selectbox("Operator", [">", "<", "="], key="grpA_op")
            min_val_a = float(np.nanmin(obs[group_feat_a]))
            max_val_a = float(np.nanmax(obs[group_feat_a]))
            cut_a = st.slider("Cutoff", min_val_a, max_val_a, float((min_val_a + max_val_a)/2), key="grpA_cut")
            groupA_kwargs = {"feature": group_feat_a, "type": "continuous", "operator": op_a, "cutoff": cut_a}
        elif group_feat_a in categorical_features:
            with groupA_op_col:
                st.write("")
            options_a = sorted([str(x) for x in obs[group_feat_a].dropna().unique()])
            cats_a = st.multiselect("Categories", options_a, key="grpA_cat")
            groupA_kwargs = {"feature": group_feat_a, "type": "categorical", "categories": cats_a}
        col_apply_A, col_reset_A = st.columns([1,1])
        with col_apply_A:
            if st.button("➕ Add Group Criterion", key="add_groupA"):
                st.session_state.groupA.append(groupA_kwargs)
        with col_reset_A:
            if st.button("♻️ Reset Group", key="reset_groupA"):
                st.session_state.groupA = []
        if st.session_state.groupA:
            st.markdown("**Active Group Criteria:**")
            for i, g in enumerate(st.session_state.groupA):
                if g["type"] == "continuous":
                    st.markdown(f"• `{g['feature']} {g['operator']} {g['cutoff']}`")
                else:
                    cats = ', '.join(map(str, g['categories']))
                    st.markdown(f"• `{g['feature']} in [{cats}]`")
        st.markdown("_Comparison group: All other filtered samples not in Group A._")

        # ANALYSIS MODULE selection (always last in sidebar)
        st.subheader("⚙️ Analysis Module")
        analysis_module = st.radio(
            "Choose an analysis", [
                "Differential Expression",
                "Differential Accessibility",
                "Enrichment",
                "Feature Comparison",
                "Survival / Progression Modeling",
                "Dimensionality Reduction",
            ],
        )

        # Module-specific controls
        if analysis_module == "Feature Comparison":
            st.subheader("Feature Comparison Options")
            var1 = st.selectbox("Variable 1", feature_list, key="feat_comp_var1")
            var2 = st.selectbox("Variable 2", feature_list, key="feat_comp_var2")
            var3 = st.selectbox("Color by (optional)", ["None"] + feature_list, key="feat_comp_var3")
            pseudobulk = st.checkbox("Pseudobulk (aggregate by donor)", key="feat_comp_pseudobulk")
            run_feat_comparison = st.button("Run Feature Comparison")
        else:
            run_feat_comparison = False

        if analysis_module == "Dimensionality Reduction":
            st.subheader("UMAP Options")
            data_type = st.selectbox("Select data type", ["RNA", "ATAC"], key="dimred_data_type")
            umap_color_by = st.multiselect("Color UMAP by", feature_list, key="umap_color_selector")
            n_top_genes = st.slider(
                "Number of Top Variable Genes", min_value=500, max_value=5000, value=2000, step=100, key="n_top_genes_slider"
            )
        if analysis_module == "Survival / Progression Modeling":
            st.subheader("Pseudo-bulk Linear Modeling")
            outcome_choices = ["Continuous Pseudo-progression Score", "CPS by Age"]
            outcome_var = st.selectbox(
                "Select progression outcome variable:",
                outcome_choices,
                index=0
            )
            st.session_state.surv_outcome = outcome_var
            run_surv_analysis = st.button("Run Survival/Progression Modeling")
        else:
            run_surv_analysis = False

    # --------- FILTERING, GROUPING, INFO BAR ---------
    obs = adata.obs
    atac_obs = atac_data.obs
    filters = st.session_state.filters if st.session_state.filters is not None else []
    groupA = st.session_state.groupA if st.session_state.groupA is not None else []
    filter_mask = build_mask_from_criteria(obs, filters)
    obs_filtered = obs[filter_mask].copy()
    atac_obs_filtered = atac_obs[filter_mask].copy()
    donor_col = "Donor ID" if "Donor ID" in obs_filtered.columns else obs_filtered.columns[0]
    n_cells = len(obs_filtered)
    n_donors = obs_filtered[donor_col].nunique() if donor_col in obs_filtered else "?"
    groupA_mask = build_mask_from_criteria(obs_filtered, groupA)
    groupA_indices = obs_filtered.index[groupA_mask]
    groupB_indices = obs_filtered.index[~groupA_mask]
    n_groupA = len(groupA_indices)
    n_groupB = len(groupB_indices)
    obs_filtered["User_selected_group"] = np.where(groupA_mask, 1, 0)
    adata.obs.loc[obs_filtered.index, "User_selected_group"] = obs_filtered["User_selected_group"]
    atac_obs_filtered["User_selected_group"] = np.where(groupA_mask, 1, 0)
    atac_data.obs.loc[atac_obs_filtered.index, "User_selected_group"] = atac_obs_filtered["User_selected_group"]

    # --------- INFO/STATUS BAR ---------
    status_container = st.container()
    with status_container:
        st.markdown(
            f"**Cell Type:** {cell_type} | **Donors:** {n_donors} | **Cells (filtered):** {n_cells} | "
            f"**Group A:** {n_groupA} | **Comparison:** {n_groupB}"
        )
        st.markdown("---")

    st.title("Multiomic Alzheimer’s Disease Dashboard (Scaffold)")

    # --- ANALYSIS EXECUTION ---
    if analysis_module == "Dimensionality Reduction" and 'run_analysis' in locals() and run_analysis:
        with st.spinner("Generating UMAP plot..."):
            data_type = st.session_state.get("dimred_data_type", "RNA")
            data_for_umap = adata if data_type == "RNA" else atac_data
            obs_for_umap = obs_filtered if data_type == "RNA" else atac_obs_filtered
            key = (st.session_state.adata_path if data_type == "RNA" else st.session_state.atac_data_path,
                   json.dumps(st.session_state.filters, sort_keys=True),
                   n_top_genes, data_type)
            if st.session_state.last_umap_run.get('key') == key:
                st.toast("Using cached UMAP result.")
                adata_umap = st.session_state.last_umap_run.get('result')
            else:
                st.toast("Recalculating UMAP...")
                umap_adata = data_for_umap[obs_for_umap.index].copy()
                adata_umap = run_umap_pipeline(umap_adata, n_top_genes=n_top_genes)
                st.session_state.last_umap_run = {'key': key, 'result': adata_umap}
            df_coords = pd.DataFrame(adata_umap.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'], index=adata_umap.obs.index)
            for color_feature in umap_color_by:
                plot_df = df_coords.copy()
                plot_df[color_feature] = adata_umap.obs[color_feature]
                fig = px.scatter(
                    plot_df, x='UMAP_1', y='UMAP_2', color=color_feature,
                    hover_name=plot_df.index
                )
                st.session_state.plots.append({
                    "fig": fig,
                    "title": f"UMAP ({data_type}) by {color_feature} ({len(plot_df)} cells)"
                })
            df_table_run = df_coords.copy()
            for color_feature in umap_color_by:
                df_table_run[color_feature] = adata_umap.obs[color_feature]
            st.session_state.table_dfs.append(df_table_run)
            st.toast(f"Added {len(umap_color_by)} new plot(s).")

    # --------- FEATURE COMPARISON (WITH PSEUDOBULK OPTION) ---------
    if analysis_module == "Feature Comparison" and run_feat_comparison:
        st.info("Generating feature comparison plot...")
        needed_atac_cols = []
        if var1.endswith("_atac") and var1 in atac_data.obs.columns:
            needed_atac_cols.append(var1)
        if var2.endswith("_atac") and var2 in atac_data.obs.columns:
            needed_atac_cols.append(var2)
        if var3 and var3 != "None" and var3.endswith("_atac") and var3 in atac_data.obs.columns:
            needed_atac_cols.append(var3)
        df = obs_filtered.copy()
        if needed_atac_cols:
            df = df.join(atac_data.obs.loc[df.index, needed_atac_cols])
        col1 = var1
        col2 = var2
        col3 = var3 if var3 != "None" else None
        var1_cont = is_continuous(col1, df)
        var2_cont = is_continuous(col2, df)
        try:
            # Pseudobulk logic
            if 'pseudobulk' in locals() and pseudobulk and var1_cont and var2_cont:
                donor_col = "Donor ID"
                if donor_col not in df.columns:
                    st.error("No donor ID column for pseudobulk aggregation.")
                else:
                    # Calculate mean and standard error for each donor
                    agg_funcs = {col1: ['mean', 'sem'], col2: ['mean', 'sem']}
                    pseudobulk_stats = df.groupby(donor_col)[[col1, col2]].agg(agg_funcs)
                    pseudobulk_df = pd.DataFrame({
                        col1: pseudobulk_stats[(col1, 'mean')],
                        col2: pseudobulk_stats[(col2, 'mean')],
                        f"{col1}_sem": pseudobulk_stats[(col1, 'sem')],
                        f"{col2}_sem": pseudobulk_stats[(col2, 'sem')],
                    }).dropna()
                    if col3 is not None and col3 != "None":
                        st.warning("Color by is ignored when pseudobulk is enabled.")
                    fig = px.scatter(
                        pseudobulk_df, x=col1, y=col2, hover_name=pseudobulk_df.index,
                        error_x=f"{col1}_sem", error_y=f"{col2}_sem"
                    )
                    mask = pseudobulk_df[[col1, col2]].dropna().index
                    r, p = spearmanr(pseudobulk_df.loc[mask, col1], pseudobulk_df.loc[mask, col2])
                    r2 = r ** 2
                    text = f"Spearman ρ = {r:.2f}<br>R² = {r2:.2f}<br>p = {p:.2e}"
                    fig = add_significance_annotation(fig, text)
                    st.session_state.plots.append({"fig": fig, "title": f"Scatter (Pseudobulk): {col1} vs {col2}"})
                    st.plotly_chart(fig, use_container_width=True)
            else:
                if var1_cont and var2_cont:
                    fig = px.scatter(
                        df, x=col1, y=col2, color=col3 if col3 else None,
                        hover_name=df.index,
                        labels={col1: col1, col2: col2, (col3 if col3 else ""): col3}
                    )
                    mask = df[[col1, col2]].dropna().index
                    r, p = spearmanr(df.loc[mask, col1], df.loc[mask, col2])
                    r2 = r ** 2
                    text = f"Spearman ρ = {r:.2f}<br>R² = {r2:.2f}<br>p = {p:.2e}"
                    fig = add_significance_annotation(fig, text)
                    st.session_state.plots.append({"fig": fig, "title": f"Scatter: {col1} vs {col2}"})
                    st.plotly_chart(fig, use_container_width=True)
                elif (var1_cont and not var2_cont) or (not var1_cont and var2_cont):
                    y = col1 if var1_cont else col2
                    x = col2 if var1_cont else col1
                    fig = px.violin(
                        df, x=x, y=y, color=col3 if col3 else None, box=True, points="all",
                        hover_name=df.index,
                        labels={x: x, y: y, (col3 if col3 else ""): col3}
                    )
                    groups = df[x].dropna().unique()
                    pvals = {}
                    for i, group1 in enumerate(groups):
                        for group2 in groups[i+1:]:
                            vals1 = df[df[x] == group1][y].dropna()
                            vals2 = df[df[x] == group2][y].dropna()
                            if len(vals1) > 1 and len(vals2) > 1:
                                stat, pval = ttest_ind(vals1, vals2, equal_var=False)
                                pvals[(group1, group2)] = pval
                    if pvals:
                        from statsmodels.stats.multitest import multipletests
                        _, corrected, _, _ = multipletests(list(pvals.values()), method='bonferroni')
                        for i, ((g1, g2), raw_p) in enumerate(pvals.items()):
                            corr_p = corrected[i]
                            signif = ""
                            if corr_p < 0.001: signif = "***"
                            elif corr_p < 0.01: signif = "**"
                            elif corr_p < 0.05: signif = "*"
                            else: signif = "ns"
                            fig.add_annotation(
                                xref="paper", yref="paper", x=0.05, y=0.90 - i*0.05,
                                showarrow=False,
                                text=f"{g1} vs {g2}: p={corr_p:.2e} {signif}",
                                font=dict(size=12, color="blue"),
                                bgcolor="white", opacity=0.8
                            )
                    st.session_state.plots.append({"fig": fig, "title": f"Violin: {y} by {x}"})
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    ct_count = pd.crosstab(df[col1], df[col2])  # Raw counts for the statistical test
                    ct_prop = pd.crosstab(df[col1], df[col2], normalize='index')  # Proportions for plotting

                    fig = px.bar(
                        ct_prop,
                        x=ct_prop.index,
                        y=ct_prop.columns,
                        barmode="stack"
                    )
                    fig.update_yaxes(range=[0, 1], title="Fraction of Cells")

                    # Statistical test (use counts!)
                    if ct_count.shape == (2, 2):
                        _, pval = fisher_exact(ct_count)
                        testname = "Fisher's exact"
                    else:
                        _, pval, _, _ = chi2_contingency(ct_count)
                        testname = "Chi-square"
                    signif = ""
                    if pval < 0.001: signif = "***"
                    elif pval < 0.01: signif = "**"
                    elif pval < 0.05: signif = "*"
                    else: signif = "ns"
                    text = f"{testname} p={pval:.2e} {signif}"
                    fig = add_significance_annotation(fig, text)
                    st.session_state.plots.append({"fig": fig, "title": f"Stacked Bar: {col1} vs {col2}"})
                    st.plotly_chart(fig, use_container_width=True)
        except Exception as e:
            st.error(f"Error generating feature comparison plot: {e}")
    elif analysis_module not in ["Feature Comparison", "Dimensionality Reduction", "Survival / Progression Modeling"]:
        st.info(f"Analysis module '{analysis_module}' is not yet implemented.")

    # --- Survival / Progression Modeling (run only, results in Table tab) ---
    if analysis_module == "Survival / Progression Modeling" and run_surv_analysis:
        with st.spinner("Aggregating by Donor and running linear regression for each UCell score..."):
            outcome_var = st.session_state.surv_outcome if "surv_outcome" in st.session_state else "Continuous Pseudo-progression Score"
            rna_ucell_cols = [col for col in obs_filtered.columns if "_UCell" in col]
            atac_ucell_cols_in_filtered = [col for col in atac_obs_filtered.columns if "_UCell" in col]
            merged_df = obs_filtered[rna_ucell_cols + ["Donor ID", outcome_var]].copy()
            if atac_ucell_cols_in_filtered:
                cols_to_join = [col for col in atac_ucell_cols_in_filtered if col not in merged_df.columns]
                if cols_to_join:
                    merged_df = merged_df.join(atac_obs_filtered[cols_to_join], how="left")
            merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]
            donor_agg = merged_df.groupby("Donor ID").mean(numeric_only=True)
            donor_agg = donor_agg[~donor_agg[outcome_var].isna()]
            results = []
            import statsmodels.api as sm
            for ucell in rna_ucell_cols + atac_ucell_cols_in_filtered:
                if ucell not in donor_agg or donor_agg[ucell].isnull().all():
                    continue
                X = donor_agg[[ucell]].copy()
                X = sm.add_constant(X)
                y = donor_agg[outcome_var]
                try:
                    model = sm.OLS(y, X, missing="drop").fit()
                    coef = model.params[ucell]
                    pval = model.pvalues[ucell]
                    r2 = model.rsquared
                    results.append({
                        "UCell_Score": ucell,
                        "Coefficient": coef,
                        "p-value": pval,
                        "R²": r2,
                    })
                except Exception:
                    continue
            results_df = pd.DataFrame(results)
            if not results_df.empty:
                results_df = results_df.sort_values(by="Coefficient", ascending=False)
                st.session_state.survival_model_table = results_df
            else:
                st.session_state.survival_model_table = pd.DataFrame()

    # --- DISPLAY BLOCK ---
    col1, col2, _ = st.columns([1.2, 1, 8])
    with col1:
        if st.button("Reset Plots"):
            st.session_state.plots = []
            st.rerun()
    with col2:
        if st.button("Reset Table"):
            st.session_state.table_dfs = []
            st.session_state.survival_model_table = pd.DataFrame()
            st.rerun()

    tab_plot, tab_table = st.tabs(["Plot", "Table"])

    with tab_plot:
        if not st.session_state.plots:
            st.info("Configure options in the sidebar and click **Run Analysis** to generate plots.")
        else:
            cols = st.columns(2)
            for i, plot_info in enumerate(st.session_state.plots):
                with cols[i % 2]:
                    st.subheader(plot_info["title"])
                    st.plotly_chart(plot_info["fig"], use_container_width=True, key=f"plot_chart_{i}")

    with tab_table:
        show_any_table = False
        if analysis_module == "Survival / Progression Modeling":
            results_df = st.session_state.get("survival_model_table", pd.DataFrame())
            if results_df is not None and not results_df.empty:
                outcome_name = st.session_state.surv_outcome if "surv_outcome" in st.session_state else "Continuous Pseudo-progression Score"
                st.subheader(f"Predictive strength of UCell scores for {outcome_name}")
                st.dataframe(results_df, use_container_width=True)
                show_any_table = True
        if not show_any_table:
            if not st.session_state.table_dfs:
                st.info("No table data available. Run an analysis to generate data.")
            else:
                st.subheader("Data from Analysis Runs")
                for i, df in enumerate(st.session_state.table_dfs):
                    st.markdown(f"**Data from Run {i+1}** ({len(df)} cells)")
                    st.dataframe(df, use_container_width=True, key=f"table_df_{i}")

else:
    st.info("Select a cell type and click **Load Data** to begin.")
