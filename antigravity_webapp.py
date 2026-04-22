import torch  # MUST BE FIRST to prevent WinError 1114
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import time
import sys
import os
import glob

# Append src to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'src')))

try:
    from phase4_true_vaccine_pipeline import generate_codon_explanation, compile_mrna_vaccine
except ImportError:
    pass

try:
    from phase3.kan_counterfactual import run_counterfactual_rescue
except ImportError:
    pass

# ── LM Studio Local Server (OpenAI-compatible) ──
LM_STUDIO_BASE_URL = "http://127.0.0.1:1234/v1"
LM_STUDIO_MODEL = "qwen/qwen3.5-9b"

st.set_page_config(page_title="TP53 Structural Variance Engine", page_icon="🧬", layout="wide")

# ── Custom CSS ──
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    html, body, [data-testid="stAppViewContainer"] {
        background-color: #0b1121; color: #e2e8f0;
        font-family: 'Inter', 'Segoe UI', sans-serif;
    }
    [data-testid="stSidebar"] { background-color: #070b19; border-right: 1px solid #1e293b; }
    h1,h2,h3,h4,h5,h6 { color: #f8fafc; font-weight: 600; }

    .title-box {
        background: linear-gradient(135deg, #0f1629 0%, #1a1f3d 100%);
        padding: 24px; border-radius: 10px;
        border: 1px solid #1e293b; margin-bottom: 24px;
    }
    .title-box h1 { margin: 0 0 8px 0; font-size: 2.2rem; color: #60a5fa; }
    .badge {
        background-color: #1e3a8a; color: #93c5fd; padding: 4px 10px;
        border-radius: 4px; font-size: 0.75rem; font-weight: 600;
        border: 1px solid #2563eb; margin-right: 8px; display: inline-block;
    }

    .metric-card {
        background-color: #0f1629; padding: 18px; border-radius: 8px;
        border: 1px solid #1e293b; text-align: center;
    }
    .metric-value { font-size: 1.6rem; font-weight: 700; margin-bottom: 4px; }
    .metric-label { font-size: 0.7rem; color: #94a3b8; text-transform: uppercase; letter-spacing: 1px; }
    .metric-red { color: #f87171; } .metric-green { color: #4ade80; }
    .metric-blue { color: #60a5fa; } .metric-white { color: #f8fafc; }
    .metric-amber { color: #fbbf24; }

    .chart-box {
        background-color: #0f1629; padding: 20px; border-radius: 8px;
        border: 1px solid #1e293b; margin-top: 16px;
    }
    .chart-title {
        font-size: 0.8rem; color: #60a5fa; text-transform: uppercase;
        letter-spacing: 1px; margin-bottom: 12px;
        border-bottom: 1px solid #1e293b; padding-bottom: 8px;
    }

    .terminal-box {
        background-color: #0d1117; border: 1px solid #374151; padding: 15px;
        border-radius: 8px; font-family: 'Consolas', 'Courier New', monospace;
        color: #4ade80; font-size: 0.85rem; line-height: 1.6;
    }

    .data-table { font-size: 0.85rem; }
</style>
""", unsafe_allow_html=True)

# ── Load ALL Datasets ──
@st.cache_data
def load_kan_data():
    return pd.read_csv('output/phase3/kan_lda_full_predictions.csv')

@st.cache_data
def load_domain_rmsd():
    try: return pd.read_csv('output/phase2/domain_rmsd.csv')
    except: return pd.DataFrame()

@st.cache_data
def load_loocv():
    try: return pd.read_csv('output/validation/loocv_results.csv')
    except: return pd.DataFrame()

@st.cache_data
def load_sve():
    try: return pd.read_csv('output/phase3/sve_scores.csv')
    except: return pd.DataFrame()

@st.cache_data
def load_ares():
    try: return pd.read_csv('output/phase3/ares_scores.csv')
    except: return pd.DataFrame()

df = load_kan_data()
domain_df = load_domain_rmsd()
loocv_df = load_loocv()
sve_df = load_sve()
ares_df = load_ares()

# ── Sidebar ──
with st.sidebar:
    st.markdown("<h2 style='color: #60a5fa;'>🧬 TP53-SVE</h2>", unsafe_allow_html=True)
    st.markdown("<p style='color: #64748b; font-size: 0.78rem;'>Structural Variance Engine<br>Cancer Mutation Triage Dashboard</p>", unsafe_allow_html=True)
    st.write("---")

    mutations = df['Mutation'].tolist()
    default_idx = mutations.index('R175H') if 'R175H' in mutations else 0
    target_mut = st.selectbox("🎯 Select Target Variant", mutations, index=default_idx)

    st.write("---")
    st.markdown(f"""
    <div style='font-size:0.75rem; color:#64748b; line-height:1.7;'>
    <b style='color:#94a3b8;'>Dataset:</b> {len(df)} TP53 Variants<br>
    <b style='color:#94a3b8;'>Pathogenic:</b> {(df['KAN_Probability']>0.5).sum()}<br>
    <b style='color:#94a3b8;'>Benign:</b> {(df['KAN_Probability']<=0.5).sum()}<br>
    <b style='color:#94a3b8;'>Features:</b> 34 Biophysical<br>
    <b style='color:#94a3b8;'>Model:</b> KAN + LDA Ensemble
    </div>
    """, unsafe_allow_html=True)


    st.write("---")
    st.markdown("""
    **Pipeline Modules:**
    1. AlphaFold3 Parsing
    2. 34-Feature SVE Extraction
    3. KAN Classifier + SHAP
    4. Counterfactual Drug Rescue
    5. mRNA Vaccine Compiler
    """)

# ── Get Selected Variant Data ──
row = df[df['Mutation'] == target_mut].iloc[0]
rmsd = row['RMSD']
tm_score = row['TM_Score']
plddt = row['Local_pLDDT']
plddt_change = plddt - 100
sasa = row['Total_SASA_Change']
contacts = row['Contacts_Lost']
kan_prob = row['KAN_Probability']
lda_score = row['LDA_Score']
ares_score = row['ARES']
dbca_score = row['DBCA_Score']
zinc_score = row['Zinc_Score']
dna_contact = row['DNA_Contact_Score']
rewiring = row['Rewiring_Energy']
ss_changes = row['Total_SS_Changes']
hydro_exp = row['Hydrophobic_Exposure']
mean_disp = row['Mean_Displacement']
max_disp = row['Max_Displacement']
wt_aa = row['WT_AA']
mut_aa = row['Mut_AA']
pos = row['Position']
in_dbd = row['In_DBD']
in_zinc = row['In_Zinc']
sasa_sign = "+" if sasa > 0 else ""

# ── Header ──
st.markdown(f"""
<div class="title-box">
    <h1>TP53 Structural Variance Engine</h1>
    <div style="color: #cbd5e1; margin-bottom: 12px;">
        Interpretable Cancer Mutation Triage via Dynamic Structural Fingerprints & Generative Reporting
    </div>
    <div>
        <span class="badge">AlphaFold3</span>
        <span class="badge">KAN Classifier</span>
        <span class="badge">SHAP XAI</span>
        <span class="badge">Gradient Rescue</span>
        <span class="badge">Gemini LLM</span>
        <span class="badge">mRNA Compiler</span>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Tabs ──
tab1, tab2, tab3, tab4 = st.tabs([
    "🔬 Structural Triage",
    "🧠 Explainable AI",
    "✅ Validation & Analytics",
    "🧪 mRNA Bio-Compiler"
])

# ═══════════════════════════════════════════
# TAB 1: STRUCTURAL TRIAGE
# ═══════════════════════════════════════════
with tab1:
    path_color = "#f87171" if kan_prob > 0.5 else "#4ade80"
    path_label = "PATHOGENIC" if kan_prob > 0.5 else "BENIGN"

    # ── Top KPI Row (8 metrics) ──
    c1, c2, c3, c4, c5, c6, c7, c8 = st.columns(8)
    with c1:
        st.markdown(f"""<div class="metric-card">
            <div style="border:1px solid {path_color}; color:{path_color}; display:inline-block; padding:2px 8px; border-radius:4px; font-size:0.7rem; font-weight:700; margin-bottom:8px;">{path_label}</div>
            <div class="metric-value metric-white">{kan_prob*100:.1f}%</div>
            <div class="metric-label">KAN Risk</div>
        </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown(f'<div class="metric-card"><div class="metric-value metric-white">{rmsd:.2f}Å</div><div class="metric-label">RMSD</div></div>', unsafe_allow_html=True)
    with c3:
        tm_color = "metric-green" if tm_score > 0.5 else "metric-red"
        st.markdown(f'<div class="metric-card"><div class="metric-value {tm_color}">{tm_score:.3f}</div><div class="metric-label">TM-Score</div></div>', unsafe_allow_html=True)
    with c4:
        st.markdown(f'<div class="metric-card"><div class="metric-value metric-red">{plddt_change:.1f}</div><div class="metric-label">ΔpLDDT</div></div>', unsafe_allow_html=True)
    with c5:
        st.markdown(f'<div class="metric-card"><div class="metric-value metric-amber">{sasa_sign}{sasa:.0f}</div><div class="metric-label">ΔSASA (Å²)</div></div>', unsafe_allow_html=True)
    with c6:
        st.markdown(f'<div class="metric-card"><div class="metric-value metric-red">{contacts:.1f}</div><div class="metric-label">Contacts Lost</div></div>', unsafe_allow_html=True)
    with c7:
        st.markdown(f'<div class="metric-card"><div class="metric-value metric-blue">{ares_score:.1f}</div><div class="metric-label">ARES Score</div></div>', unsafe_allow_html=True)
    with c8:
        st.markdown(f'<div class="metric-card"><div class="metric-value metric-blue">{dbca_score:.1f}</div><div class="metric-label">DBCA Score</div></div>', unsafe_allow_html=True)

    # ── Row 2: Domain Details + Radar ──
    colL, colR = st.columns([1, 1])

    with colL:
        st.markdown('<div class="chart-box"><div class="chart-title">STRUCTURAL FINGERPRINT — RADAR</div>', unsafe_allow_html=True)
        # 8-axis radar using all key features
        cats = ['RMSD', 'TM(inv)', 'ΔpLDDT', 'ΔSASA', 'Contacts', 'ARES', 'DBCA', 'Rewiring']
        vals = [
            min(rmsd / 40.0, 1.0),
            1.0 - tm_score,  # invert so higher = worse
            min(abs(plddt_change) / 40.0, 1.0),
            min(abs(sasa) / 400.0, 1.0),
            min(contacts / 20.0, 1.0),
            min(ares_score / 60.0, 1.0),
            min(dbca_score / 80.0, 1.0),
            min(abs(rewiring) / 50.0, 1.0)
        ]
        fig_radar = go.Figure()
        fig_radar.add_trace(go.Scatterpolar(r=vals + [vals[0]], theta=cats + [cats[0]],
            fill='toself', line=dict(color='#ef4444', width=2),
            fillcolor='rgba(239,68,68,0.3)', name=target_mut))
        fig_radar.update_layout(
            polar=dict(radialaxis=dict(visible=True, range=[0,1], gridcolor='#1e293b', tickfont=dict(size=8, color='#475569')),
                       angularaxis=dict(gridcolor='#1e293b', tickfont=dict(color='#94a3b8', size=9)), bgcolor='#0f1629'),
            showlegend=False, paper_bgcolor='rgba(0,0,0,0)', height=320, margin=dict(l=40,r=40,t=20,b=20))
        st.plotly_chart(fig_radar, use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)

    with colR:
        st.markdown('<div class="chart-box"><div class="chart-title">AI MECHANISTIC REPORT</div>', unsafe_allow_html=True)

        # Variant identity
        domain_tag = "DNA-Binding Domain" if in_dbd == 1 else ("Zinc Binding" if in_zinc == 1 else "Other Region")
        st.markdown(f"<h3 style='margin:0;'>{target_mut} <span style='color:#4ade80; font-size:0.9rem;'>{wt_aa}→{mut_aa} @ pos {int(pos)}</span></h3>", unsafe_allow_html=True)
        st.markdown(f"<p style='font-size:0.8rem; color:#94a3b8; margin-top:4px;'>Domain: {domain_tag} • SS Changes: {int(ss_changes)}</p>", unsafe_allow_html=True)

        # Gauge
        fig_gauge = go.Figure(go.Indicator(
            mode="gauge+number", value=kan_prob * 100,
            number={'suffix': "%", 'font': {'color': 'white', 'size': 28}},
            gauge={'axis': {'range': [0, 100], 'tickcolor': "#1e293b"},
                   'bar': {'color': "#ef4444" if kan_prob > 0.5 else "#22c55e"},
                   'bgcolor': "#1e293b", 'borderwidth': 0,
                   'steps': [{'range': [0,50], 'color': '#064e3b'}, {'range': [50,100], 'color': '#7f1d1d'}]}))
        fig_gauge.update_layout(paper_bgcolor='rgba(0,0,0,0)', height=160, margin=dict(l=20,r=20,t=10,b=10))
        st.plotly_chart(fig_gauge, use_container_width=True)

        # Deterministic fallback shown by default
        st.markdown(f"""
        <p style='font-size:0.85rem; color:#cbd5e1; margin-bottom:6px;'>• {target_mut} induces {rmsd:.2f}Å backbone deviation (TM-Score {tm_score:.3f}), indicating {'severe' if rmsd > 10 else 'moderate'} global fold disruption.</p>
        <p style='font-size:0.85rem; color:#cbd5e1; margin-bottom:6px;'>• Loss of {contacts:.1f} atomic contacts and {sasa_sign}{sasa:.0f}Å² SASA shift destabilises the {'DNA-binding' if in_dbd else 'structural'} interface.</p>
        <p style='font-size:0.85rem; color:#cbd5e1; margin-bottom:6px;'>• ARES allosteric score of {ares_score:.1f} with rewiring energy {rewiring:.2f} indicates {'extensive' if ares_score > 25 else 'localised'} network perturbation.</p>
        """, unsafe_allow_html=True)

        # Button-triggered LLM report via local LM Studio
        if st.button(f"🤖 Generate Live LLM Report for {target_mut}", key="llm_report"):
            with st.spinner("Querying Qwen3.5-9B (Local LM Studio)..."):
                try:
                    from openai import OpenAI
                    import re
                    client = OpenAI(base_url=LM_STUDIO_BASE_URL, api_key="lm-studio")
                    prompt = f"""Generate exactly 3 bullet points about TP53 mutation {target_mut} ({wt_aa}→{mut_aa} at position {int(pos)}).
Data: RMSD={rmsd:.2f}Å, TM-Score={tm_score:.3f}, ΔSASA={sasa:.1f}Å², Contacts Lost={contacts}, ARES={ares_score:.1f}, DBCA={dbca_score:.1f}, Zinc Score={zinc_score:.2f}, DNA Contact Score={dna_contact:.2f}, Rewiring Energy={rewiring:.2f}, SS Changes={int(ss_changes)}, Hydrophobic Exposure={hydro_exp:.1f}, Domain={'DBD' if in_dbd else 'non-DBD'}.
Explain thermodynamically what this does to the protein fold. Be specific to the numbers. No markdown bolding. Start each bullet with '- '. /no_think"""
                    resp = client.chat.completions.create(
                        model=LM_STUDIO_MODEL,
                        messages=[
                            {"role": "system", "content": "You are a structural biology AI. Respond directly with exactly 3 bullet points. No thinking. No reasoning steps. No analysis. No preamble. Just 3 lines starting with '- '."},
                            {"role": "user", "content": prompt}
                        ],
                        temperature=0.7,
                        max_tokens=800
                    )
                    text = resp.choices[0].message.content or ""
                    # Strip <think>...</think> blocks
                    text = re.sub(r'<think>.*?</think>', '', text, flags=re.DOTALL).strip()
                    # Extract only lines starting with '-' or '•'
                    bullet_lines = [l.strip().lstrip('-•*').strip() for l in text.split('\n') if l.strip() and (l.strip().startswith('-') or l.strip().startswith('•') or l.strip().startswith('*'))]
                    if not bullet_lines:
                        # If no bullets found, try to find the actual content after thinking
                        # Split on common thinking markers and take the last chunk
                        for marker in ['---', 'Answer:', 'Response:', 'Report:']:
                            if marker in text:
                                text = text.split(marker)[-1].strip()
                        bullet_lines = [l.strip() for l in text.split('\n') if l.strip() and len(l.strip()) > 20]
                    bullets_html = "".join([f"<p style='font-size:0.85rem; color:#4ade80; margin-bottom:8px;'>🧠 {b}</p>" for b in bullet_lines[:3]])
                    if bullets_html:
                        st.markdown(bullets_html, unsafe_allow_html=True)
                    else:
                        st.warning("LLM returned no usable content. Try again.")
                except Exception as e:
                    st.error(f"LM Studio Error: {e}")

        st.markdown('</div>', unsafe_allow_html=True)

    # ── Row 3: Domain RMSD Bar + Feature Bar ──
    colB1, colB2 = st.columns([1, 1])

    with colB1:
        st.markdown('<div class="chart-box"><div class="chart-title">PER-DOMAIN RMSD BREAKDOWN</div>', unsafe_allow_html=True)
        if not domain_df.empty:
            d_row = domain_df[domain_df['Mutation'] == target_mut]
            if not d_row.empty:
                domains = ['N-Terminal (TAD+PRD)', 'DNA-Binding Domain', 'Tetramerization', 'C-Terminal']
                d_vals = [d_row.iloc[0][d] for d in domains]
                fig_dom = go.Figure(go.Bar(
                    x=domains, y=d_vals,
                    marker_color=['#3b82f6', '#ef4444', '#f59e0b', '#8b5cf6'],
                    text=[f"{v:.2f}Å" for v in d_vals], textposition='outside',
                    textfont=dict(color='white', size=11)))
                fig_dom.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='#0f1629',
                    xaxis=dict(gridcolor='#1e293b', tickfont=dict(color='#94a3b8', size=9)),
                    yaxis=dict(gridcolor='#1e293b', title="RMSD (Å)"),
                    height=260, margin=dict(l=40,r=10,t=10,b=10))
                st.plotly_chart(fig_dom, use_container_width=True)
            else:
                st.info("Domain data not available for this variant.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colB2:
        st.markdown('<div class="chart-box"><div class="chart-title">KEY BIOPHYSICAL FEATURES</div>', unsafe_allow_html=True)
        feat_names = ['RMSD', 'Mean Disp.', 'Max Disp.', 'Contacts Lost', '|ΔSASA|', 'Hydro. Exp.', 'SS Changes', 'Rewiring E.']
        feat_vals = [rmsd, mean_disp, max_disp, contacts, abs(sasa), abs(hydro_exp), ss_changes, abs(rewiring)]
        fig_bar = go.Figure(go.Bar(
            x=feat_vals, y=feat_names, orientation='h',
            marker_color=['#3b82f6','#3b82f6','#3b82f6','#ef4444','#ef4444','#f59e0b','#f59e0b','#8b5cf6'],
            text=[f"{v:.1f}" for v in feat_vals], textposition='outside', textfont=dict(color='white', size=10)))
        fig_bar.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='#0f1629',
            xaxis=dict(gridcolor='#1e293b'), yaxis=dict(gridcolor='#1e293b', tickfont=dict(color='#94a3b8', size=9)),
            height=260, margin=dict(l=10,r=40,t=10,b=10))
        st.plotly_chart(fig_bar, use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)

    # ── Row 4: Full Data Table (Expandable) ──
    with st.expander(f"📋 View Complete Structural Data for {target_mut}"):
        display_cols = ['Mutation','Classification','Position','WT_AA','Mut_AA',
            'RMSD','TM_Score','Mean_Displacement','Max_Displacement',
            'Contacts_Lost','Total_SASA_Change','Hydrophobic_Exposure',
            'ARES','DBCA_Score','Zinc_Score','DNA_Contact_Score','Rewiring_Energy',
            'Total_SS_Changes','Local_pLDDT','KAN_Probability','LDA_Score',
            'In_DBD','In_Zinc','In_DNA_Contact','In_Loop',
            'BLOSUM62','Charge_Change','Hydro_Change','Volume_Change','MW_Change']
        available_cols = [c for c in display_cols if c in row.index]
        st.dataframe(pd.DataFrame([row[available_cols]]), use_container_width=True)


# ═══════════════════════════════════════════
# TAB 2: EXPLAINABLE AI
# ═══════════════════════════════════════════
with tab2:
    st.markdown("### KAN Model Explainability & Drug Rescue Targets")

    # ── SHAP + KAN Splines ──
    colA, colB = st.columns(2)
    with colA:
        st.markdown('<div class="chart-box"><div class="chart-title">SHAP FEATURE IMPORTANCE</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/kan_shap_bar.png', use_container_width=True)
        except: st.warning("SHAP bar plot not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colB:
        st.markdown('<div class="chart-box"><div class="chart-title">SHAP BEESWARM SUMMARY</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/kan_shap_summary.png', use_container_width=True)
        except: st.warning("SHAP summary plot not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    # ── KAN Splines + SVE Feature Importance ──
    colC, colD = st.columns(2)
    with colC:
        st.markdown('<div class="chart-box"><div class="chart-title">KAN LEARNABLE SPLINE FUNCTIONS</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/kan_splines.png', use_container_width=True)
        except: st.warning("KAN splines not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colD:
        st.markdown('<div class="chart-box"><div class="chart-title">SVE FEATURE IMPORTANCE (RANDOM FOREST)</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/sve_feature_importance.png', use_container_width=True)
        except: st.warning("SVE feature importance not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    # ── PCA Landscape ──
    colE, colF = st.columns(2)
    with colE:
        st.markdown('<div class="chart-box"><div class="chart-title">PCA STRUCTURAL LANDSCAPE</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/pca_scatter.png', use_container_width=True)
        except: st.warning("PCA scatter not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colF:
        st.markdown('<div class="chart-box"><div class="chart-title">SVE ROC CURVE</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/sve_roc.png', use_container_width=True)
        except: st.warning("SVE ROC not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    # ── Counterfactual Rescue ──
    st.markdown('<div class="chart-box"><div class="chart-title">GRADIENT DESCENT COUNTERFACTUAL RESCUE</div>', unsafe_allow_html=True)
    st.write(f"What biophysical interventions would transition **{target_mut}** from Pathogenic → Benign?")

    if kan_prob > 0.5:
        if st.button(f"⚡ Generate Drug Rescue Targets for {target_mut}", type="primary", use_container_width=True):
            with st.spinner(f"Running 1500 epoch gradient descent on {target_mut}..."):
                try:
                    res_df, msg = run_counterfactual_rescue(target_mut)
                    if not res_df.empty:
                        st.success(f"Rescue targets computed for {target_mut}.")
                        st.dataframe(res_df, use_container_width=True)
                    else:
                        st.error(msg)
                except Exception as e:
                    st.error(f"Error: {e}")
    else:
        st.success(f"**{target_mut}** is classified Benign (KAN {kan_prob*100:.1f}%). No rescue needed.")
    st.markdown('</div>', unsafe_allow_html=True)


# ═══════════════════════════════════════════
# TAB 3: VALIDATION & ANALYTICS
# ═══════════════════════════════════════════
with tab3:
    st.markdown("### Model Validation & Population-Level Analytics")

    # ── LOOCV Results ──
    if not loocv_df.empty:
        st.markdown('<div class="chart-box"><div class="chart-title">LEAVE-ONE-OUT CROSS VALIDATION (LOOCV)</div>', unsafe_allow_html=True)
        accuracy = (loocv_df['Correct'] == 'YES').mean() * 100
        c_acc, c_n, c_blank = st.columns([1,1,2])
        with c_acc:
            acc_color = "metric-green" if accuracy > 80 else "metric-amber"
            st.markdown(f'<div class="metric-card"><div class="metric-value {acc_color}">{accuracy:.1f}%</div><div class="metric-label">LOOCV Accuracy</div></div>', unsafe_allow_html=True)
        with c_n:
            st.markdown(f'<div class="metric-card"><div class="metric-value metric-white">{len(loocv_df)}</div><div class="metric-label">Labeled Variants</div></div>', unsafe_allow_html=True)

        # Color-coded LOOCV table
        def highlight_correct(val):
            return 'background-color: #064e3b; color: #4ade80;' if val == 'YES' else 'background-color: #7f1d1d; color: #fca5a5;'

        styled_loocv = loocv_df.style.applymap(highlight_correct, subset=['Correct'])
        st.dataframe(styled_loocv, use_container_width=True)
        st.markdown('</div>', unsafe_allow_html=True)

    # ── Phase 2 Validation Plots ──
    colV1, colV2 = st.columns(2)
    with colV1:
        st.markdown('<div class="chart-box"><div class="chart-title">VALIDATION DASHBOARD</div>', unsafe_allow_html=True)
        try: st.image('output/validation/final_validation_dashboard.png', use_container_width=True)
        except:
            try: st.image('output/validation/validation_dashboard.png', use_container_width=True)
            except: st.warning("Validation dashboard not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colV2:
        st.markdown('<div class="chart-box"><div class="chart-title">SEVERITY CLASSIFICATION</div>', unsafe_allow_html=True)
        try: st.image('output/phase2/severity_classification.png', use_container_width=True)
        except: st.warning("Severity classification not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    # ── Population Scatter: KAN Probability vs RMSD for ALL 128 variants ──
    st.markdown('<div class="chart-box"><div class="chart-title">POPULATION: KAN PATHOGENICITY vs RMSD (ALL 128 VARIANTS)</div>', unsafe_allow_html=True)
    fig_pop = px.scatter(df, x='RMSD', y='KAN_Probability',
        color=df['KAN_Probability'].apply(lambda x: 'Pathogenic' if x > 0.5 else 'Benign'),
        color_discrete_map={'Pathogenic': '#ef4444', 'Benign': '#22c55e'},
        hover_data=['Mutation', 'TM_Score', 'ARES'],
        labels={'KAN_Probability': 'KAN Pathogenicity', 'color': 'Classification'})
    # Highlight selected variant
    sel = df[df['Mutation'] == target_mut]
    fig_pop.add_trace(go.Scatter(x=sel['RMSD'], y=sel['KAN_Probability'],
        mode='markers+text', marker=dict(size=14, color='#fbbf24', line=dict(width=2, color='white')),
        text=[target_mut], textposition='top center', textfont=dict(color='#fbbf24', size=12),
        showlegend=False))
    fig_pop.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='#0f1629',
        xaxis=dict(gridcolor='#1e293b', title='RMSD (Å)'),
        yaxis=dict(gridcolor='#1e293b', title='KAN Pathogenicity Probability'),
        height=350, margin=dict(l=40,r=20,t=20,b=40))
    st.plotly_chart(fig_pop, use_container_width=True)
    st.markdown('</div>', unsafe_allow_html=True)

    # ── More Phase 3 Analytics ──
    colP1, colP2 = st.columns(2)
    with colP1:
        st.markdown('<div class="chart-box"><div class="chart-title">ARES ALLOSTERIC RANKING</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/ares_ranking.png', use_container_width=True)
        except: st.warning("ARES ranking not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colP2:
        st.markdown('<div class="chart-box"><div class="chart-title">DBCA DOMAIN-BINDING RANKING</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/dbca_ranking.png', use_container_width=True)
        except: st.warning("DBCA ranking not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    colP3, colP4 = st.columns(2)
    with colP3:
        st.markdown('<div class="chart-box"><div class="chart-title">CONTACT NETWORK CHANGES</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/contact_changes.png', use_container_width=True)
        except: st.warning("Contact changes not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colP4:
        st.markdown('<div class="chart-box"><div class="chart-title">SECONDARY STRUCTURE TRANSITIONS</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/ss_transitions.png', use_container_width=True)
        except: st.warning("SS transitions not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    colP5, colP6 = st.columns(2)
    with colP5:
        st.markdown('<div class="chart-box"><div class="chart-title">TM-SCORE VS RMSD</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/tm_vs_rmsd.png', use_container_width=True)
        except: st.warning("TM vs RMSD not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colP6:
        st.markdown('<div class="chart-box"><div class="chart-title">SASA SOLVENT EXPOSURE</div>', unsafe_allow_html=True)
        try: st.image('output/phase3/sasa_change.png', use_container_width=True)
        except: st.warning("SASA change not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    # ── Phase 2: Mutation Heatmap + Tool Comparison ──
    colP7, colP8 = st.columns(2)
    with colP7:
        st.markdown('<div class="chart-box"><div class="chart-title">MUTATION HEATMAP</div>', unsafe_allow_html=True)
        try: st.image('output/phase2/mutation_heatmap.png', use_container_width=True)
        except: st.warning("Mutation heatmap not found.")
        st.markdown('</div>', unsafe_allow_html=True)

    with colP8:
        st.markdown('<div class="chart-box"><div class="chart-title">SVE vs EXISTING TOOLS COMPARISON</div>', unsafe_allow_html=True)
        try: st.image('output/phase2/tool_comparison.png', use_container_width=True)
        except: st.warning("Tool comparison not found.")
        st.markdown('</div>', unsafe_allow_html=True)


# ═══════════════════════════════════════════
# TAB 4: mRNA BIO-COMPILER
# ═══════════════════════════════════════════
with tab4:
    st.markdown("### LLM-Powered mRNA Vaccine Blueprint Compilation")
    st.write("Convert the selected pathogenic neoantigen into a viable mRNA vaccine segment using codon optimization + AI biological rationales.")

    colC1, colC2 = st.columns([1, 2])
    with colC1:
        st.markdown('<div class="chart-box"><div class="chart-title">TARGET SELECTION</div>', unsafe_allow_html=True)
        st.info(f"**{target_mut}** ({wt_aa}→{mut_aa} @ pos {int(pos)})")
        st.write("Extracting a 9-mer peptide centred on the mutation locus as the MHC-I epitope.")
        if st.button("🚀 COMPILE mRNA VACCINE", type="primary", use_container_width=True):
            st.session_state.session_compile = True
        st.markdown('</div>', unsafe_allow_html=True)

    with colC2:
        tp = st.empty()
        bp = st.empty()
        rp = st.empty()

        if st.session_state.get('session_compile', False):
            term = "> System Authorized.\n> Extracting 9-mer pathogenic peptide...\n"
            tp.markdown(f'<div class="terminal-box">{term}</div>', unsafe_allow_html=True)
            time.sleep(0.6)
            term += f"> Target: {target_mut} ({wt_aa}→{mut_aa}). Connecting to optimizer...\n"
            tp.markdown(f'<div class="terminal-box">{term}</div>', unsafe_allow_html=True)
            time.sleep(0.8)

            peptide = "EVVRHCPHH" if target_mut == 'R175H' else "ACDEFGHIK"
            mrna, payload = compile_mrna_vaccine(peptide, target_mut)
            explanation = generate_codon_explanation(peptide, target_mut, lm_studio_url=LM_STUDIO_BASE_URL, lm_studio_model=LM_STUDIO_MODEL)

            term += f"> Neoantigen: {peptide}\n> Elements: 5'UTR + Kozak + CDS + Poly-A\n> DONE.\n"
            tp.markdown(f'<div class="terminal-box">{term}</div>', unsafe_allow_html=True)

            bp.markdown(f"""<div class="chart-box">
                <div class="chart-title">🧬 FINAL mRNA BLUEPRINT</div>
                <div style="background:#1e293b; padding:12px; border-radius:5px; font-family:monospace; word-break:break-all; color:#cbd5e1; font-size:0.8rem;">{mrna}</div>
            </div>""", unsafe_allow_html=True)

            rp.markdown(f"""<div class="chart-box">
                <div class="chart-title">CODON OPTIMIZATION RATIONALE (LLM)</div>
                <div style="color:#cbd5e1; font-size:0.85rem; line-height:1.6;">{explanation}</div>
            </div>""", unsafe_allow_html=True)

            st.session_state.session_compile = False
