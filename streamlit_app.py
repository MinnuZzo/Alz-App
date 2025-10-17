import streamlit as st
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
import networkx as nx
from pyvis.network import Network
import os

# --- Title and Intro ---
st.set_page_config(page_title="Alzheimer’s Pathway Explorer", layout="wide")
st.title("🧠 Alzheimer’s Disease Pathway Explorer")
st.markdown("Explore KEGG biomarker interactions and pathways related to Alzheimer's disease.")

# --- Biomarker List ---
biomarkers = ["APP", "BACE1", "PSEN1", "MAPT", "GSK3B",
              "APOE", "IL6", "TNF", "SOD1", "CYCS"]

selected_biomarkers = st.multiselect(
    "Select biomarkers to visualize:",
    options=biomarkers,
    default=biomarkers
)

# --- Load KEGG Pathway Data ---
xml_path = "hsa05010.xml"

if not os.path.exists(xml_path):
    with st.spinner("Fetching KEGG Alzheimer’s pathway..."):
        kgml_data = REST.kegg_get("hsa05010", "kgml").read()
        with open(xml_path, "w") as f:
            f.write(kgml_data)

try:
    pathway = KGML_parser.read(open(xml_path))
except Exception as e:
    st.error(f"❌ Failed to parse pathway file: {e}")
    st.stop()

G = nx.DiGraph()

# --- Build Graph ---
for entry in pathway.entries.values():
    if entry.type in ["gene", "enzyme", "compound"]:
        label = getattr(entry.graphics, "name", entry.name)
        G.add_node(entry.id, name=entry.name, label=label, type=entry.type)

for rel in pathway.relations:
    e1, e2 = rel.entry1, rel.entry2
    if e1 and e2:
        G.add_edge(e1.id, e2.id, type=rel.type)

# --- Fix node filtering (match biomarker names inside labels) ---
sub_nodes = [n for n, d in G.nodes(data=True)
             if any(bio.lower() in d.get("label", "").lower() for bio in selected_biomarkers)]

if not sub_nodes:
    st.warning("⚠️ No matching biomarkers found in the KEGG pathway. Showing full pathway instead.")
    subgraph = G
else:
    subgraph = G.subgraph(sub_nodes).copy()

# --- Visualization ---
net = Network(height="650px", width="100%", bgcolor="#111", font_color="white", notebook=False, directed=True)
net.from_nx(subgraph)

# ✅ PyVis 2024+ fix: safely edit node properties
for node in net.nodes:
    node_id = node["id"]
    data = subgraph.nodes[node_id]
    node["title"] = f"<b>{data.get('label', 'Unknown')}</b><br>Type: {data.get('type', 'N/A')}"
    node["label"] = data.get("label", "Unknown")

    # Color by node type
    if data.get("type") == "gene":
        node["color"] = "#4CAF50"
    elif data.get("type") == "compound":
        node["color"] = "#2196F3"
    else:
        node["color"] = "#FFC107"

html_path = "pathway_network.html"
net.save_graph(html_path)

# --- Display in Streamlit ---
with open(html_path, "r", encoding="utf-8") as f:
    html = f.read()

st.components.v1.html(html, height=700, scrolling=True)

# --- Optional Biomarker Info Section ---
st.subheader("ℹ️ Biomarker Information")
for bio in selected_biomarkers:
    try:
        data = REST.kegg_get(f"hsa:{bio}").read()
        if "DEFINITION" in data:
            desc = data.split("DEFINITION")[1].split("PATHWAY")[0].strip()
            st.markdown(f"**{bio}** — {desc}")
        else:
            st.markdown(f"**{bio}** — description unavailable.")
    except Exception as e:
        st.markdown(f"**{bio}** — info not available ({e})")
