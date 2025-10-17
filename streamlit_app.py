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
if not os.path.exists("hsa05010.xml"):
    with st.spinner("Fetching KEGG Alzheimer’s pathway..."):
        kgml_data = REST.kegg_get("hsa05010", "kgml").read()
        with open("hsa05010.xml", "w") as f:
            f.write(kgml_data)

pathway = KGML_parser.read(open("hsa05010.xml"))
G = nx.DiGraph()

# --- Build Graph ---
for entry in pathway.entries.values():
    if entry.type in ["gene", "enzyme", "compound"]:
        label = entry.graphics.name if getattr(entry, "graphics", None) and getattr(entry.graphics, "name", None) else entry.name
        G.add_node(entry.id, name=entry.name, label=label, type=entry.type)


for relation in pathway.relations:
    G.add_edge(relation.entry1.id, relation.entry2.id, type=relation.type)

# --- Filter by Biomarkers ---
sub_nodes = [n for n, d in G.nodes(data=True)
             if any(bio in d.get("name", "") for bio in selected_biomarkers)]
subgraph = G.subgraph(sub_nodes)

# --- Visualization ---
net = Network(height="600px", width="100%", bgcolor="#111", font_color="white", notebook=False)
net.from_nx(subgraph)

for n, d in subgraph.nodes(data=True):
    net.nodes[n]['title'] = f"{d.get('label')}<br>Type: {d.get('type')}"
    net.nodes[n]['label'] = d.get('label', 'Unknown')

net.save_graph("pathway_network.html")

# --- Display in Streamlit ---
with open("pathway_network.html", "r", encoding="utf-8") as f:
    html = f.read()
st.components.v1.html(html, height=650, scrolling=True)

# --- Optional Biomarker Info Section ---
st.subheader("ℹ️ Biomarker Information")
for bio in selected_biomarkers:
    try:
        data = REST.kegg_get(f"hsa:{bio}").read()
        desc = data.split("DEFINITION")[1].split("PATHWAY")[0].strip()
        st.markdown(f"**{bio}** — {desc}")
    except:
        st.markdown(f"**{bio}** — information not available.")

