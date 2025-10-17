import streamlit as st
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
import networkx as nx
from pyvis.network import Network
import os

# --- Streamlit Page Config ---
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

# --- Load or Fetch KEGG Pathway Data ---
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
# --- Improved Node Filtering ---
sub_nodes = [
    n for n, d in G.nodes(data=True)
    if any(bio.lower() in d.get("search_text", "") for bio in selected_biomarkers)
]

if not sub_nodes:
    st.warning("⚠️ No exact biomarker matches found in KGML — showing full pathway instead.")
    subgraph = G
else:
    subgraph = G.subgraph(sub_nodes).copy()
    st.success(f"✅ Found {len(sub_nodes)} nodes related to selected biomarkers.")



# --- Improved Node Filtering ---
def matches_biomarker(label, name, biomarkers):
    text = f"{label} {name}".lower()
    return any(bio.lower() in text for bio in biomarkers)

sub_nodes = [
    n for n, d in G.nodes(data=True)
    if matches_biomarker(d.get("label", ""), d.get("name", ""), selected_biomarkers)
]

if not sub_nodes:
    st.warning("⚠️ No exact biomarker matches found — showing full pathway instead.")
    subgraph = G
else:
    subgraph = G.subgraph(sub_nodes).copy()
    st.success(f"✅ Found {len(sub_nodes)} nodes related to selected biomarkers.")


# --- Build PyVis Network ---
net = Network(
    height="700px", width="100%",
    bgcolor="#111", font_color="white", notebook=False, directed=True
)
net.from_nx(subgraph)

# Color + tooltip customization
for node in net.nodes:
    node_id = node["id"]
    data = subgraph.nodes[node_id]
    node["title"] = f"<b>{data.get('label', 'Unknown')}</b><br>Type: {data.get('type', 'N/A')}"
    node["label"] = data.get("label", "Unknown")

    # Highlight selected biomarkers in red
    if any(bio.lower() in data.get("label", "").lower() for bio in selected_biomarkers):
        node["color"] = "#FF5252"  # red highlight
        node["size"] = 25
    elif data.get("type") == "gene":
        node["color"] = "#4CAF50"
        node["size"] = 18
    elif data.get("type") == "compound":
        node["color"] = "#2196F3"
        node["size"] = 15
    else:
        node["color"] = "#FFC107"
        node["size"] = 12

# --- Save and Display Network ---
html_path = "pathway_network.html"
net.save_graph(html_path)

with open(html_path, "r", encoding="utf-8") as f:
    html = f.read()

st.components.v1.html(html, height=750, scrolling=True)


# --- Biomarker Information Section ---
st.markdown("---")
st.subheader("🧬 Biomarker Information")

for bio in selected_biomarkers:
    with st.expander(f"**{bio}** — click to view KEGG details"):
        try:
            data = REST.kegg_get(f"hsa:{bio}").read()
            if "DEFINITION" in data:
                desc = data.split("DEFINITION")[1].split("PATHWAY")[0].strip()
                st.markdown(f"**Description:** {desc}")

            # Optional: show KEGG link
            st.markdown(f"[🔗 Open in KEGG](https://www.kegg.jp/dbget-bin/www_bget?hsa:{bio})")
        except Exception as e:
            st.markdown(f"❌ Information unavailable ({e})")

# --- Footer ---
st.markdown("<br><center>Built with ❤️ using Streamlit, Biopython, and PyVis</center>", unsafe_allow_html=True)
