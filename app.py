import streamlit as st
import urllib.request
import urllib.parse
import json
import re
import pandas as pd
from Bio import Entrez, SeqIO

st.set_page_config(page_title="PathoMatch", page_icon="🔬", layout="wide", initial_sidebar_state="expanded")
Entrez.email = "pathomatch.dev@gmail.com"

# --- โค้ด CSS ปรับแต่งหน้าตาเว็บ ---
st.markdown("""
    <style>
        div[data-testid="InputInstructions"] { display: none !important; }
        div[data-baseweb="input"]:focus-within {
            border-color: #3498db !important;
            box-shadow: 0 0 0 1px #3498db !important;
        }
        div[data-testid="stStatusWidget"] { display: none !important; }
    </style>
""", unsafe_allow_html=True)

# --- 1. ระบบจัดการ Session State ---
if 'search_history' not in st.session_state: st.session_state.search_history = []
if 'active_query' not in st.session_state: st.session_state.active_query = ""
if 'user_name' not in st.session_state: st.session_state.user_name = "Researcher"
if 'user_avatar' not in st.session_state: st.session_state.user_avatar = "https://cdn-icons-png.flaticon.com/512/3135/3135715.png"

# --- Session State สำหรับคำนวณ Volume อัตโนมัติ ---
if 'target_vol' not in st.session_state: st.session_state.target_vol = 500
if 'num_plates' not in st.session_state: st.session_state.num_plates = 25
if 'vol_per_plate' not in st.session_state: st.session_state.vol_per_plate = 20

def calc_target_vol():
    st.session_state.target_vol = st.session_state.num_plates * st.session_state.vol_per_plate

# --- SIDEBAR (Navigation & Profile) ---
with st.sidebar:
    profile_container = st.container()
    with st.expander("⚙️ Edit Profile"):
        new_name = st.text_input("Display Name", value=st.session_state.user_name)
        
        avatar_mode = st.radio("Avatar Source", ["Choose Preset", "Upload Custom Image"])
        
        custom_avatar_url = st.session_state.user_avatar
        if avatar_mode == "Choose Preset":
            avatar_choice = st.selectbox("Choose Avatar", ["🔬 Lab Microscope", "🧬 DNA Strand", "🌿 Botanist"])
            if avatar_choice == "🔬 Lab Microscope": custom_avatar_url = "https://cdn-icons-png.flaticon.com/512/3135/3135715.png"
            elif avatar_choice == "🧬 DNA Strand": custom_avatar_url = "https://cdn-icons-png.flaticon.com/512/2026/2026414.png"
            else: custom_avatar_url = "https://cdn-icons-png.flaticon.com/512/2913/2913990.png"
        else:
            uploaded_file = st.file_uploader("Upload Image (PNG, JPG)", type=["png", "jpg", "jpeg"])
            if uploaded_file is not None:
                # แปลงไฟล์อัปโหลดเป็น Base64 หรือใช้สตรีมเพื่อเก็บใน session state
                import base64
                bytes_data = uploaded_file.getvalue()
                encoded = base64.b64encode(bytes_data).decode()
                custom_avatar_url = f"data:image/png;base64,{encoded}"

        if st.button("Update Profile"):
            st.session_state.user_name = new_name
            st.session_state.user_avatar = custom_avatar_url
            st.rerun()

    with profile_container:
        st.image(st.session_state.user_avatar, width=90)
        st.markdown(f"### **{st.session_state.user_name}**")

    st.divider()
    app_mode = st.radio("**Select Workspace:**", ["Bio-Mining Database", "Lab Calculators & Protocols"])
    
    st.divider()
    st.markdown("**PathoMatch System**")
    st.caption("v6.4 (Custom Avatar & Clean UI)")
    st.write("© 2026 PathoMatch Project")

# ==========================================
# ฟังก์ชันและฐานข้อมูลสำหรับระบบวิเคราะห์
# ==========================================
def expand_abbreviation(q):
    q = q.strip()
    generic_map = {
        "F.": "Fusarium", "B.": "Bacillus", "T.": "Trichoderma", 
        "S.": "Streptomyces", "C.": "Colletotrichum", "A.": "Alternaria", "R.": "Rhizoctonia",
        "P.": "Pseudomonas"
    }
    if q.upper() in generic_map: return generic_map[q.upper()]
    match = re.match(r"^([a-zA-Z])\.\s*([a-zA-Z]+)", q)
    if match:
        initial = match.group(1).upper() + "."
        species = match.group(2).lower()
        if initial == "R." and species == "solani": return "Rhizoctonia solani"
        if initial == "F." and species == "solani": return "Fusarium solani"
        if initial == "P." and species in ["infestans", "palmivora", "capsici"]: return f"Phytophthora {species}"
        if initial == "P." and species in ["aeruginosa", "fluorescens", "putida"]: return f"Pseudomonas {species}"
        if initial == "P." and species in ["aphanidermatum", "ultimum"]: return f"Pythium {species}"
        if initial in generic_map: return f"{generic_map[initial]} {species}"
        return f"{initial.replace('.', '')} {species}"
    return q

def get_organism_type(query):
    q = query.lower()
    fungi_list = ["fusarium", "trichoderma", "botrytis", "alternaria", "colletotrichum", "rhizoctonia", "phytophthora", "pythium"]
    gram_pos_list = ["bacillus", "streptomyces"]
    gram_neg_list = ["pseudomonas", "escherichia", "agrobacterium", "xanthomonas", "ralstonia"]
    
    for f in fungi_list:
        if f in q: return "Fungi / Oomycete"
    for gp in gram_pos_list:
        if gp in q: return "Gram-Positive Bacteria"
    for gn in gram_neg_list:
        if gn in q: return "Gram-Negative Bacteria"
    return "Unknown Organism Type"

PROTOCOL_FUNGI = """
1. Grind 50-100 mg of fungal mycelium into fine powder using liquid nitrogen.
2. Add 500 µl of pre-warmed (65°C) 2% CTAB buffer. Mix well.
3. Incubate at 65°C for 30-60 mins, mixing occasionally.
4. Add equal volume (500 µl) of Chloroform:Isoamyl alcohol (24:1). Mix by inversion for 5 mins.
5. Centrifuge at 12,000 x g for 10 mins. Transfer the upper aqueous phase to a new tube.
6. Add 0.6 volumes of ice-cold Isopropanol to precipitate DNA. Incubate at -20°C for 1 hr.
7. Centrifuge at 12,000 x g for 10 mins. Discard supernatant.
8. Wash pellet with 500 µl of 70% Ethanol. Centrifuge, discard ethanol, and air dry pellet.
9. Resuspend DNA in 50 µl TE buffer or sterile dH2O.
"""

PROTOCOL_GRAM_NEG = """
1. Centrifuge 1-2 ml of overnight bacterial culture at 10,000 x g for 2 mins. Discard supernatant.
2. Resuspend the pellet in 200 µl of sterile dH2O or TE buffer.
3. Boil at 95-100°C for 10 minutes to easily lyse the thin cell wall.
4. Snap-cool on ice for 2 minutes.
5. Centrifuge at 13,000 x g for 5 mins to pellet cellular debris.
6. Transfer the clear supernatant (containing genomic DNA) to a new sterile tube.
7. Store at -20°C or use directly as a template for PCR.
"""

PROTOCOL_GRAM_POS = """
1. Centrifuge 1-2 ml of overnight bacterial culture at 10,000 x g for 2 mins. Discard supernatant.
2. Resuspend the pellet in 180 µl of Lysis Buffer containing Lysozyme (e.g., 20 mg/ml) to break down the thick peptidoglycan layer.
3. Incubate at 37°C for 30-60 minutes.
4. Add 20 µl of Proteinase K and 200 µl of AL Buffer (or standard lysis buffer). Mix by vortexing.
5. Incubate at 56°C for 30 minutes.
6. Add 200 µl of absolute Ethanol, mix well, and transfer to a spin column.
7. Centrifuge at 10,000 x g for 1 min. Discard flow-through.
8. Wash with standard wash buffers (AW1 and AW2), then elute DNA with 50-100 µl of Elution Buffer.
"""

# ==========================================
# MODE 2: โหมดคำนวณและเตรียมสาร (Lab Calculators)
# ==========================================
if app_mode == "Lab Calculators & Protocols":
    st.title("Lab Calculators & Protocols")
    st.markdown("Freely create, edit, and scale media preparation recipes and extraction buffers.")
    st.divider()
    
    st.markdown("#### Universal Recipe & Media Builder")
    st.caption("**Instructions:** Select a template or start blank. **Double-click on any cell in the table below to edit chemical names, amounts, or units.** Click the **'+'** icon at the bottom of the table to add new rows!")
    
    template_choice = st.selectbox("Choose a Recipe Template:", [
        "Custom (Blank)", 
        "PDA (Pre-mixed Powder)",
        "PDA (From Scratch)", 
        "LB Broth (Pre-mixed Powder)",
        "LB Broth (From Scratch)", 
        "NA Media (Pre-mixed Powder)",
        "NA Media (From Scratch)", 
        "CTAB Extraction Buffer 2%"
    ])
    
    if "PDA (Pre-mixed" in template_choice:
        default_df = pd.DataFrame([{"Chemical Name": "PDA Pre-mixed Powder", "Amount": 39.0, "Unit": "g"}])
        def_base = 1000
    elif "PDA (From Scratch" in template_choice:
        default_df = pd.DataFrame([{"Chemical Name": "Potato Extract", "Amount": 4.0, "Unit": "g"}, {"Chemical Name": "Dextrose", "Amount": 20.0, "Unit": "g"}, {"Chemical Name": "Agar", "Amount": 15.0, "Unit": "g"}])
        def_base = 1000
    elif "LB Broth (Pre-mixed" in template_choice:
        default_df = pd.DataFrame([{"Chemical Name": "LB Broth Pre-mixed Powder", "Amount": 25.0, "Unit": "g"}])
        def_base = 1000
    elif "LB Broth (From Scratch" in template_choice:
        default_df = pd.DataFrame([{"Chemical Name": "Tryptone", "Amount": 10.0, "Unit": "g"}, {"Chemical Name": "Yeast Extract", "Amount": 5.0, "Unit": "g"}, {"Chemical Name": "NaCl", "Amount": 10.0, "Unit": "g"}])
        def_base = 1000
    elif "NA Media (Pre-mixed" in template_choice:
        default_df = pd.DataFrame([{"Chemical Name": "NA Pre-mixed Powder", "Amount": 28.0, "Unit": "g"}])
        def_base = 1000
    elif "NA Media (From Scratch" in template_choice:
        default_df = pd.DataFrame([{"Chemical Name": "Peptone", "Amount": 5.0, "Unit": "g"}, {"Chemical Name": "Beef Extract", "Amount": 3.0, "Unit": "g"}, {"Chemical Name": "NaCl", "Amount": 5.0, "Unit": "g"}, {"Chemical Name": "Agar", "Amount": 15.0, "Unit": "g"}])
        def_base = 1000
    elif "CTAB" in template_choice:
        default_df = pd.DataFrame([{"Chemical Name": "CTAB powder", "Amount": 2.0, "Unit": "g"}, {"Chemical Name": "NaCl (5M stock)", "Amount": 28.0, "Unit": "ml"}, {"Chemical Name": "Tris-HCl 1M (pH 8.0)", "Amount": 10.0, "Unit": "ml"}, {"Chemical Name": "EDTA 0.5M (pH 8.0)", "Amount": 4.0, "Unit": "ml"}, {"Chemical Name": "PVP (optional)", "Amount": 1.0, "Unit": "g"}])
        def_base = 100
    else:
        default_df = pd.DataFrame([{"Chemical Name": "Ingredient 1", "Amount": 1.0, "Unit": "g"}, {"Chemical Name": "Ingredient 2", "Amount": 50.0, "Unit": "ml"}])
        def_base = 1000

    st.write("---")
    st.markdown("#### Volume Settings")
    
    vcol1, vcol2 = st.columns(2)
    with vcol1:
        base_volume = st.number_input("Recipe Base Volume (ml):", value=def_base, min_value=1)
    with vcol2:
        target_volume = st.number_input("Desired Target Volume (ml):", min_value=1, key="target_vol")
        
    with st.expander("Optional: Calculate Target Volume from Plates"):
        st.caption("Changes here will instantly update the 'Desired Target Volume' above.")
        p_col1, p_col2 = st.columns(2)
        with p_col1:
            st.number_input("Number of Plates:", min_value=1, key="num_plates", on_change=calc_target_vol)
        with p_col2:
            st.number_input("Volume per Plate (ml):", min_value=1, key="vol_per_plate", on_change=calc_target_vol)

    st.write("---")
    st.write("**Edit Recipe Components Table (Double-click to edit):**")
    edited_df = st.data_editor(default_df, num_rows="dynamic", use_container_width=True, key=template_choice)

    st.success(f"**Final Preparation Summary (Total: {target_volume} ml):**")
    
    for index, row in edited_df.iterrows():
        try:
            amt = float(row['Amount'])
            calc = (amt / base_volume) * target_volume
            st.markdown(f"- **{row['Chemical Name']}**: {calc:.2f} {row['Unit']}")
        except:
            st.markdown(f"- *Please enter a valid number for {row.get('Chemical Name', 'this item')}*")
            
    st.markdown(f"- **Distilled Water**: {target_volume} ml")

    st.divider()
    st.markdown("### Standard DNA Extraction Protocols")
    st.caption("Reference these baseline protocols based on your target organism type.")
    
    with st.expander("Fungal DNA Extraction (CTAB Method)"):
        st.markdown(PROTOCOL_FUNGI)
    with st.expander("Gram-Negative Bacteria (Boiling Method)"):
        st.markdown(PROTOCOL_GRAM_NEG)
    with st.expander("Gram-Positive Bacteria (Enzymatic Lysis / Column Method)"):
        st.markdown(PROTOCOL_GRAM_POS)

# ==========================================
# MODE 1: โหมดค้นหาฐานข้อมูล (Bio-Mining Database)
# ==========================================
else:
    st.title("🔬 PathoMatch Database")
    st.markdown("Automated pathogen, biocontrol, and genomic retrieval system")
    st.divider()

    GENUS_SUGGESTIONS = {
        "bacillus": ["Bacillus subtilis", "Bacillus thuringiensis", "Bacillus amyloliquefaciens", "Bacillus cereus"],
        "fusarium": ["Fusarium oxysporum", "Fusarium solani", "Fusarium graminearum"],
        "trichoderma": ["Trichoderma harzianum", "Trichoderma viride", "Trichoderma virens"],
        "phytophthora": ["Phytophthora infestans", "Phytophthora palmivora", "Phytophthora capsici"],
        "pythium": ["Pythium aphanidermatum", "Pythium ultimum"],
        "pseudomonas": ["Pseudomonas fluorescens", "Pseudomonas aeruginosa", "Pseudomonas putida"],
        "rhizoctonia": ["Rhizoctonia solani"],
        "botrytis": ["Botrytis cinerea"],
        "alternaria": ["Alternaria alternata"],
        "colletotrichum": ["Colletotrichum gloeosporioides", "Colletotrichum acutatum"]
    }

    EPITHET_SUGGESTIONS = {
        "solani": ["Fusarium solani", "Rhizoctonia solani"],
        "oxysporum": ["Fusarium oxysporum"],
        "cinerea": ["Botrytis cinerea"],
        "subtilis": ["Bacillus subtilis"],
        "harzianum": ["Trichoderma harzianum"],
        "infestans": ["Phytophthora infestans"],
        "palmivora": ["Phytophthora palmivora"],
        "aphanidermatum": ["Pythium aphanidermatum"],
        "aeruginosa": ["Pseudomonas aeruginosa"],
        "fluorescens": ["Pseudomonas fluorescens"],
        "alternata": ["Alternaria alternata"],
        "gloeosporioides": ["Colletotrichum gloeosporioides"]
    }

    KNOWLEDGE_BASE = {
        "fusarium": {
            "biocontrols": ["Trichoderma", "Bacillus", "Pseudomonas"], 
            "bsl": "BSL-1 or BSL-2 (Plant Pathogen)",
            "prevention": "Crop rotation, soil solarization, planting resistant varieties, and proper field drainage to avoid excess moisture.",
            "toxicity": "Produces dangerous mycotoxins (e.g., fumonisins, deoxynivalenol) which are highly toxic to humans and livestock.",
            "plant_impact": "Highly pathogenic. Causes devastating vascular wilts, root rots, and blights leading to severe crop loss."
        },
        "phytophthora": {
            "biocontrols": ["Trichoderma", "Bacillus"], 
            "bsl": "BSL-1 (Plant Pathogen)",
            "prevention": "Improve soil drainage, avoid over-watering, use disease-free seeds/tubers, and apply appropriate fungicides preventatively.",
            "toxicity": "Generally not toxic to humans or animals.",
            "plant_impact": "Extremely destructive plant pathogen (Oomycete). Causes late blights, root rots, and stem lesions (e.g., potato famine)."
        },
        "trichoderma": {
            "pathogens": ["Fusarium", "Phytophthora", "Rhizoctonia", "Pythium"], 
            "bsl": "BSL-1",
            "prevention": "Maintain proper humidity in mushroom cultivation to prevent Trichoderma green mold infections.",
            "toxicity": "Generally safe for humans and animals (GRAS). Produces enzymes toxic strictly to other fungi.",
            "plant_impact": "Beneficial endophytic symbiont. Promotes root growth, nutrient uptake, and induces systemic resistance (ISR) in plants."
        },
        "bacillus": {
            "pathogens": ["Fusarium", "Rhizoctonia", "Botrytis"], 
            "bsl": "BSL-1 (mostly)",
            "prevention": "For beneficial strains, maintain proper storage. For foodborne strains (like B. cereus), ensure proper food refrigeration.",
            "toxicity": "Most strains are safe (probiotics). Bacillus cereus can cause severe food poisoning in humans.",
            "plant_impact": "Beneficial Plant Growth-Promoting Rhizobacteria (PGPR). Helps in nitrogen fixation and protects against soil-borne pathogens."
        },
        "botrytis": {
            "biocontrols": ["Trichoderma", "Bacillus", "Streptomyces"],
            "bsl": "BSL-1",
            "prevention": "Ensure good air circulation, reduce canopy humidity, and remove decaying plant debris.",
            "toxicity": "May cause allergic reactions (winegrower's lung) upon inhaling high spore loads.",
            "plant_impact": "Causes gray mold (Botrytis blight) on hundreds of plant species, severely impacting flowers and fruits pre/post-harvest."
        }
    }

    CULTURE_KB = {
        "fusarium": "Standard Protocol: Potato Dextrose Agar (PDA) at 25–28 °C. Incubation 5–7 days.",
        "bacillus": "Standard Protocol: Nutrient Agar (NA) or LB broth at 30–37 °C for 18–24 hours.",
        "trichoderma": "Standard Protocol: PDA or Malt Extract Agar (MEA) at 25–30 °C for 5–7 days.",
    }

    with st.form(key='search_form'):
        search_query = st.text_input("🔍 Search microorganism name or abbreviation (e.g., F. oxysporum, B. subtilis):", value=st.session_state.active_query, key="input_query")
        submit_button = st.form_submit_button(label="Search Data")

    clicked_history = None
    if st.session_state.search_history:
        st.write("🕒 **Recent Searches:**")
        cols = st.columns(5) 
        for idx, past_query in enumerate(st.session_state.search_history):
            if cols[idx].button(past_query, key=f"hist_{idx}"): clicked_history = past_query

    if clicked_history:
        search_query = clicked_history
        st.session_state.active_query = clicked_history
    if submit_button:
        st.session_state.active_query = search_query

    target_query = st.session_state.active_query if st.session_state.active_query else search_query

    if submit_button or clicked_history or target_query:
        if target_query:
            target_query = expand_abbreviation(target_query)
            target_lower = target_query.lower()
            parts = target_lower.split()
            genus_part = parts[0] if parts else ""
            species_part = parts[1] if len(parts) > 1 else ""
            
            is_typo_caught = False
            
            if species_part in EPITHET_SUGGESTIONS:
                valid_combos = [s.lower() for s in EPITHET_SUGGESTIONS[species_part]]
                if target_lower not in valid_combos:
                    st.error(f"🚨 **Typo Detected!** '{genus_part.capitalize()}' is unusual for '{species_part}'. Did you mean:")
                    sug_cols = st.columns(len(EPITHET_SUGGESTIONS[species_part]))
                    for i, sug_sp in enumerate(EPITHET_SUGGESTIONS[species_part]):
                        if sug_cols[i].button(sug_sp, key=f"sug_typo_{i}"):
                            st.session_state.active_query = sug_sp
                            st.rerun()
                    st.divider()
                    is_typo_caught = True
            
            if not is_typo_caught and genus_part in GENUS_SUGGESTIONS and target_lower == genus_part:
                st.info(f"💡 **Are you looking for a specific '{genus_part.capitalize()}' species?**")
                sug_cols = st.columns(len(GENUS_SUGGESTIONS[genus_part]))
                for i, sug_sp in enumerate(GENUS_SUGGESTIONS[genus_part]):
                    if sug_cols[i].button(sug_sp, key=f"sug_gen_{i}"):
                        st.session_state.active_query = sug_sp
                        st.rerun()
                st.divider()
            
            if target_query in st.session_state.search_history: st.session_state.search_history.remove(target_query)
            st.session_state.search_history.insert(0, target_query)
            if len(st.session_state.search_history) > 5: st.session_state.search_history = st.session_state.search_history[:5]

            with st.spinner(f"Mining Literature, Taxonomy, and Genomics for '{target_query}'..."):
                
                ncbi_id = "Not Found"
                lineage = "Data not available"
                genus = target_query.split()[0] if target_query else ""
                related_species = []
                
                try:
                    handle = Entrez.esearch(db="taxonomy", term=target_query)
                    record = Entrez.read(handle)
                    handle.close()
                    if record["IdList"]:
                        ncbi_id = record["IdList"][0]
                        fetch_handle = Entrez.efetch(db="taxonomy", id=ncbi_id, retmode="xml")
                        fetch_record = Entrez.read(fetch_handle)
                        fetch_handle.close()
                        if fetch_record:
                            lineage = fetch_record[0].get("Lineage", lineage)
                    
                    if genus:
                        sib_handle = Entrez.esearch(db="taxonomy", term=f"{genus}[Next Level]", retmax=5)
                        sib_record = Entrez.read(sib_handle)
                        sib_handle.close()
                        if sib_record["IdList"]:
                            sib_fetch = Entrez.efetch(db="taxonomy", id=",".join(sib_record["IdList"]), retmode="xml")
                            sib_data = Entrez.read(sib_fetch)
                            sib_fetch.close()
                            related_species = [s["ScientificName"] for s in sib_data if s["ScientificName"].lower() != target_query.lower()]
                except: pass

                fasta_seq = "Sequence data not found or took too long to retrieve."
                try:
                    seq_query = f"{target_query}[Organism] AND (internal transcribed spacer OR 16S ribosomal RNA OR ITS)"
                    seq_search = Entrez.esearch(db="nucleotide", term=seq_query, retmax=1)
                    seq_record = Entrez.read(seq_search)
                    seq_search.close()
                    if seq_record["IdList"]:
                        seq_id = seq_record["IdList"][0]
                        seq_fetch = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
                        fasta_seq = seq_fetch.read()
                        seq_fetch.close()
                except: pass

                combined_text = ""
                wiki_image_url = None
                try:
                    url = f"https://en.wikipedia.org/w/api.php?action=query&prop=extracts|pageimages&exintro=true&explaintext=true&piprop=original|thumbnail&pithumbsize=500&redirects=1&titles={target_query.replace(' ', '%20')}&format=json"
                    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                    res = json.loads(urllib.request.urlopen(req).read())
                    pages = res["query"]["pages"]
                    for page_id in pages:
                        page_data = pages[page_id]
                        combined_text += page_data.get("extract", "") + " "
                        if "original" in page_data: wiki_image_url = page_data["original"]["source"]
                        elif "thumbnail" in page_data: wiki_image_url = page_data["thumbnail"]["source"]
                except Exception: pass

                paper_citations = []
                try:
                    epmc_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={target_query.replace(' ', '%20')}&format=json&resultType=core"
                    req = urllib.request.Request(epmc_url, headers={'User-Agent': 'Mozilla/5.0'})
                    res = json.loads(urllib.request.urlopen(req).read())
                    results = res.get("resultList", {}).get("result", [])
                    for item in results[:10]:
                        combined_text += item.get("abstractText", "") + " "
                        title = item.get("title", "Unknown Title")
                        author = item.get("authorString", "Unknown")
                        year = item.get("pubYear", "N/A")
                        doi = item.get("doi", "")
                        paper_citations.append({"author": author, "year": year, "title": title, "doi": doi})
                except Exception: pass

                disease_info, control_info, general_info, culture_info, host_info = [], [], [], [], []
                toxin_info, prev_info, human_risk_info, plant_impact_info, extract_info = [], [], [], [], []
                detected_biocontrols, detected_pathogens = set(), set()
                
                if combined_text.strip():
                    clean_text = re.sub(r'<[^>]+>', '', combined_text)
                    sentences = re.split(r'(?<=[.!?]) +', clean_text.replace('\n', ' '))
                    
                    disease_kw = ["disease", "pathogen", "rot", "wilt", "infect", "symptom", "blight", "cause"]
                    control_kw = ["control", "antagonist", "fungicide", "inhibit", "biocontrol", "protect", "suppress"]
                    culture_kw = ["agar", "medium", "media", "temperature", "°c", "incubate", "cultivation", "pda", "broth"]
                    host_kw = ["host", "crop", "plant", "isolate", "endophyte", "associated with", "found in"]
                    toxin_kw = ["toxin", "metabolite", "mycotoxin", "antibiotic", "surfactin", "iturin", "compound", "produce"]
                    prev_kw = ["prevent", "quarantine", "sanitation", "crop rotation", "hygiene", "management practice", "steriliz", "disinfect", "avoid"]
                    human_risk_kw = ["human", "animal", "mammal", "health risk", "poison", "allergic", "toxicity", "toxic to"]
                    plant_impact_kw = ["plant pathogen", "phytotoxin", "phytopathogenic", "crop loss", "yield loss", "wilt in", "rot in", "plant disease"]
                    extract_kw = ["dna extraction", "extracted using", "genomic dna", "extraction kit", "ctab method", "dna was extracted", "isolate dna", "dna isolation"]
                    
                    common_biocontrols = ["trichoderma", "bacillus", "pseudomonas", "streptomyces", "gliocladium"]
                    common_pathogens = ["fusarium", "phytophthora", "botrytis", "pythium", "rhizoctonia", "alternaria"]

                    for sentence in sentences:
                        lower_s = sentence.lower()
                        if len(sentence) < 20: continue
                        
                        for bc in common_biocontrols:
                            if bc in lower_s and bc not in target_lower: detected_biocontrols.add(f"🛡️ {bc.capitalize()}")
                        for p in common_pathogens:
                            if p in lower_s and p not in target_lower: detected_pathogens.add(f"🦠 {p.capitalize()}")

                        if any(kw in lower_s for kw in extract_kw):
                            if sentence not in extract_info: extract_info.append(sentence)

                        if any(kw in lower_s for kw in prev_kw):
                            if sentence not in prev_info: prev_info.append(sentence)
                        
                        is_toxin = any(kw in lower_s for kw in toxin_kw)
                        is_human_risk = any(kw in lower_s for kw in human_risk_kw)
                        
                        if is_human_risk and (is_toxin or "disease" in lower_s):
                            if sentence not in human_risk_info: human_risk_info.append(sentence)
                        elif is_toxin:
                            if sentence not in toxin_info: toxin_info.append(sentence)

                        if any(kw in lower_s for kw in plant_impact_kw):
                            if sentence not in plant_impact_info: plant_impact_info.append(sentence)

                        if any(kw in lower_s for kw in host_kw) and sentence not in host_info: host_info.append(sentence)
                        elif any(kw in lower_s for kw in culture_kw) and sentence not in culture_info: culture_info.append(sentence)
                        elif any(kw in lower_s for kw in control_kw) and sentence not in control_info: control_info.append(sentence)
                        elif any(kw in lower_s for kw in disease_kw) and sentence not in disease_info: disease_info.append(sentence)
                        else:
                            if sentence not in general_info: general_info.append(sentence)

                bsl_status = "BSL-1 (Assumed standard for non-human pathogens)"
                for key, proto in CULTURE_KB.items():
                    if key in target_lower: culture_info.insert(0, proto)
                    
                for key, relations in KNOWLEDGE_BASE.items():
                    if key in target_lower:
                        if "biocontrols" in relations:
                            for bc in relations["biocontrols"]: detected_biocontrols.add(f"🛡️ {bc}")
                        if "pathogens" in relations:
                            for p in relations["pathogens"]: detected_pathogens.add(f"🦠 {p}")
                        if "bsl" in relations:
                            bsl_status = relations["bsl"]
                        if "prevention" in relations:
                            prev_info.insert(0, relations["prevention"])
                        if "toxicity" in relations:
                            human_risk_info.insert(0, relations["toxicity"])
                        if "plant_impact" in relations:
                            plant_impact_info.insert(0, relations["plant_impact"])

                report_text = f"PathoMatch Data Report: {target_query}\n"
                report_text += f"NCBI ID: {ncbi_id} | BSL: {bsl_status}\n"
                report_text += f"Lineage: {lineage}\n"
                report_text += "="*50 + "\n"
                if general_info: report_text += "\n[General Info]\n- " + "\n- ".join(general_info[:5])
                if culture_info: report_text += "\n\n[Culture & Protocol]\n- " + "\n- ".join(culture_info[:5])
                if host_info: report_text += "\n\n[Host Range]\n- " + "\n- ".join(host_info[:5])
                report_text += f"\n\n[Reference Sequence]\n{fasta_seq}\n"

                st.divider()
                
                head_col1, head_col2, head_col3 = st.columns([2, 5, 2])
                with head_col1:
                    if wiki_image_url: st.image(wiki_image_url, caption=target_query, use_container_width=True)
                    else: st.info("🖼️ No image available")
                with head_col2:
                    st.markdown(f"## 🦠 **{target_query}**")
                    st.write(f"**NCBI Taxonomy ID:** {ncbi_id} | **Biosafety Level:** {bsl_status}")
                with head_col3:
                    st.download_button(label="📥 Download Report", data=report_text, file_name=f"{target_query.replace(' ', '_')}_report.txt", mime="text/plain")
                
                st.write(f"**🧬 Lineage:** {lineage}")
                    
                st.divider()
                
                tab1, tab2, tab3, tab4, tab5 = st.tabs([
                    "📊 Overview & Disease", 
                    "🧫 Culture & Extraction", 
                    "🌾 Host, Plants & Toxins", 
                    "🧬 Sequences & Relatives",
                    "📚 Literature"
                ])
                
                with tab1:
                    col_t1, col_t2 = st.columns(2)
                    with col_t1:
                        st.info("**General Biological Information:**\n\n" + "\n\n".join(general_info[:5] if general_info else ["No general information found."]))
                    with col_t2:
                        st.warning("⚠️ **Disease & Pathogenicity Overview:**\n\n" + "\n\n".join(disease_info[:5] if disease_info else ["No direct pathogenicity information found."]))
                
                with tab2:
                    st.markdown("#### 🧫 Culture & Growth Conditions:")
                    if culture_info:
                        with st.expander("📌 Growth Media & Incubation Protocol", expanded=False):
                            st.write("\n\n".join(culture_info[:5]))
                    else:
                        st.write("No specific cultivation data found.")
                        
                    st.markdown("#### 🧬 DNA Extraction Protocols")
                    org_type = get_organism_type(target_query)
                    st.write(f"**Organism Classification:** {org_type}")
                    
                    if extract_info:
                        st.success("🔬 **Methods from Recent Literature:**\n\n" + "\n\n".join(extract_info[:3]))
                        
                    if org_type == "Fungi / Oomycete":
                        with st.expander("📍 Target Protocol: Fungal CTAB Method", expanded=True):
                            st.markdown(PROTOCOL_FUNGI)
                    elif org_type == "Gram-Positive Bacteria":
                        with st.expander("📍 Target Protocol: Gram-Positive Bacteria (Enzymatic Lysis)", expanded=True):
                            st.markdown(PROTOCOL_GRAM_POS)
                    elif org_type == "Gram-Negative Bacteria":
                        with st.expander("📍 Target Protocol: Gram-Negative Bacteria (Boiling Method)", expanded=True):
                            st.markdown(PROTOCOL_GRAM_NEG)
                    else:
                        st.info("Protocol cannot be determined automatically. Please check the 'Lab Calculators' mode for full baseline protocols.")
                    
                    st.markdown("#### 🛡️ Control & Prevention Strategies:")
                    if prev_info:
                        with st.expander("🚧 Prevention & Management Practices", expanded=False):
                            st.write("\n\n".join(prev_info[:5]))
                    
                    if control_info:
                        with st.expander("📌 Biocontrol & Chemical Interventions", expanded=False):
                            st.write("\n\n".join(control_info[:5]))
                        
                    all_related_organisms = list(detected_biocontrols) + list(detected_pathogens)
                    if all_related_organisms:
                        st.divider()
                        st.markdown("💡 **Related Organisms :**")
                        bc_cols = st.columns(len(all_related_organisms) + 2) 
                        for idx, org_name in enumerate(all_related_organisms):
                            if bc_cols[idx].button(org_name, key=f"org_btn_{idx}"):
                                actual_name = org_name.split(" ", 1)[1]
                                st.session_state.active_query = actual_name
                                st.rerun()

                with tab3:
                    col_t3, col_t4 = st.columns(2)
                    with col_t3:
                        st.success("🌾 **Host Range & Distribution:**\n\n" + "\n\n".join(host_info[:4] if host_info else ["No specific host range data found."]))
                        
                        if plant_impact_info:
                            st.warning("🍂 **Plant Pathogenicity & Impact:**\n\n" + "\n\n".join(plant_impact_info[:3]))
                        else:
                            st.info("🍂 **Plant Pathogenicity & Impact:**\n\nNo major plant pathogenicity indicators found.")

                    with col_t4:
                        if human_risk_info:
                            st.error("☠️ **Human & Animal Toxicity Risk:**\n\n" + "\n\n".join(human_risk_info[:3]))
                        else:
                            st.info("✅ **Human & Animal Toxicity Risk:**\n\nNo significant direct toxicity to humans or animals reported in the retrieved data.")
                        
                        st.warning("🧪 **Other Toxins & Secondary Metabolites:**\n\n" + "\n\n".join(toxin_info[:3] if toxin_info else ["No general metabolite information found."]))

                with tab4:
                    st.markdown("#### 🧬 Reference Sequence (Marker Gene)")
                    if "Sequence data not found" not in fasta_seq:
                        blast_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&QUERY={urllib.parse.quote(fasta_seq)}"
                        st.markdown(f"[🚀 Run Sequence in NCBI BLAST]({blast_url})", unsafe_allow_html=True)

                    with st.expander("Show FASTA Sequence", expanded=False):
                        st.caption("Automatically fetched from NCBI Nucleotide (ITS for Fungi / 16S for Bacteria).")
                        st.code(fasta_seq, language="fasta")
                    
                    st.divider()
                    st.markdown("#### 🌳 Closely Related Species (Same Genus)")
                    if related_species:
                        st.write("These species share the same evolutionary lineage and may exhibit similar pathogenicity or characteristics:")
                        for sp in related_species:
                            if st.button(f"🔗 {sp}", key=f"rel_{sp}", type="tertiary"):
                                st.session_state.active_query = sp
                                st.rerun()
                    else:
                        st.write("No related species data fetched.")

                with tab5:
                    st.markdown("#### 📚 Recent Scientific Literature (Europe PMC)")
                    st.caption("Auto-mined from the latest publications related to this organism.")
                    if paper_citations:
                        for i, paper in enumerate(paper_citations):
                            doi_link = f"https://doi.org/{paper['doi']}" if paper['doi'] else "#"
                            st.markdown(f"{i+1}. {paper['author']} ({paper['year']}). **{paper['title']}** [🔗 DOI]({doi_link})")
                    else:
                        st.write("No literature citations retrieved for this query.")
