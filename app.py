
import streamlit as st
import tempfile
import os
from snp_annotator.core import SNPAnnotator

st.title("SNP annotation tool")

st.subheader("1. Upload data")
seq_file = st.file_uploader("Gene Sequence (FASTA)", type=["fasta", "txt"])
msa_file = st.file_uploader("MSA (Clustal)", type=["clustal", "aln"])

st.subheader("2. Enter mutation")
nucleotide = st.selectbox("New Nucleotide", ['A', 'C', 'G', 'T'])
position = st.number_input("Position (1-based)", min_value=1, value=1)

def save_and_get_path(uploaded):
    """Helper to save uploaded file to temp."""
    if uploaded is None: return None
    with tempfile.NamedTemporaryFile(delete=False, suffix=f"_{uploaded.name}") as f:
        f.write(uploaded.getvalue())
        return f.name

if st.button("Analyze"):
    if seq_file and msa_file:
        seq_path = save_and_get_path(seq_file)
        msa_path = save_and_get_path(msa_file)
        
        try:
            annotator = SNPAnnotator(msa_path, seq_path)
            mutated = annotator.insert_snp(position, nucleotide)
            amino = annotator.translate_sequence(mutated)
            outcome = annotator.determine_outcome(amino, position)
            verdict = annotator.get_verdict(outcome)
            
            st.success(f"Result: {verdict.upper()}")
            st.write(f"Conservation Score: {outcome}%")
            
        except Exception as e:
            st.error(f"Error: {e}")
        finally:
            if seq_path: os.remove(seq_path)
            if msa_path: os.remove(msa_path)
    else:
        st.error("Please upload both files.")
