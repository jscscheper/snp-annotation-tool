# SNP annotation tool

![Python](https://img.shields.io/badge/Python-3.7%2B-blue?logo=python)
![Biopython](https://img.shields.io/badge/Biopython-1.79-green)
![Streamlit](https://img.shields.io/badge/Streamlit-App-FF4B4B?logo=streamlit)
![License](https://img.shields.io/badge/License-MIT-lightgrey)

A tool for introducing and analyzing Single Nucleotide Polymorphisms (SNPs) in DNA sequences. It provides conservation scoring against Multiple Sequence Alignments (MSA) to predict the potential impact of mutations.

Available as both a Command line interface (CLI) and a Web application.

## Table of contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Web application](#web-application)
  - [Command line](#command-line)
- [Contact](#contact)

---

## Overview

The main goal of this tool is to introduce a single nucleotide polymorphism (SNP) at a specific position in a gene sequence. The tool evaluates the effect of the introduced SNP by measuring conservation against a Multiple Sequence Alignment (MSA) of evolutionarily related protein sequences.

Scoring verdicts:
- Deleterious: <90% conservation
- No effect: ≥90% conservation

Example data files are provided for the TP53 gene in the `data/` directory (`sequence.FASTA` and `msa.clustal`).

## Features

- Single SNP annotation: Analyze specific mutations.
- Batch processing: Process multiple mutations via CSV (CLI only).
- Conservation scoring: Evaluate impact based on evolutionary conservation.
- Dual interface: User-friendly web app and scriptable CLI.
- Format support: Supports FASTA (sequence) and Clustal (MSA) formats.

---

## Installation

1. Clone the repository
```bash
git clone https://github.com/jscscheper/snp-annotater.git
cd snp-annotater
```

2. Set up environment (optional but recommended)
```bash
python -m venv venv
# Windows
venv\Scripts\activate
# Linux/Mac
source venv/bin/activate
```

3. Install dependencies
```bash
pip install -r requirements.txt
```

---

## Usage

### Web application
For an interactive experience, use the Streamlit web app.

```bash
streamlit run app.py
```
This will launch the tool in your default web browser. You can upload your own FASTA and Clustal files or use the provided examples.

### Command line
Run the tool directly from your terminal.

General syntax:
```bash
python -m snp_annotator.cli -n {A, C, G, T} -p <POSITION> -m <MSA_FILE> -s <SEQ_FILE>
```

Arguments:
- `-n, --nucleotide`: The new nucleotide (A, C, G, T).
- `-p, --position`: 1-based position of the SNP.
- `-m, --msa`: Path to the MSA file (Clustal format).
- `-s, --dna_sequence`: Path to the DNA sequence file (FASTA format).
- `--batch`: Path to a CSV file for batch processing (alternative to `-n` / `-p`).

Example:
```bash
python -m snp_annotator.cli -n A -p 1 -m data/msa.clustal -s data/sequence.FASTA
```

---

## Contact

Author: Jamie Scheper  
Email: [jscscheper@gmail.com](mailto:jscscheper@gmail.com)
