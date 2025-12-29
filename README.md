# SNP Annotation Tool

![Python](https://img.shields.io/badge/Python-3.7%2B-blue?logo=python)
![Biopython](https://img.shields.io/badge/Biopython-1.79-green)
![License](https://img.shields.io/badge/License-Educational-lightgrey)

A command-line tool for introducing and analyzing Single Nucleotide Polymorphisms (SNPs) in DNA sequences, with conservation scoring against multiple sequence alignments.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Contact](#contact)

---

## Overview

The main goal of this tool is to introduce a single nucleotide polymorphism (SNP) on a given position in a gene sequence specified by the user. The effect of the introduced SNP is measured in conservation against a multiple sequence alignment (MSA) which includes protein sequences that are evolutionarily closely aligned with the given gene sequence.

The system measures conservation in percentage and gives a final verdict on whether an SNP is:
- **Deleterious**: <90% conservation
- **No effect**: ≥90% conservation

Example data files are provided for the TP53 gene (`data/sequence.FASTA` and `data/msa.clustal`).

---

## Features

- Single SNP annotation
- Batch processing from CSV files
- Conservation scoring against MSA
- Supports FASTA and Clustal formats

---

## Installation


**Step 1: Acquiring Files**

Either [clone][clone] or [download][download] the source files.

**Step 2: Installing Python**

The script is developed in the language Python (version 3.7.9). Please follow instructions to install Python [here][Python].

**Step 3: Installing Necessary Packages**

In addition, a necessary external package needs to be installed; [BioPython][biopython] (version 1.79) is used to read in the FASTA and MSA files. Installation can be done by first setting up a virtual environment, although not required. Thereafter, simply execute the following command on the command-line:

```bash
$ pip install biopython
```

Now you are set to use this script.

### Usage
As stated before, the user is free to use their own files. For that reason, I included one statement for direct execution (copy-paste) and a more general description.

#### General description
```bash
$ python3 snp.py -n {A, C, G, T} -p <POSITION> -m <MSA.CLUSTAL FILE> -s <GENE_SEQUENCE.FASTA FILE>
```

Where `-n` represents the nucleotide—A, C, G or T—that forms the SNP, `-p` the position of the replacement, `-m` the MSA in the form of a Clustal file and `-s` the gene sequence as a FASTA file.

#### Direct execution
```bash
$ python3 snp.py -n A -p 1 -m msa.clustal -s sequence.FASTA
```

### Contact
If any issue or question remains, please contact us at [d.j.scheper@st.hanze.nl](mailto:d.j.scheper@st.hanze.nl)

[msa]: https://www.ncbi.nlm.nih.gov/homologene/41131
[gene]: https://www.ncbi.nlm.nih.gov/nuccore/NM_016399.3
[clone]: https://jscscheper@bitbucket.org/jscscheper/snp_opdracht_bin3.git
[download]: https://bitbucket.org/jscscheper/snp_opdracht_bin3/src/master/
[python]: https://www.python.org/
[biopython]: https://biopython.org/