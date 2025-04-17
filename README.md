# ChEMBL to UniProt Drug-Target Analysis

## Project Description
This project retrieves approved drugs from ChEMBL, identifies their protein targets, and fetches associated UniProt keywords. The solution addresses the three tasks outlined in the "Offline activities of Module II" assignment.

## Solution Overview
The Python script performs the following operations:
1. Retrieves all approved drugs from ChEMBL (max_phase=4)
2. Filters for drugs approved since 2019
3. Identifies protein targets (UniProt accessions) for each drug
4. Fetches UniProt keywords for each target protein
5. Generates two output CSV files with the results

## Files Included
- `chembl_uniprot_drug_targets.py`: Main Python script
- `drug_targets.csv`: Output file mapping drugs to targets
- `target_keywords.csv`: Output file mapping targets to keywords
- `README.txt`: This documentation file

## Requirements
- Python 3.6+
- Required packages:
  - chembl_webresource_client
  - requests
  - tqdm

## Installation
```bash
pip install chembl_webresource-client requests tqdm
