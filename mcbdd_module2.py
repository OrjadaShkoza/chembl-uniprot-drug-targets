import requests
from chembl_webresource_client.new_client import new_client
from time import sleep
import csv
from tqdm import tqdm


def get_approved_drugs():
    """Retrieve all approved drugs from ChEMBL, sorted by approval year and name"""
    try:
        molecule = new_client.molecule
        approved_drugs = molecule.filter(max_phase=4).order_by('first_approval', 'pref_name')
        return list(approved_drugs)  # Convert to list to ensure we can iterate multiple times
    except Exception as e:
        print(f"Error retrieving approved drugs: {e}")
        return []


def get_drug_targets(drug):
    """Retrieve UniProt accession numbers for targets of a drug"""
    try:
        drug_chembl_id = drug['molecule_chembl_id']
        mechanism = new_client.mechanism.filter(molecule_chembl_id=drug_chembl_id)
        targets = set()

        for mech in mechanism:
            if mech['target_chembl_id']:
                target = new_client.target.filter(target_chembl_id=mech['target_chembl_id']).only('target_components')
                for t in target:
                    for component in t['target_components']:
                        if 'accession' in component and component['accession']:
                            targets.add(component['accession'])

        return [t for t in targets if t is not None and not t.startswith('ENSG')]
    except Exception as e:
        print(f"Error getting targets for drug {drug.get('pref_name', 'Unknown')}: {e}")
        return []


def get_uniprot_keywords(uniprot_id):
    """Enhanced UniProt keyword retrieval with multiple fallback methods"""
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    try:
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=20)
        response.raise_for_status()
        data = response.json()

        keywords = set()

        # 1. Standard keywords section
        if "keywords" in data:
            for kw in data["keywords"]:
                if "name" in kw:
                    keywords.add(kw["name"])

        # 2. Protein description analysis
        if "protein" in data:
            protein_data = data["protein"]
            if "recommendedName" in protein_data and "fullName" in protein_data["recommendedName"]:
                desc = protein_data["recommendedName"]["fullName"]["value"].lower()
                if "receptor" in desc:
                    keywords.add("Receptor")
                if "enzyme" in desc:
                    keywords.add("Enzyme")
                if "channel" in desc:
                    keywords.add("Channel")
                if "transporter" in desc:
                    keywords.add("Transporter")

        # 3. Gene ontology terms
        if "dbReferences" in data:
            for ref in data["dbReferences"]:
                if ref["type"] == "GO":
                    go_term = ref.get("properties", {}).get("term", "")
                    if go_term.startswith("F:"):  # Molecular function
                        keywords.add(go_term[2:].split(";")[0].strip())

        # 4. Protein families from comments
        if "comments" in data:
            for comment in data["comments"]:
                if comment["type"] == "similarity":
                    text = comment["text"][0]["value"].lower()
                    if "belongs to the" in text or "protein family" in text:
                        family = comment["text"][0]["value"].split(".")[0]
                        keywords.add(family)

        # Convert to list and clean
        keywords = [k.strip() for k in keywords if k and k.strip()]

        return keywords if keywords else ["NO_KEYWORDS_FOUND"]

    except requests.exceptions.HTTPError as e:
        if response.status_code == 404:
            return ["UNIPROT_ENTRY_NOT_FOUND"]
        return [f"HTTP_ERROR_{response.status_code}"]
    except requests.exceptions.RequestException as e:
        return ["REQUEST_ERROR"]
    except Exception as e:
        return ["PROCESSING_ERROR"]
    finally:
        sleep(0.2)  # Slightly increased delay


def main():
    try:
        # Initialize CSV writers
        with open('drug_targets.csv', 'w', newline='', encoding='utf-8') as dt_file, \
                open('target_keywords.csv', 'w', newline='', encoding='utf-8') as tk_file:

            dt_writer = csv.writer(dt_file)
            tk_writer = csv.writer(tk_file)

            # Write headers
            dt_writer.writerow(['Drug Name', 'Approval Year', 'Targets'])
            tk_writer.writerow(['UniProt ID', 'Keywords'])

            # Step 1: Retrieve all approved drugs
            print("Retrieving approved drugs from ChEMBL...")
            approved_drugs = get_approved_drugs()

            if not approved_drugs:
                print("No approved drugs found or error retrieving data")
                return

            # Step 2: Process drugs approved since 2019
            recent_drugs = []
            for drug in approved_drugs:
                try:
                    if drug.get('first_approval') and int(drug['first_approval']) >= 2019:
                        recent_drugs.append(drug)
                except (ValueError, TypeError):
                    continue

            print(f"\nFound {len(recent_drugs)} drugs approved since 2019")

            if not recent_drugs:
                print("No recent drugs found")
                return

            # Track all unique targets
            all_targets = set()
            drug_target_mapping = []

            # Process each drug with progress bar
            print("\nProcessing drugs and collecting targets...")
            for drug in tqdm(recent_drugs, desc="Processing drugs"):
                try:
                    drug_name = drug.get('pref_name', 'Unknown')
                    approval_year = drug.get('first_approval', 'Unknown')
                    targets = get_drug_targets(drug)

                    # Store drug-target mapping
                    drug_target_mapping.append((drug_name, approval_year, targets))

                    # Add targets to our tracking set
                    all_targets.update(targets)
                except Exception as e:
                    print(f"\nError processing drug {drug.get('pref_name', 'Unknown')}: {e}")
                    continue

            # Write drug-target relationships to CSV
            for drug_name, approval_year, targets in drug_target_mapping:
                dt_writer.writerow([drug_name, approval_year, ', '.join(targets)])

            # Process all targets for keywords
            print(f"\nFetching keywords for {len(all_targets)} unique targets...")

            # Convert to list for tqdm
            targets_list = list(all_targets)

            for target in tqdm(targets_list, desc="Fetching keywords"):
                try:
                    keywords = get_uniprot_keywords(target)
                    # Write to CSV - empty string if no keywords
                    tk_writer.writerow([target, ', '.join(keywords) if keywords else 'NO_KEYWORDS_FOUND'])
                except Exception as e:
                    print(f"\nError processing target {target}: {e}")
                    tk_writer.writerow([target, 'ERROR'])

        print("\nProcessing complete!")
        print("Results saved to:")
        print("- drug_targets.csv (Drug to target mappings)")
        print("- target_keywords.csv (Target to keyword mappings)")

    except Exception as e:
        print(f"\nFatal error in main execution: {e}")


if __name__ == "__main__":
    main()