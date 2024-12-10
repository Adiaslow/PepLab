import requests
from bs4 import BeautifulSoup
import pandas as pd
import time
from pathlib import Path
from tqdm import tqdm

def scrape_cyclicpepedia():
    base_url = "https://www.biosino.org/iMAC/cyclicpepedia"
    peptides = []

    print("Starting scraping process...")
    page = 1

    while True:
        print(f"\nProcessing page {page}")
        response = requests.get(f"{base_url}/list", params={"page": page} if page > 1 else None)
        soup = BeautifulSoup(response.text, 'html.parser')

        rows = soup.find('tbody').find_all('tr') if soup.find('tbody') else []
        if not rows:
            break

        for row in tqdm(rows, desc=f"Extracting peptides from page {page}"):
            try:
                # Get ID from first column
                peptide_link = row.find('td').find('a')
                if not peptide_link:
                    continue

                peptide_id = peptide_link.text.strip()
                print(f"\nProcessing peptide {peptide_id}")

                detail_url = f"{base_url}/detail?id={peptide_id}"
                detail_response = requests.get(detail_url)
                detail_soup = BeautifulSoup(detail_response.text, 'html.parser')

                # Extract data from tables
                peptide_data = {
                    'CPKB_ID': peptide_id,
                    'Family': None,
                    'Function': None,
                    'Molecular_Formula': None,
                    'SMILES': None,
                    'Amino_acid_chain': None
                }

                tables = detail_soup.find_all('table', class_='table')
                for table in tables:
                    for row in table.find_all('tr'):
                        cells = row.find_all(['th', 'td'])
                        if len(cells) >= 2:
                            header = cells[0].get_text(strip=True)
                            value = cells[1].get_text(strip=True)

                            if 'Family' in header:
                                peptide_data['Family'] = value
                            elif 'Function' in header:
                                peptide_data['Function'] = value
                            elif 'Molecular Formula' in header:
                                peptide_data['Molecular_Formula'] = value
                            elif 'SMILES' in header:
                                peptide_data['SMILES'] = value
                            elif 'amino acid chain' in header.lower():
                                peptide_data['Amino_acid_chain'] = value

                # Check for 3D structure file
                structure_link = detail_soup.find('a', href=lambda x: x and '.pdb' in str(x).lower())
                if structure_link:
                    print("Downloading structure...")
                    structure_url = base_url + structure_link['href']
                    structure_response = requests.get(structure_url)
                    if structure_response.ok:
                        Path("structures").mkdir(exist_ok=True)
                        with open(f"structures/{peptide_id}.pdb", 'wb') as f:
                            f.write(structure_response.content)
                        peptide_data['structure_file'] = f"structures/{peptide_id}.pdb"

                peptides.append(peptide_data)
                time.sleep(1)

            except Exception as e:
                print(f"Error processing {peptide_id}: {str(e)}")
                continue

        page += 1

    print(f"\nScraped {len(peptides)} peptides from {page-1} pages")
    df = pd.DataFrame(peptides)
    df.to_csv('cyclic_peptides.csv', index=False)
    return df

if __name__ == "__main__":
    df = scrape_cyclicpepedia()
