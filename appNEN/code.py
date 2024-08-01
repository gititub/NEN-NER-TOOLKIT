import requests
import json
import pandas as pd
from bs4 import BeautifulSoup


def get_gene_info_by_gene_number(gene_numbers):
    data = []
    number_list = [num.strip() for num in gene_numbers.split(',') if num.strip()]
    for number in number_list:
        url = f"https://www.ncbi.nlm.nih.gov/gene/{number}"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')
        word_element = soup.find('dd', class_='noline')
        word = word_element.contents[0] if word_element else "Not Found"
        sp_element = soup.find('dd', class_='tax')
        sp = sp_element.find('a').contents[0] if sp_element else "Not Found"
        data.append([number, word, sp, url])

    df = pd.DataFrame(data, columns=['gene', 'gene_name', 'sp', 'url'])
    return df


def get_gene_info_by_gene_name(gene, species=None):
    data = []
    gene_list = [num.strip() for num in gene.split(',') if num.strip()]
    for gene in gene_list: 
        if species is None:
            url = f"https://www.ncbi.nlm.nih.gov/gene/?term={gene}"
        else:
            sp1, sp2 = species.split(' ')
            url = f"https://www.ncbi.nlm.nih.gov/gene/?term={sp1}+{sp2}+{gene}"

        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')

        gene_elements = soup.find_all('td', class_='gene-name-id')
        for gene_element in gene_elements[:5]:
            gene_name = gene_element.a.get_text()
            number_element = gene_element.find_next('span', class_='gene-id')
            gene_number = number_element.get_text() if number_element else "Not Found"
            species_element = gene_element.find_next('td').find_next('em')
            species_name = species_element.get_text() if species_element else "Not Found"

            if species is None:
                data.append([gene_name, species_name, gene_number, url])
            else:
                data.append([gene_name, species_name, gene_number, url])
                break  # Break the loop after finding the first entry for the specified species

    df = pd.DataFrame(data, columns=['gene_name', 'sp', 'id', 'url'])
    return df


def get_gene_info_by_rsid(rsids):
    data = []
    rsid_list = [num.strip() for num in rsids.split(',') if num.strip()]
    for rsid in rsid_list:
        url = f"https://www.ncbi.nlm.nih.gov/snp/{rsid}"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, 'html.parser')
        sp_element = soup.select_one(
            "#main_content > main > div > div.summary-box.usa-grid-full > dl:nth-child(1) > dd.species_name")
        sp = sp_element.get_text() if sp_element else "Species Not Found"
        gene_element = soup.select_one(
            "#main_content > main > div > div.summary-box.usa-grid-full > dl:nth-child(2) > dd:nth-child(4) > span")
        gene = gene_element.get_text() if gene_element else "Gene Not Found"

        data.append([rsid, sp, gene, url])
    df = pd.DataFrame(data, columns=['rsid', 'sp', 'gene:variant type', 'link'])
    return df


def fetch_litvar_data(file):
    file['litvar'] = file['gene'] + " " + file['HGVS']
    query_list = file['litvar'].drop_duplicates().tolist()
    data_dict = {
        'Query': [],
        'rsid': [],
        'gene': [],
        'name': [],
        'hgvs': [],
        'pmids_count': [],
        'data_clinical_significance': []
    }
    for query in query_list:
        url = f'https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query={query}'
        response = requests.get(url)
        data = response.json()

        if data:
            for result in data:
                data_dict['Query'].append(query)
                data_dict['rsid'].append(result.get('rsid', 'N/A'))
                data_dict['gene'].append(result.get('gene', 'N/A'))
                data_dict['name'].append(result.get('name', 'N/A'))
                data_dict['hgvs'].append(result.get('hgvs', 'N/A'))
                data_dict['pmids_count'].append(result.get('pmids_count', 'N/A'))
                data_dict['data_clinical_significance'].append(result.get('data_clinical_significance', 'N/A'))

    litvar = pd.DataFrame(data_dict)
    return litvar


def litvar_data(input):
    query_list = [num.strip() for num in input.split(',') if num.strip()]
    data_dict = {
        'Query': [],
        'rsid': [],
        'gene': [],
        'name': [],
        'hgvs': [],
        'pmids_count': [],
        'data_clinical_significance': []
    }
    for query in query_list:
        url = f'https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/autocomplete/?query={query}'
        response = requests.get(url)
        data = response.json()

        if data:
            for result in data:
                data_dict['Query'].append(query)
                data_dict['rsid'].append(result.get('rsid', 'N/A'))
                data_dict['gene'].append(result.get('gene', 'N/A'))
                data_dict['name'].append(result.get('name', 'N/A'))
                data_dict['hgvs'].append(result.get('hgvs', 'N/A'))
                data_dict['pmids_count'].append(result.get('pmids_count', 'N/A'))
                data_dict['data_clinical_significance'].append(result.get('data_clinical_significance', 'N/A'))

    litvar = pd.DataFrame(data_dict)
    return litvar
