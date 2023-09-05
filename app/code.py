import requests
import json
import pandas as pd
from Bio import Entrez

def extract_pubtator(ids, output):
    id_list = [num.strip() for num in ids.split(',') if num.strip()]
    results_json = []
    results = []

    for id in id_list:
        print(f"Processing ID: {id}")
        url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmids=" + str(id)
        try:
            response = requests.get(url)
            response.raise_for_status()
            json_data = response.json()
            results_json.append(json_data)
        except requests.exceptions.RequestException as e:
            print(f"Error for {id}: {e}")
        except ValueError as e:
            print("Error parsing JSON data:", e)
            print("URL:", url)
            print("Response:", response.content)
            continue

    for pubtator in results_json:
        PMID_list = []
        PMC_list = []
        section_type_list = []
        sentences_list = []
        entity_type_list = []
        offset_list = []
        end_list = []
        identifier_list = []
        string_text_list = []
        subtype_list = []
        tmvar_list = []
        identifiers_list = []
        homologene_list = []

        PMID = id  # pubtator['id']
        PMC = pubtator['passages'][0]['infons'].get('article-id_pmc')

        for d in pubtator['passages']:
            section_type = d['infons']['type']
            for annotation in d['annotations']:
                entity_type = annotation['infons']['type']
                subtype = annotation['infons'].get('subtype')
                tmVar = annotation['infons'].get('originalIdentifier')
                identifiers = annotation['infons'].get('identifiers')
                offset = annotation['locations'][0]['offset']
                end = offset + annotation['locations'][0]['length']
                identifier = annotation['infons']['identifier']
                string_text = annotation['text']
                ncbi_homologene = annotation['infons'].get('ncbi_homologene')

                PMID_list.append(PMID)
                PMC_list.append(PMC)
                section_type_list.append(section_type)
                entity_type_list.append(entity_type)
                offset_list.append(offset)
                end_list.append(end)
                identifier_list.append(identifier)
                identifiers_list.append(identifiers)
                string_text_list.append(string_text)
                subtype_list.append(subtype)
                tmvar_list.append(tmVar)
                homologene_list.append(ncbi_homologene)

        df = pd.DataFrame({
            'PMID': PMID_list,
            'PMC': PMC_list,
            'section_type': section_type_list,
            'string_text': string_text_list,
            'offset': offset_list,
            'end': end_list,
            'entity_type': entity_type_list,
            'entity_subtype': subtype_list,
            'ncbi_homologene': homologene_list,
            'tmVar': tmvar_list,
            'identifiers_list': identifiers_list,
            'identifier': identifier_list
        })

        results.append(df)

    if len(results) == 0:
        print("No dataframes to concatenate.")

    if output == 'biocjson':
        return json.dumps(results_json, indent=4)
    elif output == 'df':
        combined_df = pd.concat(results, ignore_index=True)
        return combined_df
    else:
        print("Invalid output format. Please choose 'biocjson' or 'df'.")


def extract_pubtator_from_pmcs(ids, output):
    id_list = [num.strip() for num in ids.split(',') if num.strip()]
    list_of_pubtators = []
    error_ids = []
    results = []

    for id in id_list:
        print(f"Processing ID: {id}")
        url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmcids=" + str(id)
        try:
            response = requests.get(url)
            response.raise_for_status()
            json_data = response.json()
            list_of_pubtators.append(json_data)
        except requests.exceptions.RequestException as e:
            error_ids.append(id)
            continue  # Skip to the next ID if an error occurs
        except ValueError as e:
            print("Error parsing JSON data:", e)
            print("URL:", url)
            print("Response:", response.content)
            continue  # Skip to the next ID if an error occurs

        print(f"Total number of errors extracting Pubtator data from PMCs: {error_ids}")

    for pubtator in list_of_pubtators:
        PMID_list = []
        PMC_list = []
        section_type_list = []
        sentences_list = []
        entity_type_list = []
        offset_list = []
        end_list = []
        identifier_list = []
        string_text_list = []
        subtype_list = []
        tmvar_list = []
        identifiers_list = []
        homologene_list = []

        if pubtator['_id'].split('|')[0] is not None:
            PMID = pubtator['_id'].split('|')[0]
        else:
            PMID = pubtator['passages'][0]['infons'].get('article-id_pmid')

        if pubtator['_id'].split('|')[1] is not None:
            PMC = pubtator['_id'].split('|')[1]
        else:
            PMC = pubtator['passages'][0]['infons'].get('article-id_pmc')

        for d in pubtator['passages']:
            section_type = d['infons']['section_type']

            for annotation in d['annotations']:
                entity_type = annotation['infons']['type']
                subtype = annotation['infons'].get('subtype')
                tmVar = annotation['infons'].get('originalIdentifier')
                identifiers = annotation['infons'].get('identifiers')
                offset = annotation['locations'][0]['offset']
                end = offset + annotation['locations'][0]['length']
                identifier = annotation['infons']['identifier']
                string_text = annotation['text']
                ncbi_homologene = annotation['infons'].get('ncbi_homologene')

                PMID_list.append(PMID)
                PMC_list.append(PMC)
                section_type_list.append(section_type)
                entity_type_list.append(entity_type)
                offset_list.append(offset)
                end_list.append(end)
                identifier_list.append(identifier)
                identifiers_list.append(identifiers)
                string_text_list.append(string_text)
                subtype_list.append(subtype)
                tmvar_list.append(tmVar)
                homologene_list.append(ncbi_homologene)

        df = pd.DataFrame({
            'PMID': PMID_list,
            'PMC': PMC_list,
            'section_type': section_type_list,
            'string_text': string_text_list,
            'offset': offset_list,
            'end': end_list,
            'entity_type': entity_type_list,
            'entity_subtype': subtype_list,
            'ncbi_homologene': homologene_list,
            'tmVar': tmvar_list,
            'identifiers_list': identifiers_list,
            'identifier': identifier_list
        })
        results.append(df)

    if len(results) == 0:
        print("No dataframes to concatenate.")

    if output == 'biocjson':
        return json.dumps(list_of_pubtators, indent=4)
    elif output == 'df':
        combined_df = pd.concat(results, ignore_index=True)
        return combined_df
    else:
        print("Invalid output format. Please choose 'biocjson' or 'df'.")


def query_plain(text, output):
    url = "http://bern2.korea.ac.kr/plain"
    result = requests.post(url, json={'text': text}).json()

    extracted_data = []
    annotations = result.get('annotations', [])
    for annotation in annotations:
        annotation_id = annotation.get('id', [])
        is_neural_normalized = annotation.get('is_neural_normalized')
        prob = annotation.get('prob')
        mention = annotation.get('mention')
        normalized_name = annotation.get('normalizedName')
        mutation_type = annotation.get('mutationType')
        obj = annotation.get('obj')
        span_begin = annotation['span']['begin']
        span_end = annotation['span']['end']

        extracted_item = {
            'id': annotation_id,
            'is_neural_normalized': is_neural_normalized,
            'prob': prob,
            'mention': mention,
            'normalized_name': normalized_name,
            'mutation_type': mutation_type,
            'obj': obj,
            'span_begin': span_begin,
            'span_end': span_end
        }

        extracted_data.append(extracted_item)

    if len(extracted_data) == 0:
        print("No dataframes to concatenate.")
    df = pd.DataFrame(extracted_data)
    df['dbSNP'] = df['normalized_name'].str.extract(r'(?:rs|RS#:)(\d+)', expand=False)
    df['dbSNP'] = 'rs' + df['dbSNP']

    if output == 'biocjson':
        return json.dumps(result, indent=4)
    elif output == 'df':
        return df
    else:
        print("Invalid output format. Please choose 'biocjson' or 'df'.")



def extract_pubtator_from_pmcs_query(query, pub_date, retmax, output):
    email = "weliw001@hotmail.com"
    Entrez.email = email
    if pub_date:
        pub_date_formatted = pub_date.strftime("%Y/%m/%d")
        query += f" AND {pub_date_formatted}[Date - Publication]"
    handle = Entrez.esearch(db="pmc", term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    pmcs = record["IdList"]

    list_of_pubtators = []
    error_pmc_ids = []
    error_sum = 0
    pmc_sum = 0

    for pmc in pmcs:
        print(f"Processing PMC: PMC{pmc}")
        url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmcids=PMC" + str(pmc)
        try:
            response = requests.get(url)
            response.raise_for_status()
            json_data = response.json()
        except requests.exceptions.RequestException as e:
            error_pmc_ids.append(pmc)
            error_sum += 1
        except ValueError as e:
            print("Error parsing JSON data:", e)
            print("URL:", url)
            print("Response:", response.content)
            error_sum += 1
        else:
            list_of_pubtators.append(json_data)
            pmc_sum += 1

    print(f"Total number of errors extracting Pubtator data from PMCs: {error_sum}")
    print(f"Total number of PMC articles annotated: {pmc_sum}")

    df_list = []

    for pubtator in list_of_pubtators:
        PMID_list = []
        PMC_list = []
        section_type_list = []
        sentences_list = []
        entity_type_list = []
        offset_list = []
        end_list = []
        identifier_list = []
        string_text_list = []
        subtype_list = []
        tmvar_list = []
        identifiers_list = []
        homologene_list = []

        if pubtator['_id'].split('|')[0] is not None:
            PMID = pubtator['_id'].split('|')[0]
        else:
            PMID = pubtator['passages'][0]['infons'].get('article-id_pmid')

        if pubtator['_id'].split('|')[1] is not None:
            PMC = pubtator['_id'].split('|')[1]
        else:
            PMC = pubtator['passages'][0]['infons'].get('article-id_pmc')

        for d in pubtator['passages']:
            section_type = d['infons']['section_type']

            for annotation in d['annotations']:
                entity_type = annotation['infons']['type']
                subtype = annotation['infons'].get('subtype')
                tmVar = annotation['infons'].get('originalIdentifier')
                identifiers = annotation['infons'].get('identifiers')
                offset = annotation['locations'][0]['offset']
                end = offset + annotation['locations'][0]['length']
                identifier = annotation['infons']['identifier']
                string_text = annotation['text']
                ncbi_homologene = annotation['infons'].get('ncbi_homologene')

                PMID_list.append(PMID)
                PMC_list.append(PMC)
                section_type_list.append(section_type)
                entity_type_list.append(entity_type)
                offset_list.append(offset)
                end_list.append(end)
                identifier_list.append(identifier)
                identifiers_list.append(identifiers)
                string_text_list.append(string_text)
                subtype_list.append(subtype)
                tmvar_list.append(tmVar)
                homologene_list.append(ncbi_homologene)

        df = pd.DataFrame({
            'PMID': PMID_list,
            'PMC': PMC_list,
            'section_type': section_type_list,
            'string_text': string_text_list,
            'offset': offset_list,
            'end': end_list,
            'entity_type': entity_type_list,
            'entity_subtype': subtype_list,
            'ncbi_homologene': homologene_list,
            'tmVar': tmvar_list,
            'identifiers_list': identifiers_list,
            'identifier': identifier_list
        })

        df = df[df['identifier'].notna()]
        df_list.append(df)

    if len(df_list) == 0:
        print("No dataframes to concatenate.")

    merged_df = pd.concat(df_list, ignore_index=True)

    if output == 'biocjson':
        return json.dumps(list_of_pubtators, indent=4)
    elif output == 'df':
        return merged_df
    else:
        print("Invalid output format. Please choose 'biocjson' or 'df'.")
