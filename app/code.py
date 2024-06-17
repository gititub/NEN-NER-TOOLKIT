import requests

import json
import pandas as pd
from Bio import Entrez
from drug_named_entity_recognition import find_drugs
import spacy
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
from json import JSONEncoder


class SetEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


def extract_relations(pubtator, id):
    relations_list = pubtator.get('relations', [])
    relations_display_list = pubtator.get('relations_display', [])
    relations_data = []

    for relation, relation_display in zip(relations_list, relations_display_list):
        relation_id = relation['id']
        score = relation['infons']['score']
        role1 = relation['infons']['role1']
        role2 = relation['infons']['role2']
        relation_type = relation_display['name']

        relations_data.append(
            [id, relation_id, score, role1['type'], role1['identifier'], role2['type'], role2['identifier'],
             relation_type])

    return pd.DataFrame(relations_data,
                        columns=['ID', 'relation_id', 'score', 'role1_type', 'role1_id', 'role2_type', 'role2_id',
                                 'type'])


def extract_sibils(ids, output_format):
    if not isinstance(ids, list):
        id_list = [num.strip() for num in ids.split(',') if num.strip()]
    else:
        id_list = ids

    print("Extracting Sibils results ...")
    list_of_sibils = []
    error_sum = 0
    id_sum = 0

    for id in id_list:
        print(f"Processing ID: {id}")
        if id.startswith('PMC'):
            url = f"https://sibils.text-analytics.ch/api/fetch?ids={id}&col=pmc"
        else:
            url = f"https://sibils.text-analytics.ch/api/fetch?ids={id}&col=medline"
        try:
            response = requests.post(url)
            response.raise_for_status()
            sibil_data = response.json()
            list_of_sibils.extend(sibil_data.get("sibils_article_set", []))
            id_sum += len(sibil_data.get("sibils_article_set", []))
        except requests.exceptions.RequestException as e:
            print(f"Error for {id}: {e}")
            error_sum += 1

    print(f"Total number of errors extracting Sibils data: {error_sum}")
    print(f"Total number of articles annotated: {id_sum}")

    df_list = []
    for article in list_of_sibils:
        id = article["document"]["_id"]
        # mesh_terms = ", ".join(article["document"]["mesh_terms"]) if article["document"]["mesh_terms"] else ""
        # chemicals = ", ".join([chem.split(',')[1] for chem in article["document"]["chemicals"]]) if article["document"]["chemicals"] else ""
        annotations = []
        for annotation in article.get('annotations', []):
            annotations.append({
                "ID": id,
                "Type": annotation['type'],
                "Concept Source": annotation['concept_source'],
                "Version": annotation['version'],
                "Concept ID": annotation['concept_id'],
                "Concept Form": annotation['concept_form'],
                "Preferred Term": annotation['preferred_term'],
                "Nature": annotation['nature'],
                "Field": annotation['field'],
                "Start Index": annotation['start_index'],
                "End Index": annotation['end_index']
            })

        df_list.append(pd.DataFrame(annotations))

    if output_format == 'biocjson':
        if list_of_sibils:
            return json.dumps(list_of_sibils, indent=4)
        else:
            return f'No results found. Check if the PubMed or PMC ID is correct.'
    elif output_format == 'df':
        if not df_list:
            print("No articles with annotations found.")
            return pd.DataFrame()
        else:
            merged_df = pd.concat(df_list, ignore_index=True)
            return merged_df
    else:
        print("Invalid output format. Please specify 'df' or 'biocjson'.")
        return None


def extract_pubtator(ids, output_format):
    if not isinstance(ids, list):
        id_list = [num.strip() for num in ids.split(',') if num.strip()]
    else:
        id_list = ids
    print("Extracting PubTator results ...")
    list_of_pubtators = []
    error_sum = 0
    id_sum = 0
    all_relations_df_list = []

    for id in id_list:
        print(f"Processing ID: {id}")
        if id.startswith('PMC'):
            url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/pmc_export/biocjson?pmcids={id}"
        else:
            url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids={id}&full=true"
        try:
            response = requests.get(url)
            response.raise_for_status()
            json_data = response.json()

            # Extracting PubTator data from new JSON structure
            for item in json_data["PubTator3"]:
                passages = item.get("passages", [])
                if passages:
                    list_of_pubtators.append(item)
                    id_sum += 1
                    # Extract relations DataFrame and append to the list
                    relations_df = extract_relations(item, id)
                    all_relations_df_list.append(relations_df)

        except requests.exceptions.RequestException as e:
            print(f"Error for {id}: {e}")
            error_sum += 1
        except ValueError as e:
            print("Error parsing JSON data:", e)
            print("URL:", url)
            print("Response:", response.content)
            error_sum += 1

    print(f"Total number of errors extracting Pubtator data: {error_sum}")
    print(f"Total number of articles annotated: {id_sum}")

    df_list = []
    texts_list = []
    for pubtator in list_of_pubtators:
        data = []
        full_text = []
        if pubtator['_id'].split('|')[0] is not None:
            PMID = pubtator['_id'].split('|')[0]
        else:
            PMID = pubtator['passages'][0]['infons'].get('article-id_pmid')

        if pubtator['_id'].split('|')[1] is not None:
            PMC = pubtator['_id'].split('|')[1]
        else:
            PMC = pubtator['passages'][0]['infons'].get('article-id_pmc')

        for d in pubtator['passages']:
            section_type = d['infons'].get('section_type')
            text = d['text']
            full_text.append({'PMID': PMID, 'PMC': PMC, 'Section': section_type, 'Text': text})
            for annotation in d['annotations']:
                entity_type = annotation['infons']['type']
                subtype = annotation['infons'].get('subtype')
                # tmVar = annotation['infons'].get('originalIdentifier')
                # identifiers = annotation['infons'].get('identifiers')
                offset = annotation['locations'][0]['offset']
                end = offset + annotation['locations'][0]['length']
                identifier = annotation['infons'].get('identifier')
                string_text = annotation['text']
                ncbi_homologene = annotation['infons'].get('ncbi_homologene')
                database = annotation['infons'].get('database')
                # Additional data for 'normalized' information
                normalized = annotation['infons'].get('normalized')
                normalized_id = annotation['infons'].get('normalized_id')
                biotype = annotation['infons'].get('biotype')
                name = annotation['infons'].get('name')
                accession = annotation['infons'].get('accession')

                data.append([
                    PMID, PMC, section_type, string_text, offset, end, entity_type, subtype,
                    identifier, ncbi_homologene, database, normalized,
                    normalized_id, biotype, name, accession
                ])

        df = pd.DataFrame(data,
                          columns=['PMID', 'PMC', 'section_type', 'string_text', 'offset', 'end',
                                   'entity_type', 'entity_subtype', 'identifier', 'ncbi_homologene',
                                   'database', 'normalized', 'normalized_id', 'biotype', 'name', 'accession'])

        texts_list.extend(full_text)
        df_full_text = pd.DataFrame(texts_list)
        if not df_full_text.empty:
            df_grouped = df_full_text.groupby(['PMID', 'Section'], sort=False).agg({'Text': ' '.join}).reset_index()
            df_cleaned = df_grouped.dropna(axis=1, how='all')

        df = df[df['identifier'].notna()]
        df_list.append(df)

    try:
        all_relations_df = pd.concat(all_relations_df_list, ignore_index=True)
        all_relations_df[['correlation_type', 'entity1', 'entity2']] = all_relations_df['type'].str.split('|',
                                                                                                          expand=True)
        all_relations_df = all_relations_df.drop('type', axis=1)
    except:
        pass

    if output_format == 'biocjson':
        if list_of_pubtators:
            return json.dumps(list_of_pubtators, indent=4), json.dumps(list_of_pubtators, indent=4), json.dumps(
                texts_list, indent=4)
        else:
            return f'No results found. Check if the PubMed or PMC ID is correct.'
    elif output_format == 'df':
        if not df_list:
            print("No articles with annotations found.")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        else:
            merged_df = pd.concat(df_list, ignore_index=True)
            return merged_df, all_relations_df, df_cleaned


def query(query, retmax, pub_date=None):
    email = "weliw001@hotmail.com"
    Entrez.email = email
    if pub_date:
        pub_date_formatted = pub_date.strftime("%Y/%m/%d")
        query += f" AND {pub_date_formatted}[Date - Publication]"
    handle = Entrez.esearch(db="pmc", term=query, retmax=retmax)  # db="pubmed"
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    return ids


def bern_extract_pmids(pmids, output):
    results = []
    pmid_list = [num.strip() for num in pmids.split(',') if num.strip()]
    json_data = None

    for pmid in pmid_list:
        try:
            json_data, df = process_pmid(pmid)
            if df is not None:
                results.append(df)
        except Exception as e:
            print(f"Error processing PMID {pmid}: {str(e)}")

    if output == 'biocjson':
        if json_data:
            return json.dumps(json_data, indent=4)
        else:
            return f'No results found. Check if the PubMed ID is correct.'
    elif output == 'df':
        if results:
            bern = pd.concat(results, ignore_index=True)
            bern_cleaned = bern.dropna(axis=1, how='all')
            return bern_cleaned


def process_pmid(pmid):
    print(f"Processing PMID {pmid} with BERN2...")
    url = f"http://bern2.korea.ac.kr/pubmed/{pmid}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            json_data = response.json()
            df = json_to_df(json_data)
            df_cleaned = df.dropna(axis=1, how='all')
            pd.set_option('display.max_colwidth', 30)
            if json_data:
                return json_data, df_cleaned
            else:
                return f'No results found. Check if the PubMed ID is correct.'
        else:
            print(f"Request for PMID {pmid} failed with status code:", response.status_code)
    except Exception as e:
        print(f"Error processing PMID {pmid}: {str(e)}")


def json_to_df(json_data):
    pmid_list = []
    id_list = []
    is_neural_normalized_list = []
    prob_list = []
    mention_list = []
    obj_list = []
    begin_list = []
    end_list = []
    norm_list = []
    mut_type_list = []

    # Iterate over the JSON dictionaries and extract the information
    for item in json_data:
        if item['pmid'] is not None:
            pmid = item['pmid']
        else:
            pmid = item.get("pmid")
        annotations = item.get("annotations", [])
        for annotation in annotations:
            annotation_id = annotation["id"][0]
            is_neural_normalized = annotation["is_neural_normalized"]
            prob = annotation.get("prob")
            mention = annotation["mention"]
            mutation_type = annotation.get('mutationType')
            norm = annotation.get('normalizedName')
            obj = annotation["obj"]
            span = annotation.get("span", {})
            begin = span.get("begin")
            end = span.get("end")

            # Append the extracted information to the respective lists
            pmid_list.append(pmid)
            id_list.append(annotation_id)
            is_neural_normalized_list.append(is_neural_normalized)
            prob_list.append(prob)
            mention_list.append(mention)
            obj_list.append(obj)
            begin_list.append(begin)
            end_list.append(end)
            norm_list.append(norm)
            mut_type_list.append(mutation_type)

    # Create a DataFrame using the extracted information
    df = pd.DataFrame({
        "pmid": pmid_list,
        "id": id_list,
        # "is_neural_normalized": is_neural_normalized_list,
        "prob": prob_list,
        "mention": mention_list,
        "normalized name": norm_list,
        "mutation type": mut_type_list,
        "obj": obj_list,
        "span_begin": begin_list,
        "span_end": end_list
    })
    df[['PubChem', 'chEBI', 'DrugBank']] = df.apply(apply_db_from_wikipedia, axis=1, result_type='expand')
    df_cleaned = df.dropna(axis=1, how='all')
    return df_cleaned


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

    if output == 'biocjson':
        if result:
            return json.dumps(result, indent=4)
        else:
            return f'No results found.'
    elif output == 'df':
        df = pd.DataFrame(extracted_data)
        if not df.empty:
            df['dbSNP'] = df['normalized_name'].str.extract(r'(?:rs|RS#:)(\d+)', expand=False)
            df['dbSNP'] = 'rs' + df['dbSNP']
            df[['PubChem', 'chEBI', 'DrugBank']] = df.apply(apply_db_from_wikipedia, axis=1, result_type='expand')
            df_cleaned = df.dropna(axis=1, how='all')
            return df_cleaned


def apply_db_from_wikipedia(row):
    if row['obj'] == 'drug':
        pubchem, chebi, drugbank = db_from_wikipedia(row['mention'])
        return pubchem, chebi, drugbank
    else:
        return None, None, None


def db_from_wikipedia(mention):
    url = f"https://en.wikipedia.org/wiki/{mention}"
    response = requests.get(url)
    pubchem = ""
    chebi = ""
    drugbank = ""

    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        infobox = soup.find('table', class_='infobox')
        if infobox:
            id_elements = infobox.find_all('a', href=True)
            for id_element in id_elements:
                if 'pubchem.ncbi.nlm.nih.gov' in id_element['href']:
                    pubchem = id_element.contents[0]
                elif 'www.drugbank.ca' in id_element['href']:
                    drugbank = id_element.contents[0]
                elif 'www.ebi.ac.uk' in id_element['href']:
                    chebi = id_element.contents[0]

    return pubchem, chebi, drugbank


def count_characters(input_text):
    character_count = len(input_text)
    if character_count > 4990:
        return "You exceeded the character limit"
    else:
        return character_count


def plain_drugs(input_data, output):
    nlp = spacy.blank("en")
    data_dicts = []

    if isinstance(input_data, pd.DataFrame):
        df = input_data
    elif isinstance(input_data, str):
        df = pd.DataFrame({'Text': [input_data]})
    else:
        return f'No results found. Check if your PMC ID or PubMed ID is correct.'
        # raise ValueError("Invalid input data type. Please provide a DataFrame or a plain text string.")

    for index, row in df.iterrows():
        text = row['Text']
        doc = nlp(text)
        drugs = find_drugs([t.text for t in doc], is_ignore_case=True)
        drugs2 = find_drugs(text.split(" "), is_ignore_case=True)

        if drugs:
            for data in drugs:
                drug_dict = {}
                if 'ID' in row:
                    drug_dict['ID'] = row['ID']
                else:
                    drug_dict['ID'] = None
                if 'Section' in row:
                    drug_dict['section'] = row['Section']
                else:
                    drug_dict['section'] = None
                drug_dict['start'] = data[1]
                drug_dict['end'] = data[2]
                drug_info = data[0]

                if 'name' in drug_info:
                    drug_dict['name'] = drug_info['name']
                else:
                    drug_dict['name'] = None
                if 'synonyms' in drug_info:
                    drug_dict['synonyms'] = list(drug_info['synonyms'])
                else:
                    drug_dict['synonyms'] = None
                if 'mesh_id' in drug_info:
                    drug_dict['mesh_id'] = drug_info['mesh_id']
                else:
                    drug_dict['mesh_id'] = None
                if 'drugbank_id' in drug_info:
                    drug_dict['drugbank_id'] = drug_info['drugbank_id']
                else:
                    drug_dict['drugbank_id'] = None
                if 'medline_plus_id' in drug_info:
                    drug_dict['medline_plus_id'] = drug_info['medline_plus_id']
                else:
                    drug_dict['medline_plus_id'] = None
                if 'nhs_url' in drug_info:
                    drug_dict['nhs_url'] = drug_info['nhs_url']
                else:
                    drug_dict['nhs_url'] = None

                data_dicts.append(drug_dict)

    if output == 'biocjson':
        return json.dumps(data_dicts, indent=4, cls=SetEncoder)

    elif output == 'df':
        df = pd.DataFrame(data_dicts)
        if not df.empty:
            df['PubChem'], df['chEBI'], df['DrugBank'] = zip(*df['name'].apply(db_from_wikipedia))
            df['synonyms'] = df['synonyms'].apply(lambda x: ', '.join(x))
            df_cleaned = df.dropna(axis=1, how='all')
            return df_cleaned


def synvar_ann(ids, output):
    id_list = [num.strip() for num in ids.split(',') if num.strip()]
    results_json = []

    data_dict = {
        'pmid': [],
        'genes': [],
        'drugs': []
    }

    for id in id_list:
        print(f"Processing ID: {id}")
        url = f"https://variomes.text-analytics.ch/api/fetchDoc?ids=" + str(id)

        response = requests.get(url)
        if response.ok:
            result = response.json()
            results_json.append(result)

            data_dict['pmid'].append(id)
            gene_data = result.get('publications', [])[0].get('details', {}).get('facet_details', {}).get('genes', [])
            gene_info = [f"{gene.get('preferred_term')}({gene.get('id')})" for gene in gene_data]
            data_dict['genes'].append(gene_info)
            drug_data = result.get('publications', [])[0].get('details', {}).get('facet_details', {}).get('drugs', [])
            drug_info = [f"{drug.get('preferred_term')}({drug.get('id')})" for drug in drug_data]
            data_dict['drugs'].append(drug_info)

    synvar_df = pd.DataFrame(data_dict)

    if output == 'biocjson':
        if results_json:
            return json.dumps(results_json, indent=4)
        else:
            return f'No results found. Check if the PubMed InoD is correct.'
    elif output == 'df':
        if not synvar_df.empty:
            # synvar_df = synvar_df.explode('genes', ignore_index=True)
            # synvar_df = synvar_df.explode('drugs', ignore_index=True)
            df_cleaned = synvar_df.dropna(axis=1, how='all')
            pd.set_option('display.max_colwidth', 30)
            return df_cleaned
