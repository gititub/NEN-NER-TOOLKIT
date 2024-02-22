#python src/ptc_extract_pmids.py [file_path_pmids] [output_filename]

import argparse
import json
import requests
import pandas as pd
import time
import nltk

def split_into_sentences(text):
    sentences = nltk.sent_tokenize(text)
    return sentences

# Function to find the sentence containing the specified text
def find_sentence_containing_text(sentences, target_text):
    sentence_list = []
    for sentence in sentences:
        if target_text in sentence:
            sentence_list.append(sentence)
    return sentence_list if sentence_list else None

def extract_relations(pubtator, pmid):
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
            [pmid, relation_id, score, role1['type'], role1['identifier'], role2['type'], role2['identifier'], relation_type])

    return pd.DataFrame(relations_data,
                        columns=['PMID', 'relation_id', 'score', 'role1_type', 'role1_id', 'role2_type', 'role2_id', 'type'])

def extract_pubtator(pmids, output):
    print("Extracting PubTator results from PMIDs...")
    list_of_abstracts_pubtator = []
    error_count = 0
    all_relations_df_list = []

    for pmid in pmids:
        print(f"Processing PMID: {pmid}")
        url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids={pmid}"
        try:
            response = requests.get(url)
            response.raise_for_status()
            json_data = response.json()
            if not json_data.get('passages'):
                print(f"No Pubtator results found for PMID: {pmid}")
                continue
        except requests.exceptions.RequestException as e:
            print(f"Error for {pmid}: {e}")
            error_count += 1
        except ValueError as e:
            print("Error parsing JSON data:", e)
            print("URL:", url)
            print("Response:", response.content)
            error_count += 1
        else:
            list_of_abstracts_pubtator.append(json_data)
            # Extract relations DataFrame and append to the list
            relations_df = extract_relations(json_data, pmid)
            all_relations_df_list.append(relations_df)

    df_list = []
    for pubtator in list_of_abstracts_pubtator:
        pmid = pubtator['id']
        PMC = pubtator['passages'][0]['infons'].get('article-id_pmc')
        data = []
        for d in pubtator['passages']:
            text_list = []
            if d['infons'].get('type') == 'title':
                title = d['text']
                title_sentences = split_into_sentences(title)
                text_list.extend(title_sentences)
            elif d['infons'].get('type') == 'abstract':
                abstract = d['text']
                abstract_sentences = split_into_sentences(abstract)
                text_list.extend(abstract_sentences)
            section_type = d['infons']['type']
            for annotation in d['annotations']:
                entity_type = annotation['infons']['type']
                subtype = annotation['infons'].get('subtype')
                #tmVar = annotation['infons'].get('originalIdentifier')
                #identifiers = annotation['infons'].get('identifiers')
                offset = annotation['locations'][0]['offset']
                end = offset + annotation['locations'][0]['length']
                identifier = annotation['infons']['identifier']
                string_text = annotation['text']
                ncbi_homologene = annotation['infons'].get('ncbi_homologene')
                database = annotation['infons'].get('database')
                # Additional data for 'normalized' information
                normalized = annotation['infons'].get('normalized')
                normalized_id = annotation['infons'].get('normalized_id')
                biotype = annotation['infons'].get('biotype')
                name = annotation['infons'].get('name')
                accession = annotation['infons'].get('accession')
                # Extract the sentence containing the string text
                sentence_containing_text = find_sentence_containing_text(text_list, string_text)

                data.append([
                    text_list, pmid, PMC, section_type, string_text, offset, end, entity_type, subtype,
                    identifier, ncbi_homologene, database, normalized,
                    normalized_id, biotype, name, accession, sentence_containing_text
                ])

        df = pd.DataFrame(data, columns=['sentences', 'PMID', 'PMC', 'section_type', 'string_text', 'offset', 'end',
                                         'entity_type', 'entity_subtype', 'identifier', 'ncbi_homologene',
                                         'database', 'normalized', 'normalized_id', 'biotype', 'name', 'accession',
                                         'sentence'])

        df = df[df['identifier'].notna()]
        df_list.append(df)

    all_relations_df = pd.concat(all_relations_df_list, ignore_index=True)

    if output == 'biocjson':
        return list_of_abstracts_pubtator, all_relations_df, error_count
    elif output == 'df':
        if not df_list:
            print("No abstracts with annotations found.")
            return pd.DataFrame(), all_relations_df, error_count
        else:
            merged_df = pd.concat(df_list, ignore_index=True)
            return merged_df, all_relations_df, error_count
    else:
        print("Invalid output format. Please choose 'biocjson' or 'df'.")
        return None, all_relations_df, error_count


# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Path to the input file containing PMIDs")
parser.add_argument("output_filename", help="Output filename")
args = parser.parse_args()

# Set input arguments
input_file = args.input_file
output_filename = args.output_filename

# Read PMIDs from the input file
with open(input_file, 'r') as f:
    pmids = [line.strip() for line in f]

start_time = time.time()

if output_filename.endswith(".json"):
    list_of_abstracts_pubtator, all_relations_df, error_count = extract_pubtator(pmids, 'biocjson')
    with open(output_filename, 'w') as output_file:
        json.dump(list_of_abstracts_pubtator, output_file, indent=2)
    print(f"Biocjson data saved to {output_filename}")
elif output_filename.endswith(".tsv"):
    merged_df, all_relations_df, error_count = extract_pubtator(pmids, 'df')
    merged_df.to_csv(output_filename, sep='\t', index=False)
    print(f"DataFrame saved to {output_filename}")
    relations_filename = output_filename.replace(".tsv", "_relations.tsv")
    all_relations_df.to_csv(relations_filename, sep='\t', index=False)
    print(f"All Relations DataFrames saved to {relations_filename}")
else:
    print("Invalid output format. Please choose 'biocjson' or 'df'.")

end_time = time.time()
total_time = end_time - start_time

print("Total time taken to extract annotations with PTC:", total_time, "seconds")
print("Total PMIDs annotated:", len(pmids))
print("Total errors encountered:", error_count)

