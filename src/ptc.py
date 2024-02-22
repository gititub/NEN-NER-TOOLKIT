
import argparse
from Bio import Entrez
import json
import pandas as pd
import requests
import time

class Pubtator():

    @staticmethod
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
                [id, relation_id, score, role1['type'], role1['identifier'], role2['type'], role2['identifier'], relation_type])

        return pd.DataFrame(relations_data,
                            columns=['ID', 'relation_id', 'score', 'role1_type', 'role1_id', 'role2_type', 'role2_id', 'type'])

    @staticmethod
    def extract_pubtator(ids, output_format):
        """Python function to extract full-text Pubtator results from a list of PMCIDs or PubMedIDs.
        Returns results in biocjson format if output = "biocjson" and returns results
        in pandas DataFrame format if output = "df".
        """
        print("Extracting PubTator results ...")
        list_of_pubtators = []
        error_sum = 0
        id_sum = 0
        all_relations_df_list = []

        for id in ids:
            print(f"Processing ID: {id}")
            if id.startswith('PMC'):
                url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/pmc_export/biocjson?pmcids={id}"
            else:
                url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson?pmids={id}&full=true"
            try:
                response = requests.get(url)
                response.raise_for_status()
                json_data = response.json()
                if not json_data.get('passages'):
                    print(f"No PubTator results found for ID: {id}")
                    continue
            except requests.exceptions.RequestException as e:
                print(f"Error for {id}: {e}")
                error_sum += 1
            except ValueError as e:
                print("Error parsing JSON data:", e)
                print("URL:", url)
                print("Response:", response.content)
                error_sum += 1
            else:
                list_of_pubtators.append(json_data)
                id_sum += 1
                # Extract relations DataFrame and append to the list
                relations_df = Pubtator.extract_relations(json_data, id)
                all_relations_df_list.append(relations_df)

        print(f"Total number of errors extracting full-text Pubtator data from PubMedIDs: {error_sum}")
        print(f"Total number of articles annotated: {id_sum}")

        df_list = []
        for pubtator in list_of_pubtators:
            data = []
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

                for annotation in d['annotations']:
                    entity_type = annotation['infons']['type']
                    subtype = annotation['infons'].get('subtype')
                    #tmVar = annotation['infons'].get('originalIdentifier')
                    #identifiers = annotation['infons'].get('identifiers')
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

            df = df[df['identifier'].notna()]
            df_list.append(df)

        try:
            all_relations_df = pd.concat(all_relations_df_list, ignore_index=True)
            all_relations_df[['correlation_type', 'entity1', 'entity2']] = all_relations_df['type'].str.split('|', expand=True)
            all_relations_df = all_relations_df.drop('type', axis=1)
        except:
            pass

        if output_format == 'biocjson':
            return list_of_pubtators, id_sum, error_sum
        elif output_format == 'df':
            if not df_list:
                print("No articles with annotations found.")
                return pd.DataFrame(), id_sum, error_sum
            else:
                merged_df = pd.concat(df_list, ignore_index=True)
                return merged_df, all_relations_df, id_sum, error_sum
        else:
            return None, None, id_sum, error_sum

    @staticmethod
    def query(query, retmax, pub_date=None):
        email = "weliw001@hotmail.com"
        Entrez.email = email
        if pub_date:
            query += f" AND {pub_date}[Date - Publication]"
        handle = Entrez.esearch(db="pmc", term=query, retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        return ids

    @staticmethod
    def run_extraction(ids, output_file):
        if output_file.endswith(".json"):
            output_data, id_sum, error_sum = Pubtator.extract_pubtator(ids, 'biocjson')
            with open(output_file, 'w') as output_file:
                json.dump(output_data, output_file, indent=2)
            print(f"Biocjson data saved to {output_file}")
        elif output_file.endswith(".tsv"):
            merged_df, all_relations_df, id_sum, error_sum = Pubtator.extract_pubtator(ids, 'df')
            merged_df.to_csv(output_file, sep='\t', index=False)
            print(f"DataFrame saved to {output_file}")
            relations_filename = output_file.replace(".tsv", "_relations.tsv")
            all_relations_df.to_csv(relations_filename, sep='\t', index=False)
            print(f"All Relations DataFrames saved to {relations_filename}")
        else:
            print("Invalid output format. Please choose 'biocjson' or 'df'.")

        return id_sum, error_sum

    @staticmethod
    def parse_arguments():
        parser = argparse.ArgumentParser(description="Extract PubTator data")

        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument("-i", "--input_file", help="Path to the input file containing PMCIDs or PubedIDs")
        input_group.add_argument("-q", "--query", help="PubMed query")

        parser.add_argument("-o", "--output_file", help="Output filename", required=True)
        parser.add_argument("-max", "--retmax", type=int, help="Maximum number of articles to retrieve", required=True)
        parser.add_argument("-date", "--pub_date", help="Publication date (YYYY/MM/DD)")

        args = parser.parse_args()
        return args

def main():
    # Parse command-line arguments
    args = Pubtator.parse_arguments()

    # Read IDs from the input file if provided
    if args.input_file:
        with open(args.input_file, 'r') as f:
            ids = [line.strip() for line in f]
    elif args.query:
        ids = Pubtator.query(args.query, args.retmax, args.pub_date)

    start_time = time.time()
    # Run extraction
    id_sum, error_sum = Pubtator.run_extraction(ids, args.output_file)

    end_time = time.time()
    total_time = end_time - start_time

    print(f"Total time taken to extract {id_sum} annotations with PTC: {total_time} seconds")
    print(f"Total number of errors encountered: {error_sum}")

if __name__ == "__main__":
    main()
