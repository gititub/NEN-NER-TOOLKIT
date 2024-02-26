import os
import argparse
import time
import warnings
from bs4 import BeautifulSoup
import requests
import shutil
import uuid

warnings.filterwarnings("ignore")


def save_text_to_files(pmc_list, max_length, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for pmc in pmc_list:
        pmc = pmc.strip()
        try:
            url = f'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{pmc}/ascii'
            response = requests.get(url)
            print(f'Response status for {pmc}: {response.status_code}')
            if response.status_code == 200:
                content = response.content.decode()
                chunks = []
                while content:
                    last_period_index = content.rfind('.', 0, max_length)
                    if last_period_index == -1:
                        last_period_index = max_length
                    chunk = content[:last_period_index + 1]  # Include the last period in the chunk
                    content = content[last_period_index + 1:].strip()
                    chunks.append(chunk)

                for chunk_number, chunk in enumerate(chunks, start=1):
                    output_file_path = os.path.join(output_directory, f"{pmc}({chunk_number}).txt")
                    with open(output_file_path, 'w', encoding='utf-8') as output_file:
                        output_file.write(chunk)
            else:
                url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/pmc_export/biocjson?pmcids={pmc}"
                response = requests.get(url)
                print(f'Retry: Response status for {pmc}: {response.status_code}')
                json_data = response.json()
                text_content = []
                if not json_data.get('passages'):
                    print(f"No text found for ID: {pmc}")
                    continue
                else:
                    for passage in json_data['passages']:
                        text = passage['text']
                        text_content.append(text)
                    joined_text = '.'.join(text_content)
                    chunks_pubtator3 = []
                    while joined_text:
                        last_period_index = joined_text.rfind('.', 0, max_length)
                        if last_period_index == -1:
                            last_period_index = max_length
                        chunk_pubtator3 = joined_text[:last_period_index + 1]
                        joined_text = joined_text[last_period_index + 1:].strip()
                        chunks_pubtator3.append(chunk_pubtator3)

                    for chunk_number, chunk_pubtator3 in enumerate(chunks_pubtator3, start=1):
                        output_file_path_pubtator3 = os.path.join(output_directory,
                                                                  f"{pmc}_pubtator3({chunk_number}).txt")
                        with open(output_file_path_pubtator3, 'w', encoding='utf-8') as output_file_pubtator3:
                          output_file_pubtator3.write(chunk_pubtator3)

        except Exception as e:
            print("Error occurred while processing {}: {}".format(pmc, e))
            
    file_list = os.listdir(output_directory)
    num_subdirectories = len(file_list) // 120 + 1
    for i in range(num_subdirectories):
         subdirectory_name = f'subdirectory_{uuid.uuid4()}'
         subdirectory_path = os.path.join(output_directory, subdirectory_name)
         os.makedirs(subdirectory_path)
         start_index = i * 120
         end_index = min((i + 1) * 120, len(file_list))
         files_to_move = file_list[start_index:end_index]
         # Move the files to the subdirectory
         for file_name in files_to_move:
             source_path = os.path.join(output_directory, file_name)
             destination_path = os.path.join(subdirectory_path, file_name)
             shutil.move(source_path, destination_path)
         print(f'Subdirectory "{subdirectory_name}" created with {len(files_to_move)} files.')

    print("Text extraction and saving completed.")

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Path to the input list of pmcs")
parser.add_argument("output_directory", help="Path to the output directory")
args = parser.parse_args()

# Set input arguments
input_file = args.input_file
output_directory = args.output_directory

# Read PMCs from the input file
if input_file.startswith('PMC'):
    pmc_list = [num.strip() for num in input_file.split(',') if num.strip()]
else:
    with open(input_file, 'r') as f:
        pmc_list = [line.strip() for line in f]

start_time = time.time()
max_length = 4900

save_text_to_files(pmc_list, max_length, output_directory)

end_time = time.time()
total_time = end_time - start_time
print("Total time taken:", total_time, "seconds")
