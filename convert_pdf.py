import os
import argparse
import time
import warnings
import shutil
from PyPDF2 import PdfReader
from tqdm import tqdm

warnings.filterwarnings("ignore")

def convert_pdf_to_text(pdf_file_path):
    pdf_reader = PdfReader(pdf_file_path)
    text = ""
    for page in pdf_reader.pages:
        text += page.extract_text()
    return text.strip()

def convert_directory_to_text(input_directory, output_directory, max_length):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    pdf_files = [file for file in os.listdir(input_directory) if file.endswith(".pdf")]
    file_number = 1

    for pdf_file in tqdm(pdf_files, desc="Converting PDFs"):
        pdf_file_path = os.path.join(input_directory, pdf_file)
        pdf_text = convert_pdf_to_text(pdf_file_path)

        while len(pdf_text) > max_length:
            split_text = pdf_text[:max_length]
            pdf_text = pdf_text[max_length:]

            output_file_path = os.path.join(output_directory, f"{pdf_file}({file_number}).txt")
            with open(output_file_path, "w") as output_file:
                output_file.write(split_text)
            file_number += 1

        if len(pdf_text) > 0:
            output_file_path = os.path.join(output_directory, f"{pdf_file}({file_number}).txt")
            with open(output_file_path, "w") as output_file:
                output_file.write(pdf_text)
            file_number += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_directory", help="Path to the input directory containing PDF files")
    parser.add_argument("output_directory", help="Path to the output directory for plain text files")
    args = parser.parse_args()

    input_directory = args.input_directory
    output_directory = args.output_directory
    max_length = 4900

    start_time = time.time()

    convert_directory_to_text(input_directory, output_directory, max_length)

    end_time = time.time()
    total_time = end_time - start_time
    print("Total time taken:", total_time, "seconds")
