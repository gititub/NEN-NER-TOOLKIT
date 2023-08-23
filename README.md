# test
This repository provides the code for automatic system to retrieve annotations of biomedical concepts such as genes and mutations in PubMed abstracts, PMC full-text articles or plain text.


## Installation

```
git clone https://github.com/gititub/test.git
cd test; pip install -r requirements.txt
```

## NER and NEN for PubMed abstracts with PubTator

$ python ptc_extract_pmids.py [file_path_pmids] [output_format] [output_filename]

Run example: 
```
python ptc_extract_pmids.py pmids.tsv biocjson output_pmids.json
```
```
python ptc_extract_pmids.py pmids.tsv df output_df.tsv
```
## NER and NEN for PubMed abstracts with BERN2

```
python bern_extract_pmids.py list_of_pmids.txt output_file.tsv
```
```
python bern_extract_pmids.py list_of_pmids.txt output_file.json
```
## NER and NEN for PMC full-text articles with PubTator

$ python ptc_extract_ann.py [file_path_pmcs] [output_format] [output_filename]

Run example:
```
python ptc_extract_pmc.py pmcs.txt biocjson output_pmc.json
```
```
python ptc_extract_pmc.py pmcs.txt df output_pmc.tsv
```

## NER and NEN for plain text with BERN2

$ python bern_extract_ann.py [input_dir] [output_dir] [json/df]

Plain text is limited to 5000 characters. To expedite the process, you can distribute the texts (each no longer than 5000 characters) into subdirectories within the designated input directory (with no more than 120 files per subdirectory). This process is parallelized to concurrently process the subdirectories. The outcome is acquired in a file bearing the same name as the subdirectory, following the specified format (JSON or TSV) – one result file per subdirectory. You can download the full text automatically divided into subdirectories yourself, from a PMC list or from a folder with pdf files.

From PMC list:

$ python download_fulltext_bern.py [filepath_pmcs] [output_dir]

Run example: 
```
python download_fulltext_bern.py pmcs.txt bern_ft
```
From PDF files:

$ python convert_pdf.py [pdf_directory] [output_directory]

Run example: 
```
python convert_pdf.py gkw943.pdf bern_ft
```

## NER and NEN from query

Pubmed abstracts

$ python ptc_extract_pmids_query.py [query] [output_format] [output_filename]

Run example: 
```
python ptc_extract_pmids_query.py 'melanoma' df output_df.tsv
```
PMC articles

#python ptc_extract_pmc_query.py [query] [output_format] [output_filename]

Run example: 
```
python ptc_extract_pmc_query.py 'breast cancer' biocjson output_pmc2.json
```
