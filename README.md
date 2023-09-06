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

$ python bern_extract_pmids.py [file_path_pmids] [output_filename]

```
python bern_extract_pmids.py pmids2.csv output_file.tsv
```
```
python bern_extract_pmids.py pmids.tsv output_file.json
```
## NER and NEN for PMC full-text articles with PubTator

$ python ptc_extract_pmc.py [file_path_pmcs] [output_format] [output_filename]

Run example:
```
python ptc_extract_pmc.py pmcs.txt biocjson output_pmc.json
```
```
python ptc_extract_pmc.py pmcs.txt df output_pmc.tsv
```

## NER and NEN for plain text with BERN2

Plain text is limited to 5000 characters. To expedite the process, you can distribute the texts (each no longer than 5000 characters) into subdirectories within the designated input directory (with no more than 120 files per subdirectory). This process is parallelized to concurrently process the subdirectories. The outcome is acquired in a file bearing the same name as the subdirectory, following the specified format (JSON or TSV) – one result file per subdirectory. 

You can download the full text automatically divided into subdirectories yourself, from a PMC list or from a folder with pdf files.

From PMC list:

$ python download_fulltext_bern.py [filepath_pmcs] [output_dir]

Run example: 

```
python download_pmc_fulltext.py pmcs.txt bern_ft
```
Or using Chromedriver:
```
python download_fulltext_bern.py pmcs.txt bern_ft
```
From PDF files:

$ python convert_pdf.py [pdf_directory] [output_directory]

Run example: 
```
python convert_pdf.py pdf_files bern_ft_pdf
```
Then, you can run BERN2:

$ python bern_extract_ann.py [input_dir] [output_dir] [json/df]

Run example: 
```
python bern_extract_ann.py bern_ft_pdf bern_ft_pdf_results df
```

## NER and NEN from query

This command will search for PMC articles related to a query, for example *biotin* or *Hodgkin&Lymphoma* (it is not case sensitive), using Bio.Entrez (Spaces may be replaced by '&' sign). Retrieves and processes PubTator annotations and save the results in the specified output file in biocjson or tsv format.

Entrez.ESearch searches and retrieves primary IDs and term translations, and optionally retains results for future use in the user’s environment. The fourth command-line argument is the number of IDs to retrieve.

The last command-line argument is the oldest publication date. If you don't want to filter by date, simply omit the --pub_date argument.

**Pubmed abstracts**

$ python ptc_extract_pmids_query.py [query] [output_format] [output_filename] [max retrievals] --pub_date ["YYYY/MM/DD"]

Run example: published after January 1, 2023
```
python ptc_extract_pmids_query.py biotin df output_df.tsv 50 --pub_date "2022/01/01"
```

**PMC articles**

$ python ptc_extract_pmc_query.py [query] [output_format] [output_filename] [max retrievals] --pub_date ["YYYY/MM/DD"]

Run example: 
```
python ptc_extract_pmc_query.py BRAF biocjson output_pmc.json 35
```
```
python ptc_extract_pmc_query.py Hodgkin&Lymphoma df output_lymphoma.tsv 25 --pub_date "2021/01/01"
```
## PMC ID converter

Here, pmids is your input file containing the list of pmids (tsv, csv and txt format allowed), output_directory is the directory where you want to save the output file, and _pmc is the suffix you want to add to the output file's name. This script will read the input pmids, retrieve PMC IDs, and save the output file in the specified directory with the desired name.

$ python pmc_from_pmid.py [pmids] [output_directory] [_pmc]
 
Run example: 
```
python pmc_from_pmid.py pmids.tsv '.' _pmc
```
## NER&NEN-App

<img
  src="https://github.com/gititub/test/blob/main/Screenshot.png"
  alt="Alt text"
  title="NER&NEN app"
  style="display: inline-block; margin: 0 auto; max-width: 30px">

Run NER-App in Linux: 
```
cd app; pip install -r requirements.txt
shiny run --reload
```
You can also run NER-App in Windows


Run NEN-App in Linux: 
```
cd appNEN; pip install -r requirements.txt
shiny run --reload
```

<img
  src="https://github.com/gititub/test/blob/main/Screenshot2.png"
  alt="Alt text"
  title="NER&NEN app"
  style="display: inline-block; margin: 0 auto; max-width: 30px">

## NER with SynVar and LitVar

To run example use test.tsv or test2.csv as input file, or use your own data with 3 columns: pmid, gene, HGVS. Returns two files in the specified output directory, one with LitVar normalization and the second one with SynVar normalization.

$ python normalize.py [input_file] [output_directory] 

```
python normalize.py test2.csv '.'
```
