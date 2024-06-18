
This repository provides the code for automatic system to retrieve annotations of biomedical concepts such as genes, mutations and chemicals, from PubMed abstracts, PMC full-text articles, or plain text. The system utilizes API tools such as [PubTator3][def2], [BERN2][def], [SIBiLS][def5], and [Variomes][def4] for Named Entity Recognition (NER), Named Entity Normalization (NEN), and related entities (RE). Additionally, it incorporates [Auto-CORPus][def1] and [GWAS-Miner][def0] that allows retrieval of information about diseases, variants, and significance values (p-values) from full-text articles, including data from inner tables. [LitVar][def6] and [SynVar][def7] are used for NEN. The workflow is integrated into the development of a web application using Shiny for Python.

## Biomedical Entity Recognition, Normalization and Relation Extraction

[Biomedical Entity Recognition, Normalization and Relation Extraction](https://2maxto-amelia-martinez0sequera.shinyapps.io/nerversetoolkit2/)


<img
  src="https://github.com/gititub/test/blob/main/rsc/APP.png"
  alt="Alt text"
  style="display: block; width:400px">


1. Make a query:  
- PMCID or PubMedID to [PubTator][def2] (One or more, comma separated).
  For example, "36064841, PMC9797458" (They can be mixed together).  
- PMCID or PubMedID to [PubTator - Relations][def2] (One or more, comma separated).
  For example, "36064841"
- [PubTator - Relations][def2] from a word query (replace space with ‘&’) + Publication Date + max. Retrievals
- Plain Text to [BERN2][def] (max. 5000 characters)
- PubMedID to [BERN2][def] (one or more, comma separated)
- Plain Text to [Drug Named Entity Recognition][def3]
- PMCID or PubMedID to [Drug Named Entity Recognition][def3] (one or more, comma separated. They can be mixed together)
- PubMedID to [Variomes][def4]: Retrieves genes and drugs from abstracts.
- PMCID or PubMedID to [SIBiLS][def5] (Swiss Institute of Bioinformatics Literature Services)

2. Select output type: DataFrame or BioCjson.  
3. Download results

<img
  src="https://github.com/gititub/test/blob/main/rsc/app.png"
  alt="Alt text"
  style="display: block; width:400px">


### Run the APP in Linux
```
cd app;shiny run --reload
```
You can also run NER-App in Windows.  


## Biomedical Entity Normalization

[Biomedical Entity Normalization](https://nerversetoolkit.shinyapps.io/normamedtoolbox1/)

<img
  src="https://github.com/gititub/test/blob/main/rsc/NEN.png"
  alt="Alt text"
  style="display: block; width: 400px">

1.	Make a query  
- [LitVar][def6] Normalization: e.g. BRAFp.V600E  (one or more, comma separated) or upload CSV/TSV/TXT  file with two mandatory columns,'gene' and 'HGVS. → dbSNP rs ID  s
- [SynVar][def7] Normalization : e.g. 19915144, MEK1(p.Q56P) or upload CSV file with three mandatory columns: 'pmid', gene' and 'HGVS'.
- Gene ncbi Normalization to gene ID (one by one, only a gene name or gene + specie)  
- Gene ID → Gene Name (one or more, comma separated)  
- Rs id → Gene Info (one or more, comma separated)
- PubMed ID to PMC ID (one or more, comma separated, or upload CSV/TSV/TXT file with tab separated IDs)*
- PMC ID to PubMed ID (one or more, comma separated, or upload CSV/TSV/TXT file with tab separated IDs)*
  
  (*) It is not necessary a column name and there can be more than one column, it just will read the first column.  
  
2.	Download results  


### Run the APP in Linux
```
cd appNEN;shiny run --reload
```


# Run Full Pipeline 

<img
  src="https://github.com/gititub/test/blob/main/rsc/workflow2.png"
  alt="Alt text"
  style="display: block; height:450px; width:800px">


## Installation

```
git clone https://github.com/gititub/NEN-NER-TOOLKIT.git
cd NEN-NER-TOOLKIT

```
Create a new Conda environment from environment.yml file.
```
conda env create -f myenv.yml
conda activate env_toolkit
pip install -r requirements.txt
```
Here is a step-by-step guide to installing Conda and creating your own Conda environment:

https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

  1. Download Anaconda: Head over to the Anaconda website (https://www.anaconda.com/download/) and download the latest version.

  2. In your terminal window, run: bash Anaconda-latest-Linux-x86_64.sh

## Run Full Pipeline from the Command Line on Linux:

$ ./ann.sh [input file or directory] [json/tsv]

This code first checks if the input is a directory and then process each file within the directory while applying the appropriate logic based on the file content, PubMed ID or PMC ID. Based on which function is being executed, the output files will have distinct names: "pmids", "PMC", "bern" (only PubMedIDs), "gwasminer" (only PMCIDs) and "ptc" will be appended to the output file name. Additionally, the script distinguishes between two output formats: 'biocjson' or 'dataframe' depending on whether the second argument is 'json' or 'tsv' respectively.

Results will be saved in a directory named "results_[datetime]".

Run example of a directory: 
```
chmod +x *.sh
./ann.sh example tsv
```
Or of a single file:
```
./ann.sh example/pmids_sample.tsv json
```
```
./ann.sh example/pmcs_sample.txt tsv
```

ℹ️ **Accepted file input extensions for all commands include .txt, .tsv, or .csv. The data must be organized into a single column of elements. If it's a number, interpret it as a PubMedID; if it starts with "PMC," interpret it as a PMC ID. The first row may also include a column name, which can be either 'pmid' or 'PMC'.**

If the output filename concludes with '.tsv', you will receive the results as a DataFrame. However, if it concludes with '.json', the results will be provided in the bioCjson format.


## NER, NEN and RE for full-text articles with [PubTator3][def2] from PubMedIDs or PMCIDs.

$ python src/ptc.py -i [file_path_pmids/pmcids] -o [output_filename]

Run example: 
```
python src/ptc.py -i example/pmids_sample.tsv -o output.json
```
```
python src/ptc.py -i example/pmcs_sample.txt -o output.tsv
```
In case you choose 'tsv', it returns two dataframes: `output.tsv` with entities, and `output_relations.tsv` with correlations between entities. 


## NER, NEN and Relations with [PubTator3][def2] from a query.

This command will search for PMC articles related to a query, for example *biotin* or *Hodgkin+Lymphoma* (it is not case sensitive), using Bio.Entrez (Spaces may be replaced by '+' sign). Retrieves and processes PubTator annotations and save the results in the specified output file in biocjson or tsv format.

The fourth command-line argument is the number of IDs to retrieve. The last command-line argument is the oldest publication date. If you don't want to filter by date, simply omit the `-date` argument. This last argument is not mandatory.


$ python src/ptc_extract_pmc_query.py -q [query] -o [output_filename] -max [max retrievals] -date ["YYYY/MM/DD"]

Run example: published after January 1, 2022

```
python src/ptc.py -q biotin -o biotin.tsv -max 20 -date "2022/01/01"
```

```
python src/ptc.py -q multiple+sclerosis -o MS.json -max 10

```

## NER and NEN for PubMed abstracts with [BERN2][def]


$ python src/bern_extract_pmids.py [file_path_pmids] [output_filename]

```
python src/bern_extract_pmids.py example/pmids_sample2.csv output_bern.tsv
```
```
python src/bern_extract_pmids.py example/pmids_sample.tsv output_bern.json
```

## NER and NEN for plain text with [BERN2][def]

Plain text is limited to 5000 characters. To speed up the process, texts (each no longer than 5000 characters) are distributed into subdirectories within the designated input directory (with no more than 120 files per subdirectory). This process is parallelized to concurrently process the subdirectories. The outcome is acquired in a file bearing the same name as the subdirectory, following the specified format (JSON or TSV) – one result file per subdirectory. 

You can download the full text automatically divided into subdirectories yourself, from a PMC list or from a folder with pdf files.

From PMC list:

$ python src/download_pmc_fulltext.py [file_path_pmcs] [output_dir]

Run example: 
```
python src/download_pmc_fulltext.py example/pmcs_sample.txt bern_ft
```
Or from a one or more PMCs:
```  
python src/download_pmc_fulltext.py PMC2907921,PMC8885421 bern_ft2
```

Additionally, you can use **Selenium** and **ChromeDriver**. ℹ️ The ChromeDriver’s primary function is to start Google Chrome. Without them, it is impossible to automate any website and run Selenium. To use ChromeDriver, you need to first download it from the Chromium website and then install it. https://chromedriver.chromium.org/downloads  

$ python src/download_fulltext_bern.py [file_path_pmcs] [output_dir]

```
python src/download_fulltext_bern.py example/pmcs_sample.txt bern_ft
```
From **PDF files**:

$ python src/convert_pdf.py [pdf_directory] [output_directory]

Run example: 
```
python src/convert_pdf.py example/pdf_files bern_ft_pdf
```
Then, you can run [BERN2][def]:

⚠️ This process might take several minutes. 


$ python src/bern_extract_ann.py [input_dir] [output_dir] [json/df]

Run example: 
```
python src/bern_extract_ann.py bern_ft_pdf bern_ft_pdf_results df
```


## Normalize Variants with [SynVar][def7] and [LitVar][def6]

To run example use `pmid_gene_HGVS.csv/tsv` as input file, or use your own data with 3 columns: pmid, gene, HGVS. Returns two files in the specified output directory, one with [LitVar][def6]normalization and the second one with [SynVar][def7] normalization and gene+drug NER.

$ python src/normalize.py [input_file] [output_directory] 

```
python src/normalize.py example/pmid_gene_HGVS.csv '.'
```

## ID converter

### Convert PubMed ids to PMCIDs

Here, PubMedIDs is your input file containing the list of PubMedIDs (tsv, csv and txt format allowed), `output_directory` is the directory where you want to save the output file, and _pmc is the suffix you want to add to the output file's name. This script will read the input PubMedIDs, retrieve PMCIDs, and save the output file in the specified directory with the desired name.

$ python src/pmc_from_pmid.py [pmids] [output_directory] [_pmc]
 
Run example: 
```
python src/pmc_from_pmid.py example/pmids_sample.tsv '.' _pmc
```
### Convert PMC ids to PubMedIDs

Run example: 
```
python src/pmid_from_pmc.py example/pmcs_sample.txt '.' _pmid
```


[def]: http://bern2.korea.ac.kr/
[def0]: https://github.com/Thomas-Rowlands/GWAS-Miner
[def1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8885717/
[def2]: https://www.ncbi.nlm.nih.gov/research/pubtator3/api
[def3]: https://fastdatascience.com/drug-named-entity-recognition-python-library/
[def4]: https://variomes.text-analytics.ch/apis
[def5]: https://sibils.text-analytics.ch/doc/api/search/
[def6]: https://www.ncbi.nlm.nih.gov/research/litvar2/
[def7]: https://synvar.text-analytics.ch/