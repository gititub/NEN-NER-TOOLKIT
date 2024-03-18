# app.py
import json
import os
import pandas as pd
import shinyswatch
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shiny.types import ImgData
from code import count_characters, extract_pubtator, query_plain, \
    query, plain_drugs, download_from_PMC, download_from_PubMed, \
    bern_extract_pmids, synvar_ann, download_data, apply_db_from_wikipedia

app_ui = ui.page_fluid(
    shinyswatch.theme.superhero(),
    # ui.include_css("style.css"),
    ui.br(),
    ui.row(
        ui.column(
            1,
            ui.br(),
        ),
        ui.column(
            2,
            ui.output_image("imagen", height="90px"),
            ui.br(),
            ui.br(),
        ),
        ui.column(
            6,
            ui.h1("Biomedical Entity Recognition, Normalization and Relation Extraction"),
        ),
        ui.column(
            2,
            ui.output_image("imagen2", height="90px"),
            ui.tags.a("About & FAQ",
                      href="https://github.com/gititub/NER-NEN-TOOLKIT"),
        ),
    ),
    ui.layout_sidebar(
        ui.panel_sidebar(
            ui.input_select(
                "input_type",
                "Select:",
                {
                    "PTC": "PMCID or PubMedID (PubTator)",
                    "query": "Word in PubMed Central (PubTator)",
                    "relations": "Relations from PMCID or PubMedID (PubTator)",
                    "relations_query": "Relations from Word in PubMed Central",
                    "plain_text": "Plain Text (BERN2)",
                    "pmid_bern": "PubMedID (BERN2)",
                    "plain_drugs": "Plain Text (Drug NER)",
                    "id_drugs": "PMCID or PubMedID (Drug NER)",
                    "pmid_synvar": "PubMedID (Variomes)",
                    "download_text": "See Article by ID",
                },
                selected='Plain Text',
            ),
            ui.input_select(
                "output_type",
                "Output Type",
                {
                    "df": "Dataframe",
                    "biocjson": "BioCjson"
                },
                selected="BioCjson",
            ),
            ui.h6('Set Up Only For Query: '),
            ui.input_date('x', 'Date From:', value='2022-01-01'),
            ui.input_numeric('retmax', 'Max. Retrievals:', value=25),
            ui.br(),
            ui.input_switch("all_results", "Show All", True),
            class_="mb-3",
        ),
        ui.panel_main(
            ui.div(

                ui.div(
                    ui.output_text_verbatim("txt2"),
                ),
                ui.div(
                    {"class": "card mb-5"},
                    ui.div(
                        {"class": "card-body"},
                        ui.input_text_area("id", "Input Text:",
                                           placeholder="PMC or pmid (one or more, comma separated), word or plain text (less than 5000 characters for BERN2)",
                                           width='1200px', height='200px'),
                    ),
                    ui.div(
                        {"class": "card-footer"},
                        ui.output_text('ch'),
                    ),
                ),
                ui.div(
                    ui.div(

                        ui.input_action_button("action", "Submit", class_="btn-primary"),
                        ui.download_button("download", "Download results"),
                    ),
                ),
                ui.div(

                    ui.output_text_verbatim("txt"),
                    ui.output_data_frame("table"),

                ),
            ),
        ),
    ),
)

def server(input, output, session):
    @output
    @render.image
    def imagen():
        return {
            "src": 'literatura.png',
            "style": "width: 100px; max-height: 100px;",
        }

    @output
    @render.image
    def imagen2():
        return {
            "src": '5907623.png',
            "style": "width: 80px; max-height: 80px;",
        }

    @output
    @render.text
    def txt2():
        if input.input_type() == 'query':
            return "e.g. Hardy&Weinberg"
        elif input.input_type() == 'plain_text':
            return "e.g. Progressive multifocal leukoencephalopathy (PML) is a rare demyelinating disorder of the brain caused by reactivation of the JC virus (JCV), a polyomavirus that infects at least 60% of the population but is asymptomatic or results in benign symptoms in most people. "
        elif input.input_type() == 'pmid_bern':
            return "e.g. 27432226, 22383897"
        elif input.input_type() == 'plain_drugs':
            return "e.g. Ruxolitinib also is approved for treatment of patients with polycythemia vera who have had an inadequate response to or are intolerant of hydroxyurea."
        elif input.input_type() == 'relations_query':
            return "e.g. multiple&sclerosis"
        else:
            return "e.g. 27432226, 22383897 or 27432226, PMC5010513, PMC2907921"

    @output
    @render.text
    def ch():
        return f"Character counter: {count_characters(input.id())}"

    @reactive.Calc
    def result():
        if input.input_type() == "PTC":
            result = extract_pubtator(input.id(), input.output_type())[0]
        elif input.input_type() == "relations":
            result = extract_pubtator(input.id(), input.output_type())[1]
        elif input.input_type() == 'plain_text':
            result = query_plain(input.id(), input.output_type())
        elif input.input_type() == 'pmid_bern':
            result = bern_extract_pmids(input.id(), input.output_type())
        elif input.input_type() == 'query':
            ids = query(input.id(), input.retmax(), input.x())
            result = extract_pubtator(ids, input.output_type())[0]
        elif input.input_type() == 'relations_query':
            ids = query(input.id(), input.retmax(), input.x())
            result = extract_pubtator(ids, input.output_type())[1]
        elif input.input_type() == 'plain_drugs':
            result = plain_drugs(input.id(), input.output_type())
        elif input.input_type() == 'pmid_synvar':
            result = synvar_ann(input.id(), input.output_type())
        elif input.input_type() == "download_text":
            result = extract_pubtator(input.id(), input.output_type())[2]
        else:
            input_text = extract_pubtator(input.id(), 'df')[2]
            result = plain_drugs(input_text, input.output_type())
        return result

    @output
    @render.data_frame
    @reactive.event(input.action)
    def table():
        height = 450 if input.all_results() else None
        if isinstance(result(), pd.DataFrame):
            if input.all_results():
                return render.DataGrid(
                    result(),
                    width="100%",
                    height=height,
                    filters=True,
                )
            else:
                return render.DataGrid(
                    result().head(12),
                    width="100%",
                    height="100%",
                    filters=True,
                )

    @output
    @render.text
    @reactive.event(input.action)
    def txt():
        if input.output_type() == 'biocjson':
            if result():
                return result()
            else:
                return f"No results found. Try again."
        else:
            if not isinstance(result(), pd.DataFrame):
                return f"No results found. Try again."


    @render.download()
    def download():
        if input.output_type() == 'df':
            result().to_csv('results.tsv', sep='\t', index=False)
            path = os.path.join(os.path.dirname(__file__), "results.tsv")
            return path

        elif input.output_type() == 'biocjson':
            with open('results.json', 'w') as json_file:
                json.dump(result(), json_file, indent=4)
                path = os.path.join(os.path.dirname(__file__), "results.json")
                return path

    @reactive.Effect
    def _():
        req(input.action())
        ui.notification_show("Retrieving Entities!", duration=5, type="message")

app = App(app_ui, server)
