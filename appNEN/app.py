# app.py
import json
import pandas as pd
import shinyswatch
import os
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shiny.types import FileInfo
from code import fetch_litvar_data, litvar_data, get_gene_info_by_gene_number, get_gene_info_by_gene_name, \
    get_gene_info_by_rsid


app_ui = ui.page_fluid(
    shinyswatch.theme.minty(),
    # shinyswatch.theme.darkly(),
    # shinyswatch.theme.sketchy(),
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
        ),
        ui.column(
            6,
            ui.h2("Biomedical Entity Normalization"),
        ),
    ),
    ui.br(),
    ui.input_select(
        "input_type",
        "Select function:",
        {
            "norm": "Variant Normalization",
            "gene": "Gene Normalization",
            "gene_name": "Gene Name",
            "gene_info": "rs ID Gene Info",
        },
        selected="Variant Normalization",
    ),
    ui.output_text_verbatim("txt"),
    ui.input_text_area("id", "Write query:",
                       placeholder='e.g.rs5030858 or BRAF',
                       width='800px', height='200px'),
    ui.row(
        ui.column(
            4,
            ui.input_file("file", "Choose CSV or TSV File",
                          accept=[".csv", ".tsv"],
                          multiple=False),
        ),
        ui.column(
            2,
            ui.br(),
            ui.input_action_button(
                "clear", "Clear", class_="btn-primary"
            ),
        ),
    ),
    ui.br(),
    ui.input_switch("all_results", "View all results", True),
    ui.input_action_button(
        "action", "Submit", class_="btn-primary"
    ),
    ui.download_button("download", "Download results"),
    ui.output_data_frame("table"),
)


def server(input, output, session):
    @output
    @render.image
    def imagen():
        return {
            "src": 'aprender.png',
            "style": "width: 80px; max-height: 80px;",
        }

    @output
    @render.text
    def txt():
        if input.input_type() == "gene":
            return "e.g. BRAF or braf homo sapiens"
        elif input.input_type() == 'gene_name':
            return "e.g. 7157,657,4234"
        elif input.input_type() == 'gene_info':
            return "e.g. rs5030858"
        else:
            return "e.g. BRAF p.V600E (or upload CSV file with two mandatory columns,'gene' and 'HGVS)"

    @reactive.Calc
    def result():
        if input.input_type() == "gene":
            result = get_gene_info_by_gene_name(input.id())
        elif input.input_type() == 'gene_name':
            result = get_gene_info_by_gene_number(input.id())
        elif input.input_type() == 'gene_info':
            result = get_gene_info_by_rsid(input.id())
        else:
            if input.file():
                f: list[FileInfo] = input.file()
                file = pd.read_csv(f[0]["datapath"], header=0, sep='\t')
                result = fetch_litvar_data(file)
            else:
                result = litvar_data(input.id())
        return result

    @reactive.Effect
    @reactive.event(input.clear)
    def _():
        if input.file():
            return input.file() == None

    @output
    @render.data_frame
    @reactive.event(input.action)
    def table():
        if input.all_results():
            return render.DataGrid(
                result(),
                width="100%",
                height="100%",
                filters=True,
            )
        else:
            return result().head(15)

    @session.download()
    def download():
        result().to_csv('results.tsv', sep='\t', index=False)
        path = os.path.join(os.path.dirname(__file__), "results.tsv")
        return path

    @reactive.Effect
    def _():
        req(input.action())
        ui.notification_show("Go!!", duration=5, type="message")


app = App(app_ui, server)
