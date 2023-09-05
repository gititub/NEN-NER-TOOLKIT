# app.py
import json
import pandas as pd
import shinyswatch
import os
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from shiny.types import FileInfo
from code import fetch_litvar_data, litvar_data, get_gene_info_by_gene_number, get_gene_info_by_gene_name, \
    get_gene_info_by_rsid


def ui_card(title, *args):
    return (
        ui.div(
            {"class": "card mb-4"},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body"}, *args),
        ),
    )

app_ui = ui.page_fluid(
    shinyswatch.theme.minty(),
    # shinyswatch.theme.darkly(),
    # shinyswatch.theme.sketchy(),
    # ui.include_css("style.css"),
    ui.h2("Biomedical Entity Normalization"),
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
    ui.input_file("file", "Choose CSV File", accept=[".csv"], multiple=False),
    ui.input_switch("all_results", "View all results", True),
    ui.input_switch("filters", "Filter complete results", True),
    ui.input_action_button(
        "action", "Submit", class_="btn-primary"
    ),
    ui.download_button("download", "Download results"),
    ui.output_data_frame("table"),
)

def server(input, output, session):
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
            return "e.g. BRAFp.V600E (or upload CSV file with two mandatory columns,'gene' and 'HGVS)"

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

    @output
    @render.data_frame
    @reactive.event(input.action)
    def table():
        if input.all_results():
            return render.DataGrid(
                result(),
                width="100%",
                height="100%",
                filters=input.filters()  # True,
            )
        else:
            return result.head(15)

    @session.download()
    def download():
        result().to_csv('results.tsv', sep='\t', index=False)
        path = os.path.join(os.path.dirname(__file__), "results.tsv")
        return path

app = App(app_ui, server)
