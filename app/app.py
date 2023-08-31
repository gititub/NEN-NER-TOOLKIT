# app.py
import json
import shinyswatch
import os
from shiny import App, Inputs, Outputs, Session, reactive, render, req, ui
from code import extract_pubtator, extract_pubtator_from_pmcs, query_plain, extract_pubtator_from_pmcs_query


def ui_card(title, *args):
    return (
        ui.div(
            {"class": "card mb-4"},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body"}, *args),
        ),
    )


app_ui = ui.page_fluid(
    shinyswatch.theme.superhero(),
    # ui.include_css("style.css"),
    ui.h3("Biomedical Entity Recognition"),
    ui.input_select(
        "input_type",
        "id type",
        {
            "PMC": "PMC",
            "pmid": "pmid",
            "plain_text": "plain text",
            "query": "query",
        },
        selected='PMC',
    ),
    ui.input_text_area("id", "PMC, pmid or plain text",
                       placeholder='eg.PMC2882923 or 34216518 for pmids or biotin for query',
                       width='800px', height='200px'),
    ui.input_date('x', 'Date from'),
    ui.input_numeric('retmax', 'Max. retrievals', value=25),
    ui.input_select(
        "output_type",
        "Output type",
        {
            "df": "Dataframe",
            "biocjson": "BioCjson"
        },
        selected="BioCjson",
    ),
    ui.input_switch("all_results", "View all results", True),
    ui.input_switch("filters", "Filter complete results", True),
    ui.input_action_button(
        "action", "Submit", class_="btn-primary"
    ),
    ui.download_button("download", "Download results"),
    ui.output_text_verbatim("txt"),
    ui.output_data_frame("table"),
)


def server(input, output, session):
    @reactive.Effect
    @output
    @render.text
    @reactive.event(input.action)
    def txt():
        if input.output_type() == 'biocjson':
            if input.input_type() == "PMC":
                return extract_pubtator_from_pmcs(input.id(), input.output_type())
            elif input.input_type() == 'plain_text':
                return query_plain(input.id(), input.output_type())
            elif input.input_type() == 'pmid':
                return extract_pubtator(input.id(), input.output_type())
            else:
                return extract_pubtator_from_pmcs_query(input.id(), input.x(),
                                                        input.retmax(),
                                                        input.output_type())

    @output
    @render.data_frame
    @reactive.event(input.action)
    def table():
        if input.output_type() == 'df':
            if input.input_type() == "PMC":
                result = extract_pubtator_from_pmcs(input.id(), input.output_type())
                if input.all_results():
                    return render.DataGrid(
                        result,
                        width="100%",
                        height="100%",
                        filters=input.filters()  # True,
                    )
                else:
                    return result.head(15)
            elif input.input_type() == 'plain_text':
                result = query_plain(input.id(), input.output_type())
                if input.all_results():
                    return render.DataGrid(
                        result,
                        width="100%",
                        height="100%",
                        filters=input.filters()  # True,
                    )
                else:
                    return result.head(15)
            elif input.input_type() == 'pmid':
                result = extract_pubtator(input.id(), input.output_type())
                if input.all_results():
                    return render.DataGrid(
                        result,
                        width="100%",
                        height="100%",
                        filters=input.filters()  # True,
                    )
                else:
                    return result.head(15)
            else:
                result = extract_pubtator_from_pmcs_query(input.id(), input.x(),
                                                        input.retmax(),
                                                        input.output_type())
                if input.all_results():
                    return render.DataGrid(
                        result,
                        width="100%",
                        height="100%",
                        filters=input.filters()  # True,
                    )
                else:
                    return result.head(15)


    @session.download()
    def download():
        if input.output_type() == 'df':
            if input.input_type() == "PMC":
                result = extract_pubtator_from_pmcs(input.id(), input.output_type())
                result.to_csv('results.tsv', sep='\t', index=False)
                path = os.path.join(os.path.dirname(__file__), "results.tsv")
                return path
            elif input.input_type() == 'plain_text':
                result = query_plain(input.id(), input.output_type())
                result.to_csv('results.tsv', sep='\t', index=False)
                path = os.path.join(os.path.dirname(__file__), "results.tsv")
                return path
            elif input.input_type() == 'pmid':
                result = extract_pubtator(input.id(), input.output_type())
                result.to_csv('results.tsv', sep='\t', index=False)
                path = os.path.join(os.path.dirname(__file__), "results.tsv")
                return path
            else:
                result = extract_pubtator_from_pmcs_query(input.id(), input.x(),
                                                          input.retmax(),
                                                          input.output_type())
                result.to_csv('results.tsv', sep='\t', index=False)
                path = os.path.join(os.path.dirname(__file__), "results.tsv")
                return path

        elif input.output_type() == 'biocjson':
            if input.input_type() == "PMC":
                result = extract_pubtator_from_pmcs(input.id(), input.output_type())
                with open('results.json', 'w') as json_file:
                    json.dump(result, json_file, indent=4)
                path = os.path.join(os.path.dirname(__file__), "results.json")
                return path
            elif input.input_type() == 'plain_text':
                result = query_plain(input.id(), input.output_type())
                with open('results.json', 'w') as json_file:
                    json.dump(result, json_file, indent=4)
                path = os.path.join(os.path.dirname(__file__), "results.json")
                return path
            elif input.input_type() == 'pmid':
                result = extract_pubtator(input.id(), input.output_type())
                with open('results.json', 'w') as json_file:
                    json.dump(result, json_file, indent=4)
                path = os.path.join(os.path.dirname(__file__), "results.json")
                return path
            else:
                result = extract_pubtator_from_pmcs_query(input.id(), input.x(),
                                                          input.retmax(),
                                                          input.output_type())
                with open('results.json', 'w') as json_file:
                    json.dump(result, json_file, indent=4)
                path = os.path.join(os.path.dirname(__file__), "results.json")
                return path


app = App(app_ui, server)
