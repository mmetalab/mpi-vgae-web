from dash import dcc, html, Input, Output, State, callback
import dash_bootstrap_components as dbc

from pathlib import Path
import pandas as pd

from utils import dataManager as dm
from utils import layoutFunctions as lf
from utils import callbackFunctions as cf


# ------------------------------------------------------------------------------
# Initialize utility objects and useful functions
# ------------------------------------------------------------------------------

# id function that helps manage component id names. It pre-pends
# the name of the page to a string so that writing ids specific for each page is trtrSeasier 
id = cf.id_factory('database')          

# Full path of the data folder where to load raw data
dataFolder = Path(__file__).parent.parent.absolute() / 'data'

df = pd.read_csv(dataFolder/'mpidatabase/MPIDB_May2024.csv')

# ------------------------------------------------------------------------------
# LAYOUT
# ------------------------------------------------------------------------------
layout = dbc.Container([
    lf.make_CitationOffCanvas(id),
    lf.make_AboutUsOffCanvas(id),
    lf.make_MetricInfoModal(id),
    dbc.Row(lf.make_NavBar()),                  # Navigation Bar
    dbc.Row(lf.make_MPIDBHeader(id)),             # Big header
    dbc.Row([lf.make_Subtitle('Enzymatic reaction database Information')]),
    lf.make_DatabaseInfo(id,df),
    dbc.Row([lf.make_CC_licenseBanner(id)]),
    dbc.Row([],style={"margin-top": "500px"}),
])


# ------------------------------------------------------------------------------
# CALLBACKS
# ------------------------------------------------------------------------------

@callback(
    Output(component_id=id('offCanv_cite'), component_property='is_open'),
    Input(component_id=id('btn_citeHeader'),component_property='n_clicks'),
    Input(component_id='citeDropdown', component_property='n_clicks'),
    State(component_id=id('offCanv_cite'), component_property='is_open'),
    prevent_initial_call=True
)
def invertCiteMenuVisibility(n_clicks, n_clicks_dropdown, is_open):
    if n_clicks or n_clicks_dropdown:
        return not is_open
    return is_open


@callback(
    Output(component_id=id('offCanv_abtUs'), component_property='is_open'),
    Input(component_id='aboutUsDropdown', component_property='n_clicks'),
    State(component_id=id('offCanv_abtUs'), component_property='is_open'),
    prevent_initial_call=True
)
def invertAboutusMenuVisibility(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open


@callback(
    Output(component_id=id('modal_info'), component_property='is_open'),
    Input(component_id=id('btn_info_anat'),component_property='n_clicks'),
    Input(component_id=id('btn_info'),component_property='n_clicks'),
    State(component_id=id('modal_info'), component_property='is_open'),
)
def invertModalInfoVisibility(n_clicks_anat, n_clicks, is_open):
    if n_clicks or n_clicks_anat:
        return not is_open
    return is_open

@callback(
    Output(component_id=id('moreInfoCollapse'), component_property='is_open'),
    Input(component_id=id('moreInfoIcon'), component_property='n_clicks'),
    State(component_id=id('moreInfoCollapse'), component_property='is_open'),
    prevent_initial_call=True
)
def invertMoreInfoVisibility(n_clicks, is_open):
    if n_clicks:
            return not is_open
    return is_open
