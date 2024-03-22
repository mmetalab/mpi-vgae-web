from dash import dcc, html, Input, Output, State, callback, dash_table
import dash_bootstrap_components as dbc
import base64
import io
from pathlib import Path
import pandas as pd

from utils import dataManager as dm
from utils import layoutFunctions as lf
from utils import callbackFunctions as cf

# ------------------------------------------------------------------------------
# Initialize utility objects and useful functions
# ------------------------------------------------------------------------------

# id function that helps manage component id names. It pre-pends
# the name of the page to a string so that writing ids specific for each page is easier 
id = cf.id_factory('interactions')          

# Full path of the data folder where to load raw data
dataFolder = Path(__file__).parent.parent.absolute() / 'data'

# Load the Atlas dataFrame with all structures, acronyms, colors etc
structuresDf = cf.loadStructuresDf(dataFolder/'structures.json')

# ------------------------------------------------------------------------------
# Load the necessary data
# ------------------------------------------------------------------------------

# Metrics data for WFA and PV
wfa = dm.readMetricsDataForGenes(dataFolder/'originalData/data_SD1.xlsx')
pv = dm.readMetricsDataForGenes(dataFolder/'originalData/data_SD2.xlsx')
# Load Genes data
geneDict = dm.readGenesCorrelationSupplData(dataFolder/'originalData/data_SD4.xlsx')
# genesDf = df = pd.read_excel(dataFolder/'originalData/data_SD4.xlsx', header=0, index_col=0)

# Load ISH data
ish_en =  pd.read_csv(dataFolder/'gene_expression_ABA_energy.csv', index_col=0)
ish_en.columns = pd.to_numeric(ish_en.columns)

# ------------------------------------------------------------------------------
# Perform some preprocessing
# ------------------------------------------------------------------------------

genome_dict = {'Homo sapiens':"Homo_sapiens",
'Mus musculus':"Mus_musculus",
'Rattus norvegicus':"Rattus_norvegicus",
'Escherichia coli':"Escherichia_coli",
'Bos taurus':"Bos_taurus",
'Pseudomonas aeruginosa':"Pseudomonas_aeruginosa",
'Arabidopsis thaliana':"Arabidopsis_thaliana",
'Saccharomyces cerevisiae':"Saccharomyces_cerevisiae",
'Drosophila melanogaster':"Drosophila_melanogaster",
'Caenorhabditis elegans':"Caenorhabditis_elegans"}

# ------------------------------------------------------------------------------
# LAYOUT
# ------------------------------------------------------------------------------


layout = dbc.Container([
    lf.make_CitationOffCanvas(id),
    lf.make_AboutUsOffCanvas(id),
    lf.make_GeneInfoModal(id),
    dbc.Row(lf.make_NavBar()),                           # Navigation Bar
    dbc.Row(lf.make_InteractionHeader(id)),            # Big header
    #
    dbc.Row([lf.make_Subtitle('Metabolite-Protein Featurization and Visualization')]),
    dbc.Row([
        dbc.Col([
            lf.make_InteractionUploadMenu(id),
        ],xs=12,lg=5),
        dbc.Col([
            lf.make_InteractionTable(id),
        ],xs=20,lg=5)
    ], className = 'align-items-center'),
    html.Br(),
    dbc.Row([lf.make_Subtitle('Metabolite-Protein Interaction Prediction')]),
    dbc.Row([
        dbc.Col(lf.make_MPISelectionMenu(id, genome_dict, geneDict['wfa_en']),
            xs=12,lg=4, className='mt-5'
        ),
        # dbc.Col(lf.make_MPIInteractionResult(id),
        #     xs=12,lg=4, className='mt-5'
        # ),
        dbc.Col(
            dbc.Spinner(
                dcc.Graph(
                    figure=cf.make_MPIScatter(),
                    id=id('mpiPlot'), config={'displaylogo':False}, className='mt-3'),
                color='primary'
            )
        )
    ]),

    dbc.Row([lf.make_CollapsableTable(id)]),
    dbc.Row([lf.make_CC_licenseBanner(id)]),
    dbc.Row([],style={"margin-top": "500px"}),
])


@callback(
    Output(component_id=id('store-metabolite-data-upload'),component_property='data'),
    Input(component_id=id('upload-button'),component_property='n_clicks'),
    Input(component_id=id('upload-metabolite-data'),component_property='contents'),
    Input(component_id=id('metabolite-prefill-textarea'),component_property='value')
)

def update_metabolite(n_clicks, list_of_contents, value):
    if n_clicks > 0:
        df = pd.DataFrame(columns = ['Metabolite Name', 'HMDB ID','SMILES'])
        if list_of_contents is not None:
            tt = [parse_contents_df(c) for c in list_of_contents]
            df = tt[0]
        else:
            for i in value.split('\n'):
                if i == '': continue
                i = i.split(',')
                row = pd.DataFrame({'Metabolite Name': i[0], 'HMDB ID': i[1], 'SMILES': i[2]}, index=[0])
                df = pd.concat([df, row]).reset_index(drop=True)
        return df.to_json(orient='split')

def parse_contents_df(contents):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in contents:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in contents:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    return df

@callback(
    Output(component_id=id('table-metabolite-data-upload'),component_property='children'),
    Input(component_id=id('store-metabolite-data-upload'),component_property='data'),
    Input(component_id=id('upload-button'),component_property='n_clicks')
)

def update_metabolite_data(metabolite_data, n_clicks):
    if n_clicks > 0:
        df = pd.read_json(metabolite_data, orient='split')
        tab = dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True, size='sm')
        return tab


@callback(
    Output(component_id=id('store-protein-data-upload'),component_property='data'),
    Input(component_id=id('upload-button'),component_property='n_clicks'),
    Input(component_id=id('upload-protein-data'),component_property='contents'),
    Input(component_id=id('protein-prefill-textarea'),component_property='value')
)

def update_protein(n_clicks, list_of_contents, value):
    if n_clicks > 0:
        df = pd.DataFrame(columns = ['UniprotID','Protein Name','Gene Name','Organism','Sequence'])
        if list_of_contents is not None:
            tt = [parse_contents_df(c) for c in list_of_contents]
            df = tt[0]
        else:
            for i in value.split('\n'):
                if i == '': continue
                i = i.split(',')
                row = pd.DataFrame({'UniprotID': i[0], 'Protein Name': i[1], 'Gene Name': i[2], 'Organism': i[3], 'Sequence': i[4]}, index=[0]) 
                df = pd.concat([df, row]).reset_index(drop=True)
        return df.to_json(orient='split')


@callback(
    Output(component_id=id('table-protein-data-upload'),component_property='children'),
    Input(component_id=id('store-protein-data-upload'),component_property='data'),
    Input(component_id=id('upload-button'),component_property='n_clicks')
)

def update_protein_data(protein_data, n_clicks):
    if n_clicks > 0:
        df = pd.read_json(protein_data, orient='split')
        tab = dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True, size='sm')
        return tab


@callback(
    Output(component_id=id('store-molecular-feature'),component_property='data'),
    Input(component_id=id('feature-button'),component_property='n_clicks'),
    Input(component_id=id('drpD_genome_Select'),component_property='value'),
    Input(component_id=id('store-protein-data-upload'),component_property='data'),
    Input(component_id=id('store-metabolite-data-upload'),component_property='data'),
)

def update_molecular_feature(n_clicks,genome,protein_data, metabolite_data):
    if n_clicks > 0:
        if protein_data is not None and metabolite_data is not None:
            df_protein = pd.read_json(protein_data, orient='split')
            df_metabolite = pd.read_json(metabolite_data, orient='split')
            umap_df,node_feats,feature_df,g,dg = cf.make_molecule_feature_scatter(df_protein, df_metabolite,genome)
            tt = cf.make_MPIPrediction(node_feats,feature_df,g,dg,df_protein)
            # tab = dbc.Table.from_dataframe(umap_df, striped=True, bordered=True, hover=True, size='sm')
        return umap_df.to_json(orient='split')



@callback(
    Output(component_id=id('mpiPlot'), component_property='figure'),
    # Output(component_id=id('collps_Tab'), component_property='children'),
    State(component_id=id('mpiPlot'), component_property='figure'),
    Input(component_id=id('store-molecular-feature'),component_property='data'),
    # Input(component_id=id('drpD_geneSelect'), component_property='value'),
    # Input(component_id=id('drpD_metricSelector'), component_property='value'),
)
def updateMPI(fig, data):
    data = pd.read_json(data, orient='split')
    print(data.head())
    # Update the table with Gene info
    # g, geneName = cf.getGeneInfoTable(selMetric, selGene, geneDict)
    # tab = dbc.Table.from_dataframe(data, striped=True, bordered=True, hover=True)
    # metricData = cf.getMetricDf(selMetric, wfa, pv)
    # aggreDf = cf.combineGenesDf(selGene, metricData, ish_en, structuresDf)
    fig = cf.update_MPIScatter(fig,data)

    return fig

'''

@callback(
    Output(component_id=id('corrPlot'), component_property='figure'),
    Output(component_id=id('collps_Tab'), component_property='children'),
    State(component_id=id('corrPlot'), component_property='figure'),
    Input(component_id=id('drpD_geneSelect'), component_property='value'),
    Input(component_id=id('drpD_metricSelector'), component_property='value'),
)
def updateGenecorr(fig, selGene, selMetric):
    # Update the table with Gene info
    g, geneName = cf.getGeneInfoTable(selMetric, selGene, geneDict)
    tab = dbc.Table.from_dataframe(g, striped=True, bordered=True, hover=True)

    metricData = cf.getMetricDf(selMetric, wfa, pv)
    aggreDf = cf.combineGenesDf(selGene, metricData, ish_en, structuresDf)
    fig = cf.update_GenesScatter(fig, aggreDf, structuresDf, geneName)

    return fig, tab


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
    Output(component_id=id('moreInfoCollapse'), component_property='is_open'),
    Input(component_id=id('moreInfoIcon'), component_property='n_clicks'),
    State(component_id=id('moreInfoCollapse'), component_property='is_open'),
    prevent_initial_call=True
)
def invertMoreInfoVisibility(n_clicks, is_open):
    if n_clicks:
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
    Output(component_id=id('collps_Tab'), component_property='is_open'),
    Output(component_id=id('btn_openTabDiffuse'), component_property='children'),
    Output(component_id=id('btn_openTabDiffuse'), component_property='color'),
    Input(component_id=id('btn_openTabDiffuse'),component_property='n_clicks'),
    State(component_id=id('collps_Tab'), component_property='is_open'),
    prevent_initial_call=True
)
def invertTabVisibility( _ , previousState):
    newState = not previousState
    if newState:
        text = 'Collapse Gene Info'
        color = 'info'
    else:
        text = 'Open Gene Info'
        color = 'primary'
    return newState, text, color



@callback(
    Output(component_id=id('modal_info'), component_property='is_open'),
    Input(component_id=id('btn_info'),component_property='n_clicks'),
    State(component_id=id('modal_info'), component_property='is_open'),
)
def invertModalInfoVisibility(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open
'''