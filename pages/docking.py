from dash import dcc, html, Input, Output, State, callback
import dash_bootstrap_components as dbc

from pathlib import Path
import pandas as pd

from utils import dataManager as dm
from utils import layoutFunctions as lf
from utils import callbackFunctions as cf

import py3Dmol

# ------------------------------------------------------------------------------
# Initialize utility objects and useful functions
# ------------------------------------------------------------------------------

# id function that helps manage component id names. It pre-pends
# the name of the page to a string so that writing ids specific for each page is easier 
id = cf.id_factory('docking')          

# Full path of the data folder where to load raw data
dataFolder = Path(__file__).parent.parent.absolute() / 'data'

# ------------------------------------------------------------------------------

from dash import Dash, html, dash_table, Input, Output, callback
import dash_bio as dashbio
from dash_bio.utils import PdbParser, create_mol3d_style
import pandas as pd

# Parse AutoDock Vina output file to extract docked poses
def parse_vina_output(output_file):
    # Parse the output file and extract the docked poses
    # Return a list of poses (molecules or PDB blocks)

    # Example code to parse output file:
    poses = []
    with open(output_file, 'r') as f:
        # Read lines from the output file and extract poses
        # Example: parse each pose and convert it into RDKit Mol object or PDB block
        for line in f:
            if line.startswith('MODEL'):
                pose = ''  # Start collecting the pose data
            elif line.startswith('ENDMDL'):
                poses.append(pose)  # Append the pose to the list
            else:
                pose += line  # Add line to the current pose data

    return poses

# vina_output_file = dataFolder/'mpidatabase/DB00014_AF-P01704-F1-model.pdbqt'  # Adjust file path
# poses = parse_vina_output(vina_output_file)

docking_db_file = dataFolder/'mpidatabase/Docking_DB.csv'
docking_db = pd.read_csv(docking_db_file)
protein_dict_uniprot_id = dict(zip(docking_db['Uniprot ID'].tolist(),docking_db['Uniprot ID'].tolist()))
protein_dict_pdb_id = dict(zip(docking_db['PDB'].tolist(),docking_db['Uniprot ID'].tolist()))
protein_dict_protein_name = dict(zip(docking_db['Protein Name'].tolist(),docking_db['Uniprot ID'].tolist()))
protein_dict = {**protein_dict_protein_name,**protein_dict_uniprot_id,**protein_dict_pdb_id}
uniprot_Id_dict = dict(zip(docking_db['Uniprot ID'].tolist(),docking_db['PDB'].tolist()))

metabolite_dict_hmdb_id = dict(zip(docking_db['HMDB ID'].tolist(),docking_db['HMDB ID'].tolist()))
metabolite_dict_metabolite_name = dict(zip(docking_db['Metabolite Name'].tolist(),docking_db['HMDB ID'].tolist()))
metabolite_dict = {**metabolite_dict_metabolite_name,**metabolite_dict_hmdb_id}


# ------------------------------------------------------------------------------
# LAYOUT
# ------------------------------------------------------------------------------
layout = dbc.Container([
    lf.make_CitationOffCanvas(id),
    lf.make_AboutUsOffCanvas(id),
    lf.make_GeneInfoModal(id),
    dbc.Row(lf.make_NavBar()),                  # Navigation Bar
    dbc.Row(lf.make_DockingHeader(id)),             # Big header

    dbc.Row([lf.make_Subtitle('Explore MPI result')]),
    dbc.Row([
        dbc.Col(lf.make_MPIVisualizationMenu(id, metabolite_dict, protein_dict),
            xs=12,lg=4, className='mt-5'
        ),
        dbc.Col(
            dbc.Spinner(
                cf.make_DockPlot(id)
                ),
        )
        
    ]),
    html.Br(),
    dbc.Row([lf.make_CC_licenseBanner(id)]),
    dbc.Row([],style={"margin-top": "500px"}),
])


@callback(
    Output(component_id=id('poses'), component_property='data'),
    Output(component_id=id('pdb'), component_property='data'),
    Input(component_id=id('drpD_metSelect'), component_property='value'),
    Input(component_id=id('drpD_ProteinSelect'), component_property='value'),
    Input(component_id=id('dock-button'), component_property='n_clicks')
)
def updatePDB(selMet, selProtein, n_clicks):
    if n_clicks > 0:
        # Update the table with Gene info
        # Load AutoDock Vina output file and parse docked poses
        uniprotid = protein_dict[selProtein]
        pdbid = uniprot_Id_dict[uniprotid]
        hmdbid = metabolite_dict[selMet]
        print(hmdbid,uniprotid,pdbid)
        pdb1 = './data/mpidatabase/DockingDB/PDB/%s.pdb' % pdbid
        vina_output_file = './data/mpidatabase/DockingDB/Docking/%s_%s_%s.pdbqt' % (hmdbid,uniprotid,pdbid)  # Adjust file path
        print(vina_output_file)
        poses = parse_vina_output(vina_output_file)
        return poses, pdb1

# Callback to update 3D molecular viewer
@callback(
    Output(component_id=id('mol-view'), component_property='children'),
    # State(component_id=id('mol-view'), component_property='children'),
    Input(component_id=id('poses'), component_property='data'),
    # Input(component_id=id('pose-selector'), component_property='value'),
    Input(component_id=id('pdb'), component_property='data')
)
def update_viewer(poses,pdb):
    view = py3Dmol.view(width=800, height=600)
    view.addModel(open(pdb, 'r').read(),'pdb')
    view.setStyle({'cartoon': {'color':'spectrum'}})
    view.addModel(poses[0], 'pdb')
    view.setStyle({'model':1},{'stick':{'colorscheme':'greenCarbon'}})
    view.zoomTo()
    # view.show()
    return html.Iframe(srcDoc=view._make_html(), width='100%', height='600')
    # return html.Div([html.H3(f'PDB ID: {selected_pose}'), html.Div(view.create_model(), style={'display': 'inline-block'})])

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