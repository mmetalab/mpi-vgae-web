from dash import Dash, dcc, html, Input, Output
import dash_bootstrap_components as dbc
from pages import introduction, genes, database, blankPage, interactions

import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import pandas as pd
import os
import pickle
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit import Chem
import dgl
import torch
import umap
import networkx as nx
import torch.nn as nn
import torch.nn.functional as F

def cosine_similarity(arr1, arr2):
    dot_product = np.dot(arr1, arr2)
    norm_arr1 = np.linalg.norm(arr1)
    norm_arr2 = np.linalg.norm(arr2)
    similarity = dot_product / (norm_arr1 * norm_arr2)
    return similarity

def protein_mapper(proteins_df,node_feats):
    result_dict = dict(zip(proteins_df['Name'].tolist(),proteins_df['Feature'].tolist()))
    node_dict = dict(zip(list(node_feats.dbid),list(node_feats.pca_128)))
    mapp_dict = {}
    for k,v in result_dict.items():
        cos = 0
        for p,q in node_dict.items():
            t = cosine_similarity(v, q)
            if cos > t:
                cos = t
                mapp_dict[k] = [p,cos]
    return mapp_dict

def generate_edge_pair(metabolites,proteins,node_feats,adj_true):
    node_dict = dict(zip(list(node_feats.dbid),list(node_feats.index)))
    mets_index = []
    for i in metabolites:
        if i in node_dict:
         # print(i)
         mets_index.append(node_dict[i])
    protein_index = []
    for i in proteins:
        if i in node_dict:
         # print(i)
         protein_index.append(node_dict[i])
    test_pair = []
    true_pair = []
    for i in mets_index:
        for j in protein_index:
         test_pair.append([i,j])
         true_pair.append(adj_true[i,j])
    return node_dict,np.array(test_pair),true_pair

class MLPPredictor(nn.Module):
    def __init__(self, h_feats):
        super().__init__()
        self.W1 = nn.Linear(h_feats * 2, h_feats)
        self.W2 = nn.Linear(h_feats, 1)

    def apply_edges(self, edges):
        """
        Computes a scalar score for each edge of the given graph.
        Parameters
        ----------
        edges :
            Has three members ``src``, ``dst`` and ``data``, each of
            which is a dictionary representing the features of the
            source nodes, the destination nodes, and the edges
            themselves.
        Returns
        -------
        dict
            A dictionary of new edge features.
        """
        h = torch.cat([edges.src['h'], edges.dst['h']], 1)
        return {'score': self.W2(F.relu(self.W1(h))).squeeze(1)}

    def forward(self, g, h):
        with g.local_scope():
            g.ndata['h'] = h
            g.apply_edges(self.apply_edges)
            return g.edata['score']

def generate_result(pred_result,true,test_pair,node_feats):
    result = pd.DataFrame(columns=['Metabolite','Protein','Prediction Score','Existing'])
    for i in range(len(test_pair)):
        m_ind, p_ind = test_pair[i][0],test_pair[i][1]
        mets = node_feats.loc[m_ind]['node']
        prt = node_feats.loc[p_ind]['node']
        score = '%.5f' % torch.sigmoid(torch.tensor(pred_result[i])).numpy()
        if true[i] == 1:
            temp = 'Yes'
        else:
            temp = 'No'
        row = [mets,prt,score,temp]
        result.loc[i] = row
    return result

app = Dash(__name__,
    title="MPI-VGAE Web",
    external_stylesheets=[dbc.icons.FONT_AWESOME]
)

indexLayout = html.Div([
    dcc.Location(id='url', pathname='/introduction', refresh=False),
    html.Div(id='page-content')
])

# Create a "complete" layout for validating all callbacks. Otherwise when dash tries
# to validate them, most of them will thorw an error since they are linked to
# components that are not currently on the displayed page and so are not part of the 
# current layout
app.validation_layout = html.Div([
    indexLayout,
    introduction.layout,
    database.layout,
    interactions.layout,
    genes.layout,
    blankPage.layout
])

# This is the actual layout of the app
app.layout = indexLayout

@app.callback(
    Output('page-content', 'children'),
    Input('url', 'pathname')
)
def display_page(pathname):
    if pathname == '/introduction':
        return introduction.layout
    elif pathname == '/database':
        return database.layout
    elif pathname == '/interactions':
        return interactions.layout
    elif pathname == '/genes':
        return genes.layout
    else:
        return blankPage.layout


# This server object will be loaded by the WSGI script to be served as a webapp
# in a production server
server = app.server


# This will only be executed during debug when run locally, since WSGI does not 
# run this as __main__ but only takes the "server" variable
if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8080)
    # app.run_server(debug=True)

# gunicorn mpivgae:server