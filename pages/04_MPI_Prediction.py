######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
import pickle
from PIL import Image
import dgl
import torch
import torch.nn as nn
import torch.nn.functional as F
from st_aggrid import AgGrid, GridOptionsBuilder
import extra_streamlit_components as stx
import networkx as nx
# hide_st_style = """
#             <style>
#             #MainMenu {visibility: hidden;}
#             footer {visibility: hidden;}
#             header {visibility: hidden;}
#             <![]()yle>
#             """

# st.markdown(hide_st_style, unsafe_allow_html=True)

st.write("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Arimo');
html, body, [class*="css"]  {
   font-family: 'Arimo';
}
</style>
""", unsafe_allow_html=True)

######################
# Page Title
######################

image = Image.open('./images/mpi-prediction.png')
st.image(image, use_column_width=True)

#Load data

df_m = st.session_state['df_m'] 
df_p = st.session_state['df_p']
g = st.session_state['mpi']
node_feats = st.session_state['node_feats']
dg = st.session_state['dg']

try:
   feature_df = st.session_state['feats_df']
   st.write('Molecular features are loaded.')
except:
   st.write('Please note, molecular features are not generated yet.')

def convert_df(df):
   return df.to_csv(index=False).encode('utf-8')

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

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

model = torch.load('./Models/all_mpi_model.pth')
pred = torch.load('./Models/all_mpi_model_pred.pth')
adj_true = nx.adjacency_matrix(g,nodelist=node_feats['node'].tolist())
# test_pos_g = dgl.graph((test_pos_u, test_pos_v), num_nodes=dg.number_of_nodes(),device=device)

mpi_button = st.button("Click here to predict MPI") # Give button a variable name
if mpi_button: # Make button a condition.
   st.text("Start predicting MPI")

   metabolites = feature_df[feature_df['Type']=='Metabolite']['HMDB ID'].tolist()
   proteins = feature_df[feature_df['Type']=='Protein']['Name'].tolist()
   print('number of metabolites and proteins')
   print(len(metabolites),len(proteins))

   # metabolites = node_feats['dbid'].tolist()[:10]
   # proteins = node_feats['dbid'].tolist()[-10:]
   
   node_dict,test_pair,true_pair = generate_edge_pair(metabolites,proteins,node_feats,adj_true)
   chk = len(test_pair.shape)
   if chk == 1: 
      st.text("No available edges found in the MPI network")
   else:
      real_apply_g = dgl.graph((test_pair[:,0], test_pair[:,1]), num_nodes=dg.number_of_nodes(), device=device)
      print(test_pair)
      with torch.no_grad():
         pos_score = pred(real_apply_g, model)

      test_result = generate_result(pos_score,true_pair,test_pair,node_feats)
      print(test_result.head())

      st.text("Finished")
      st.write('Prediction result')

      gd = GridOptionsBuilder.from_dataframe(test_result)
      gd.configure_side_bar(filters_panel=True)
      gd.configure_pagination(enabled=True,paginationAutoPageSize=False,paginationPageSize=5)
      gd.configure_default_column(groupable=True)
      AgGrid(
         test_result,
         gridOptions=gd.build(),
      )

      # st.dataframe(test_result.head())
      csv_result = convert_df(test_result)
      st.download_button(
         "Press to download result",
         csv_result,
         "ccs_prediction.csv",
         "text/csv",
         key='download-csv'
      )