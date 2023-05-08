######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
import pickle
import networkx as nx
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit import Chem
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import plotly.express as px
import dgl
import torch

# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem import DataStructs
# from collections import OrderedDict
# from sklearn import preprocessing
# from descriptastorus.descriptors import rdNormalizedDescriptors
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

image = Image.open('./images/feature-generation.png')

st.image(image, use_column_width=True)

st.subheader('Generate molecular features.')

#Load data
df_m = st.session_state['df_m']
df_p = st.session_state['df_p']

st.set_option('deprecation.showPyplotGlobalUse', False)

st.markdown('<div style="text-align: justify;"> Select an organism to compute the MPI.\
</div>', unsafe_allow_html=True)


mode_dict = {'Homo sapiens':"Homo_sapiens",
'Mus musculus':"Mus_musculus",
'Rattus norvegicus':"Rattus_norvegicus",
'Escherichia coli':"Escherichia_coli",
'Bos taurus':"Bos_taurus",
'Pseudomonas aeruginosa':"Pseudomonas_aeruginosa",
'Arabidopsis thaliana':"Arabidopsis_thaliana",
'Saccharomyces cerevisiae':"Saccharomyces_cerevisiae",
'Drosophila melanogaster':"Drosophila_melanogaster",
'Caenorhabditis elegans':"Caenorhabditis_elegans",
'Global':'All'}

mode = st.selectbox(
    'Due to computing limitations, select one organism at a time.',
    list(mode_dict.keys()))

st.write('The server will generate features for the following organism:',mode)

mpi_file_name = './MPI-network/'+'pca_mpi_'+str(mode_dict[mode])+'.pkl'
df_file_name = './Features/'+'pca_feature_df_'+str(mode_dict[mode])+'.pkl'
g = pickle.load(open(mpi_file_name, "rb" ))
node_feats = pickle.load(open(df_file_name, "rb" ))

dg = dgl.from_networkx(g)
features = np.stack(node_feats['pca_128'].values, axis=0)
features = torch.from_numpy(features)
dg.ndata['h'] = features
dg.ndata['h'].type()



dict_mets = dict(zip(df_m['Metabolite Name'].tolist(),df_m['SMILES'].tolist()))
dict_protein = dict(zip(df_p['UniprotID'].tolist(),df_p['Sequence'].tolist()))

def feats_convert(smile):
    mol = Chem.MolFromSmiles(smile)
    counts = mol.GetNumAtoms()
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    array = np.zeros((0, ), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, array)
    return array

mets_pca = pickle.load(open("./Features/mets_pca.pkl",'rb'))
protein_vec = pickle.load(open("./Features/protein_vector.p", "rb" ))
protein_pca = pickle.load(open("./Features/protein_pca.pkl",'rb'))

mp_button = st.button("Click here to run generate features") # Give button a variable name
if mp_button: # Make button a condition.
    st.text("Start generate molecular features")
    feature_df = pd.DataFrame(columns=['Name','Type','Feature'])
    index = 0
    for k,v in dict_mets.items():
        f = feats_convert(v)
        pca_f = mets_pca.transform(f.reshape(1, -1))
        row = [k,'Metabolite',pca_f]
        feature_df.loc[index] = row
        index +=1

    for k,v in dict_protein.items():
        if k in protein_vec:
            f = protein_vec[k]
            pca_f = protein_pca.transform(f.reshape(1, -1))
            row = [k,'Protein',pca_f]
            feature_df.loc[index] = row
            index +=1 
    
    feature_df = pd.merge(feature_df,df_m,left_on='Name', right_on='Metabolite Name',how='left')
    feature_df = pd.merge(feature_df,df_p,left_on='Name', right_on='UniprotID',how='left')
    
    cols = ['Name','Type','HMDB ID','Feature']
    feature_df = feature_df[cols]

    tsne = TSNE(n_components=2, random_state=0)
    Y = feature_df["Type"].tolist()
    Z = feature_df["Name"].tolist()
    X = np.stack(feature_df['Feature'].values, axis=0)
    X = X.reshape(-1,128)
    X_2d = tsne.fit_transform(X)
    tsne_df = pd.DataFrame({'Name':Z,
                            'tSNE-1':X_2d[:,0],
                            'tSNE-2':X_2d[:,1],
                            'Label':Y})
    fig = px.scatter(
        tsne_df,
        x="tSNE-1",
        y="tSNE-2",
        color="Label"
    )
    st.text("Finished")
    st.plotly_chart(fig, theme="streamlit", use_container_width=True)
    st.session_state['feats_df'] = feature_df
    st.session_state['df_m'] = df_m
    st.session_state['df_p'] = df_p
    st.session_state['mpi'] = g
    st.session_state['dg'] = dg
    st.session_state['node_feats'] = node_feats

    st.download_button(label='Download features',
            data= pickle.dumps(feature_df),
            file_name='feature_df.pkl')