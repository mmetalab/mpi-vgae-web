######################
# Import libraries
######################
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
from st_aggrid import AgGrid, GridOptionsBuilder
import extra_streamlit_components as stx


# from rdkit import Chem
# from rdkit.Chem import Draw
# 隐藏streamlit默认格式信息
# hide_st_style = """
#             <style>
#             # MainMenu {visibility: hidden;}
#             footer {visibility: hidden;}
#             # header {visibility: hidden;}
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

image = Image.open('./images/load-data.png')
st.image(image, use_column_width=True)


######################
# Side Panel
######################

# st.sidebar.header('User scRNA-seq data (please upload files)')
uploaded_file = st.sidebar.file_uploader("Choose a CSV/TSV file with **metabolite** information", type="csv")
use_example_file = st.sidebar.checkbox(
    "Use example metabolite file", False, key='1', help="Use in-built example file to demo the app"
)

# If CSV is not uploaded and checkbox is filled, use values from the example file
# and pass them down to the next if block

st.header('Read metabolite data')

#Load data
if uploaded_file is not None:
    st.write("Using uploaded file")
    df_m = pd.read_csv(uploaded_file,sep=',')

elif use_example_file:
    st.write("An example file was loaded")
    file = './Data/example_metabolite_data.csv'
    df_m = pd.read_csv(file,sep=',')
 
else:
    st.write("An example file was preloaded")
    file = './Data/example_metabolite_data.csv'
    df_m = pd.read_csv(file,sep=',')

# Download example file
with open('./Data/example_metabolite_data.csv', 'rb') as f:
    s = f.read()
st.sidebar.download_button(
    label="Download example metabolite file",
    data=s,
    file_name='example_metabolite_data.csv',
    # mime='tsv',
)

st.header('Characteristics of loaded metabolite data')

gd = GridOptionsBuilder.from_dataframe(df_m)
gd.configure_side_bar(filters_panel=True)
gd.configure_pagination(enabled=True,paginationAutoPageSize=False,paginationPageSize=5)
gd.configure_default_column(groupable=True)
AgGrid(
    df_m,
    gridOptions=gd.build(),
)

mol_count = df_m.shape[0]
col1, col2 = st.columns(2)
col1.metric("Number of metabolites", mol_count)

uploaded_file = st.sidebar.file_uploader("Choose a CSV/TSV file with **protein** information", type="csv")
use_example_file = st.sidebar.checkbox(
    "Use example protein file", False, key='2', help="Use in-built example file to demo the app"
)

# If CSV is not uploaded and checkbox is filled, use values from the example file
# and pass them down to the next if block

st.header('Read protein data')

#Load data
if uploaded_file is not None:
    st.write("Using uploaded file")
    df_p = pd.read_csv(uploaded_file,sep=',')

elif use_example_file:
    st.write("An example file was loaded")
    file = './Data/example_protein_data.csv'
    df_p = pd.read_csv(file,sep=',')
 
else:
    st.write("An example file was preloaded")
    file = './Data/example_protein_data.csv'
    df_p = pd.read_csv(file,sep=',')

# Download example file
with open('./Data/example_protein_data.csv', 'rb') as f:
    s = f.read()
st.sidebar.download_button(
    label="Download example protein file",
    data=s,
    file_name='example_protein_data.csv',
    # mime='tsv',
)


# st.write('**Select an example molecule to show the structure**')
# example_mols = dict(zip(list(df_m.Name)[:10],list(df_m.SMI)[:10]))
# option = st.selectbox(
#     '',
#     list(example_mols.keys()))

# st.write('You selected:', option)
# if option:
#     compound_smiles = example_mols[option]
#     m = Chem.MolFromSmiles(compound_smiles)
#     Draw.MolToFile(m,'mol.png')
#     st.image('mol.png')


st.header('Characteristics of loaded protein data')

gd = GridOptionsBuilder.from_dataframe(df_p)
gd.configure_side_bar(filters_panel=True)
gd.configure_pagination(enabled=True,paginationAutoPageSize=False,paginationPageSize=5)
gd.configure_default_column(groupable=True)

AgGrid(
    df_p,
    gridOptions=gd.build(),
)

mol_count = df_p.shape[0]
col1, col2 = st.columns(2)
col1.metric("Number of proteins", mol_count)

st.session_state['df_m'] = df_m
st.session_state['df_p'] = df_p