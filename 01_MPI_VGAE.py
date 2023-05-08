######################
# Import libraries
######################

import streamlit as st
from PIL import Image
hide_st_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            <![]()yle>
            """

st.markdown(hide_st_style, unsafe_allow_html=True)

######################
# Page Title
######################
st.write("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Arimo');
html, body, [class*="css"]  {
   font-family: 'Arimo';
}
</style>
""", unsafe_allow_html=True)

image = Image.open('./images/title-page.png')
st.image(image, use_column_width=True)

image = Image.open('./images/Figure1.tiff')
st.image(image, use_column_width=True)
st.markdown('<div style="text-align: justify;">Enzymatic reactions are crucial to explore the mechanistic\
function of metabolites and proteins in cellular processes and to understand the etiology of diseases.\
The increasing number of interconnected metabolic reactions allows the development of in silico deep\
learning-based methods to discover new enzymatic reaction links between metabolites and proteins to\
further expand the landscape of existing metabolite-protein interactome. Computational approaches to\
predict the enzymatic reaction link by metabolite-protein interaction (MPI) prediction are still very limited.\
In this study, we developed a Variational Graph Autoencoders (VGAE) based framework to predict MPI in genome-scale\
heterogeneous enzymatic reaction networks across ten organisms. By incorporating molecular features of metabolites\
and proteins as well as neighboring information in the MPI networks, our MPI-VAGE predictor achieved the best predictive\
performance compared to other machine learning methods. Moreover, when applying the MPI-VGAE framework to reconstruct hundreds\
of metabolic pathways, functional enzymatic reaction networks, and a metabolite-metabolite interaction network, our method showed\
the most robust performance among all scenarios. To the best of our knowledge, this is the first MPI predictor by VGAE for enzymatic\
reaction link prediction. Furthermore, we implemented the MPI-VGAE framework to reconstruct the disease-specific MPI network based on the\
disrupted metabolites and proteins in Alzheimer’s disease and colorectal cancer, respectively. A substantial number of novel enzymatic reaction\
links were identified. We further validated and explored the interactions of these enzymatic reactions using molecular docking.\
These results highlight the potential of the MPI-VGAE framework for the discovery of novel disease-related enzymatic reactions\
and facilitate the study of the disrupted metabolisms in diseases. </div>', unsafe_allow_html=True)
st.markdown('<p>&nbsp</p>', unsafe_allow_html=True)