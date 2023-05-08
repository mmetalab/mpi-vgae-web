######################
# Import libraries
######################
import numpy as np
import pandas as pd
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

image = Image.open('./images/about.png')

st.image(image, use_column_width=True)

# st.markdown('<div style="text-align: justify;">Here we developed a XgbMPNN model\
# that integrates message passing attention network (MPNN) and Xgboosting algorithm\
# to predict CCS values of metabolites and lipids based on the chemical structures.\
# A comprehensive CCS dataset containing over 6000 lipids and metabolites were used\
# to train and optimize the model. We compared the XgbMPNN with some conventional\
# machine learning models, and it showed superior performance in CCS value prediction with\
# median relative error (MRE) of 1.2% and 2.1% in lipids and metabolites, respectively.\
# This graph neural network-based model boosts the CCS value prediction performance and improves\
# the efficiency for small molecule identification based on IM-MS. </div>', unsafe_allow_html=True)
# st.markdown('<p>&nbsp</p>', unsafe_allow_html=True)