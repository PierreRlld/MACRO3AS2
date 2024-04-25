#=====================================================#
import streamlit as st
import plotly.tools as tls
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
#=====================================================#

#? ---
#? cd main
#? streamlit run streamlit_app.py
#? ---

st.set_page_config(layout='wide', initial_sidebar_state='expanded')

with open(Path(__file__).parents[0] / "style.css") as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
    
#=====================================================#
#=====================================================#

st.sidebar.header('Model parameters')

st.sidebar.info("Classifier can take up to 1 minute")

country = st.sidebar.selectbox('Country', ('EU27','France','Germany','Spain'),index=None) 
unclass = st.sidebar.selectbox('Show unclassified', (True,False), index=0)
year = st.sidebar.selectbox('Chart start year date:',[i for i in range(2023,2000,-1)],index=8)

st.sidebar.markdown('''---''')     

order = st.sidebar.selectbox('VAR order', ("auto","fixed"), index=0)
if order=="auto":
    st.sidebar.markdown('<p style="font-size:80%;">Max lag search was set at 24</p>', unsafe_allow_html=True) 
elif order=="fixed":
    order = st.sidebar.text_input('Input VAR desired order:')
    try:
        order = int(order)
        assert 1<=order<=24
    except:
        order = None
        st.sidebar.error('Enter valid integer') #icon="🚨"
        
robust = st.sidebar.selectbox('Robust methods', (True,False), index=1)       
if robust==True:
    #st.sidebar.markdown('''---''') 
    shap_rob_plot = st.sidebar.selectbox('Plot Shapiro robust method', ['j1','j2','j3','param'], index=3)
else:
    shap_rob_plot = None
 
st.sidebar.markdown('''
---
Model : `two countries`
''')
#=====================================================#
#=====================================================#



#=====================================================#
#=====================================================#
st.markdown('# Applied Macroeconomic Modelling results report')