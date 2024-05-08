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

st.sidebar.info("Model *two countries* applied to Germany")

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
test test test
''')
#=====================================================#
#=====================================================#

def plot_stack(df,col=None,method="shapiro",robust=None,unclassified:bool=True,year:int=2015):
            
    def yrindex(ind):
        x = []
        for el in ind:
            if el.month==1:
                x.append(el.year)
            else:
                x.append("")
        return x

    if method=="shapiro":
        color = ['cornflowerblue','indianred']
        if robust==None:
            #cols = self.meth["base"].copy()
            cols = col["base"].copy()
        else:
            #cols = self.meth[robust].copy()
            cols = col[robust].copy()
        leg = ['Demand','Supply']
    elif method=="perso":
        cols = ["dem_pers","sup_pers","dem_trans","sup_trans","dem_abg","sup_abg"]
        leg = ["Persistent demand","Persistent supply","Transitory demand","Transitory supply","Ambiguous demand","Ambiguous supply"]
        color = ['mediumseagreen','brown','royalblue','orange',"yellowgreen","darksalmon"]
    else:
        if robust=="complex":
            cols = ["dem_pers","sup_pers","dem_trans","sup_trans"]
            leg = ["Persistent demand","Persistent supply","Transitory demand","Transitory supply"]
            color = ['mediumseagreen','brown','royalblue','orange']
        else:
            cols = ["dem","sup"]
            leg = ['Demand','Supply']
            color = ['cornflowerblue','indianred']

    if unclassified==True:
        cols.append('unclassified')
        leg.append('Unclassified')
        color.append('darkgray')
                
    f = df[df.index.year>=year][cols+['total']]
    f = f.rename(columns={cols[i]:leg[i] for i in range(len(cols))})
    f = f.reset_index()
    
    fig = px.bar(f,y=leg,x="index",labels={"index":"year","value":"%YoY HICP"},color_discrete_sequence=color)
    fig.add_trace(trace=go.Scatter(x=f["index"],y=f["total"], name="Total",mode='lines+markers', line_color="black"))
    #https://stackoverflow.com/questions/58188816/change-line-color-in-plotly
    return fig

#=====================================================#
#=====================================================#
st.markdown('# Applied Macroeconomic Modelling results report')