import pandas as pd
import plotly_express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import datetime
 
# using now() to get current time
current_time = datetime.datetime.now()
print(current_time)

# making dataframe 
df = pd.read_csv("test.csv") 

#group by target and kinematic_type, then do cumulative sum of charge
print(df.groupby(['target', 'kin_type'])['charge'].cumsum())
charge_csum = df.groupby(['target', 'kin_type'])['charge'].cumsum()

#fig = go.Figure()
fig = make_subplots(specs=[[{"secondary_y": True}]])

fig.add_traces(list(px.line(df, x='run', y=charge_csum, title='cumulative charge', color='target', line_dash='kin_type', markers=True).select_traces()))
fig.add_traces(list(px.bar(df, x='run', y='charge', color='target', pattern_shape="kin_type").select_traces()))
fig.update_yaxes(title_text="charge [mC]", secondary_y=False)

fig.update_layout(hovermode='x unified')
#fig.update_xaxes(tickmode='linear', ticklabelmode="period")
fig.update_layout(bargap=0)

#fig.add_trace(px.bar(df, x='run', y='charge', color='target', pattern_shape="kin_type").data[1])
#fig.update_layout(xaxis={'visible': False, 'showticklabels': False})
fig.show()


#min_run=min(df['run'])
#max_run=max(df['run'])

#df_LD2_MF = df[(df['target'].str.contains('LD2')) & (df['kin_type'].str.contains('MF')) ]
#df_LD2_SRC = df[(df['target'].str.contains('LD2')) & (df['kin_type'].str.contains('SRC')) ]

# get cumulative charges per target and per kin type (MF or SRC)
#print(df_LD2_MF['charge'].cumsum())

#total_bins=max_run - min_run + 1

#plotting the histogram
#hist = px.histogram(df, x="run", y="charge", color='target', nbins=total_bins, hover_data=df.columns)
#line = px.line(df_LD2_MF, x="run", y="charge", color="target")

#hist.update_layout(bargap=0.2)

#fig = make_subplots(specs=[[{"secondary_y": True}]])
#fig.add_histogram(hist)
# showing the plot
#fig.show()
