import hotspot
import pickle
import pandas as pd
import plotly.express as px
import plotly as py
from plotnine import *
import os,sys,re

in_hotspot=sys.argv[1]  #hotspot pickle object
outdir=sys.argv[2]  #outdir

subdir=['module_scores','module_genes',"module_3d_plot","module_2d_plot"]
for i in subdir:
    if not os.path.exists(outdir+"/"+i):
        os.mkdir(outdir+"/"+i)

def pickleload(file):
    ff = open(file,'rb')
    obj = pickle.load(ff)
    ff.close()
    return obj
hs=pickleload(in_hotspot)

#Firstly conpute module scores
module_scores = hs.calculate_module_scores()

#Export genes for each module
#Export module scores and plots
for module in range(1, hs.modules.max()+1):
    scores = hs.module_scores[module]
    plotly_df=pd.DataFrame({'x':hs.latent.iloc[:, 0],'y':hs.latent.iloc[:, 1],'z':hs.latent.iloc[:, 2],'scores':scores})
    plotly_df.to_csv(outdir+"/module_scores/module_"+str(module)+"_scores.csv")
    fig = px.scatter_3d(plotly_df, x='x', y='y', z='z',
                  color='scores',color_continuous_scale=px.colors.sequential.Viridis)
    fig.update_layout(scene_aspectmode='data')
    py.offline.plot(fig,filename=outdir+"/module_3d_plot/module_"+str(module)+"_3d.html")
    
    #plot in 2d
    plotly_df['slice_num']=(plotly_df['y'].astype(int)/7+1).astype(int)

    p=ggplot(aes(x="x", y="y"), plotly_df)+geom_point(aes(color='scores'),size=1)+facet_wrap('slice_num',ncol=5)+theme_void()
    ggsave(p,outdir+"/"+"module_2d_plot/module_"+str(module)+"_2d.png",width=15,height=15)
    ggsave(p,outdir+"/"+"module_2d_plot/module_"+str(module)+"_2d.pdf",width=15,height=15)

    #Export genes
    results = hs.results.join(hs.modules)
    results = results.loc[results.Module == module]
    results = results.sort_values('Z', ascending=False)
    results.to_csv(outdir+"/module_genes/module_"+str(module)+"_genes.csv")
