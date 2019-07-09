taylor_data = {}

for key in model_keys:
    taylor_data.update({key:[model.std(ddof=1), np.corrcoef(data, model)[0,1]]})


###################################################################
from pynmd.plotting import taylor
markersize = 6
fig = plt.figure(6,figsize=(9,9))
fig.clf()

refstd = data.std(ddof=1)

# Taylor diagram
dia = taylor.TaylorDiagram(refstd, fig=fig, rect=111, label="Reference")
 
colors = plt.matplotlib.cm.jet(np.linspace(0,1,len(taylor_data.keys())))
 
# Add samples to Taylor diagram
for imodel in range(len(taylor_data.keys())):
    key      = model_keys[imodel]
    #key      = taylor_data.keys()[imodel]
    stddev   = taylor_data[key][0]
    corrcoef = taylor_data[key][1]
    marker   = ps.marker[imodel]
    dia.add_sample(stddev, corrcoef,ref=False, marker=marker, ls='', c=ps.colors[imodel],
               markersize=markersize,label=key)
     
# add refrence point for data     
dia.add_sample(refstd, 1.0 ,ref=True, marker='*', ls='',  c='k',
               markersize=markersize*1.5,label='Ref.')
 
# Add RMS contours, and label them
contours = dia.add_contours(levels=8,data_std=refstd,colors='0.5')
plt.clabel(contours, inline=1,fmt = '%.3g', fontsize=10)
 
# Add a figure legend
if True:
    leg2=fig.legend(dia.samplePoints,
               [ p.get_label() for p in dia.samplePoints ],
               numpoints=1, prop=dict(size='small'), loc='upper right',
               ncol=2)
 
    frame=leg2.get_frame()
    frame.set_edgecolor('None')
    #frame.set_facecolor('None')
    frame.set_facecolor('w')

plt.title( 'HWM data' + '  N=' + str(len(model)), position=(0.1, 1.04))
plt.subplots_adjust(left=left1, bottom=bottom1, right= right1, top= top1,
      wspace=wspace1, hspace=hspace1)
plt.savefig(out_dir+ '/taylor_HWM' + ftype,dpi=dpi)
plt.close('all')

