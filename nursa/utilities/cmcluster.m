function [cm, cgo] = cmcluster(cm, showheatmap)

standardizeparam = 'none';
clusterparam = 'all';
distanceparam = 'cosine';
linkageparam = 'average';

cgo = cm2clustergram(cm, standardizeparam, clusterparam, distanceparam, linkageparam);

cm = conmatmap(cm, cgo.RowLabels, cgo.ColumnLabels');

if showheatmap
    HeatMap(cm.matrix, 'Colormap', redbluecmap);
end


