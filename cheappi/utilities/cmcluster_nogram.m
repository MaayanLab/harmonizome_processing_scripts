function cm = cmcluster_nogram(cm, showheatmap)
%{
distanceparam = 'cosine';
linkageparam = 'average';

if cm.numterms < 2500
    try
        d = pdist(cm.matrix, distanceparam);
        lnk = linkage(d, linkageparam);
        idxorder = optimalleaforder(lnk, d, 'Criteria', 'adjacent', 'Transformation', 'linear');
        orderedrows = cm.term(idxorder);
        disp('computed optimal leaf order of rows');
    catch e
        disp(getReport(e));
        figure(666);
        clf;
        [~, ~, idxorder] = dendrogram(lnk,0);
        close(gcf);
        orderedrows = cm.term(idxorder);
        disp('computed dendrogram default order of rows');
    end
else
    d = pdist(cm.matrix, distanceparam);
    lnk = linkage(d, linkageparam);
	figure(666);
	clf;
	[~, ~, idxorder] = dendrogram(lnk,0);
	close(gcf);
	orderedrows = cm.term(idxorder);
	disp('computed dendrogram default order of rows');
end

if cm.numentries < 2500
    try
        d = pdist(cm.matrix', distanceparam);
        lnk = linkage(d, linkageparam);
        idxorder = optimalleaforder(lnk, d, 'Criteria', 'adjacent', 'Transformation', 'linear');
        orderedcols = cm.entry(idxorder);
        disp('computed optimal leaf order of columns');
    catch e
        disp(getReport(e));
        figure(666);
        clf;
        [~, ~, idxorder] = dendrogram(lnk,0);
        close(gcf);
        orderedcols = cm.entry(idxorder);
        disp('computed dendrogram default order of columns');
    end
else
    d = pdist(cm.matrix', distanceparam);
    lnk = linkage(d, linkageparam);
	figure(666);
	clf;
	[~, ~, idxorder] = dendrogram(lnk,0);
	close(gcf);
	orderedcols = cm.entry(idxorder);
	disp('computed dendrogram default order of columns');
end

cm = conmatmap(cm, orderedrows, orderedcols);

if showheatmap
    HeatMap(cm.matrix, 'Colormap', redbluecmap);
end
%}

sm = cm2sm_cosine(cm, false);
orderedrows = sm.term;
clear sm;
sm = cm2sm_cosine(cmtranspose(cm), false);
orderedcols = sm.term;
clear sm;
cm = conmatmap(cm, orderedrows, orderedcols);
if showheatmap
    HeatMap(cm.matrix, 'Colormap', redbluecmap);
end
%}

