function sm = cm2sm_cosine(cm, show)

sm = struct;

sm.term = cm.term;
sm.termname = cm.termname;

if isfield(cm, 'termid');
    sm.termid = cm.termid;
    sm.termidname = cm.termidname;
end

if isfield(cm, 'termdesc');
    sm.termdesc = cm.termdesc;
    sm.termdescname = cm.termdescname;
end

sm.numterms = cm.numterms;

mags = sqrt(diag(cm.matrix*cm.matrix'));
sm.matrix = (cm.matrix*cm.matrix')./(mags*mags');
sm.matrix(sm.matrix > 1) = 1;
sm.matrix(sm.matrix < -1) = -1;

sm = smcluster(sm);

if show
    
    HeatMap(sm.matrix, 'RowLabels', sm.term, 'ColumnLabels', sm.term, 'Colormap', redbluecmap);
%     HeatMap(sm.scaledmatrix, 'RowLabels', sm.term, 'ColumnLabels', sm.term, 'Colormap', redbluecmap);
    
end


