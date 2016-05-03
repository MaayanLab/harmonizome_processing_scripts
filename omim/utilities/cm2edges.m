function edges = cm2edges(cm)

[i, j] = find(cm.matrix ~= 0);
k = find(cm.matrix ~= 0); % or k = sub2ind([cm.numterms cm.numentries], i, j);

edges = struct;

edges.source = cm.term(i);
edges.sourcename = cm.termname;

if isfield(cm, 'termdesc')
    edges.sourcedesc = cm.termdesc(i);
    edges.sourcedescname = cm.termdescname;
end

if isfield(cm, 'termid')
    edges.sourceid = cm.termid(i);
    edges.sourceidname = cm.termidname;
end

edges.target = cm.entry(j);
edges.targetname = cm.entryname;

if isfield(cm, 'entrydesc')
    edges.targetdesc = cm.entrydesc(j);
    edges.targetdescname = cm.entrydescname;
end

if isfield(cm, 'entryid')
    edges.targetid = cm.entryid(j);
    edges.targetidname = cm.entryidname;
end

edges.numedges = numel(edges.source);

if abs((max(cm.matrix(k)) - min(cm.matrix(k)))/min(cm.matrix(k))) >= 0.001
    edges.weight = cm.matrix(k);
    [~, si] = sort(edges.weight, 'descend');
    edges = edgesorder(edges, si);
end


