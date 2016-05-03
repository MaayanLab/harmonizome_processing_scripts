function edges = cm2edges_keepzeros(cm)

[J, I] = meshgrid(1:1:cm.numentries, 1:1:cm.numterms);
j = J(:);
i = I(:);
clear J I;

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

edges.weight = cm.matrix(:);
[~, si] = sort(edges.weight, 'descend');
edges = edgesorder(edges, si);


