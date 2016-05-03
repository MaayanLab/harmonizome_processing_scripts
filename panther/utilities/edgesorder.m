function edges = edgesorder(edges, indexorder)

edges.source = edges.source(indexorder);

edges.target = edges.target(indexorder);

if isfield(edges, 'weight')
    edges.weight = edges.weight(indexorder);
end

if isfield(edges, 'sourcedesc')
    edges.sourcedesc = edges.sourcedesc(indexorder);
end

if isfield(edges, 'sourceid')
    edges.sourceid = edges.sourceid(indexorder);
end

if isfield(edges, 'targetdesc')
    edges.targetdesc = edges.targetdesc(indexorder);
end

if isfield(edges, 'targetid')
    edges.targetid = edges.targetid(indexorder);
end

if isfield(edges, 'overlap')
    edges.overlap = edges.overlap(indexorder);
end

if isfield(edges, 'label')
    edges.label = edges.label(indexorder);
end

if isfield(edges, 'isforlearning')
    edges.isforlearning = edges.isforlearning(indexorder);
end

if isfield(edges, 'targetdesc2')
    edges.targetdesc2 = edges.targetdesc2(indexorder);
end

if isfield(edges, 'predictedlabel')
    edges.predictedlabel = edges.predictedlabel(indexorder);
end

if isfield(edges, 'label123')
    edges.label123 = edges.label123(indexorder);
end

if isfield(edges, 'predictedlabel123')
    edges.predictedlabel123 = edges.predictedlabel123(indexorder);
end


