function edges = edgesdiscard(edges, discard)

edges.source(discard) = [];

if isfield(edges, 'sourcedesc')
    edges.sourcedesc(discard) = [];
end

if isfield(edges, 'sourceid')
    edges.sourceid(discard) = [];
end

edges.target(discard) = [];

if isfield(edges, 'targetdesc')
    edges.targetdesc(discard) = [];
end

if isfield(edges, 'targetid')
    edges.targetid(discard) = [];
end

if isfield(edges, 'weight')
    edges.weight(discard,:) = [];
end

if isfield(edges, 'overlap')
    edges.overlap(discard) = [];
end

if isfield(edges, 'label')
    edges.label(discard) = [];
end

if isfield(edges, 'isforlearning')
    edges.isforlearning(discard) = [];
end

if isfield(edges, 'features')
    edges.features(discard,:) = [];
end

if isfield(edges, 'predictedlabel')
    edges.predictedlabel(discard,:) = [];
end

if isfield(edges, 'predictedlabel123')
    edges.predictedlabel123(discard) = [];
end

if isfield(edges, 'targetdesc2')
    edges.targetdesc2(discard) = [];
end

if isfield(edges, 'label123')
    edges.label123(discard) = [];
end

edges.numedges = numel(edges.source);


