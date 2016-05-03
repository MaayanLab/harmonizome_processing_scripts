function edges = edgesunique(edges)

if isfield(edges, 'weight')
    
    [~, si] = sort(edges.weight, 'descend');
    
    edges = edgesorder(edges, si);
    
end

edge = sourcetargetcat(edges.source, edges.target);

[~, ui, ~] = unique(edge, 'stable');

discard = true([edges.numedges 1]);
discard(ui) = false;

edges = edgesdiscard(edges, discard);


