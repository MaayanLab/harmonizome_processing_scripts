function outedges = edges2symedges(inedges)

outedges = struct;

outedges.source = [inedges.source; inedges.target];
outedges.sourcename = inedges.sourcename;

outedges.target = [inedges.target; inedges.source];
outedges.targetname = inedges.targetname;

if isfield(inedges, 'sourcedesc') && isfield(inedges, 'targetdesc')
    
    outedges.sourcedesc = [inedges.sourcedesc; inedges.targetdesc];
    outedges.sourcedescname = inedges.sourcedescname;
    
    outedges.targetdesc = [inedges.targetdesc; inedges.sourcedesc];
    outedges.targetdescname = inedges.targetdescname;
    
end

if isfield(inedges, 'sourceid') && isfield(inedges, 'targetid')
    
    outedges.sourceid = [inedges.sourceid; inedges.targetid];
    outedges.sourceidname = inedges.sourceidname;
    
    outedges.targetid = [inedges.targetid; inedges.sourceid];
    outedges.targetidname = inedges.targetidname;
    
end

outedges.numedges = numel(outedges.source);

if isfield(inedges, 'weight')
    
    outedges.weight = abs([inedges.weight; inedges.weight]);
    
    [~, si] = sort(outedges.weight, 'descend');
    outedges = edgesorder(outedges, si);
    
end

edge = sourcetargetcat(outedges.source, outedges.target);
[~, ui, ~] = unique(edge, 'stable');
discard = ~ismember((1:1:outedges.numedges)', ui);
outedges = edgesdiscard(outedges, discard);


