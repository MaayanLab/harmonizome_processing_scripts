function outedges = edgesvertcat(inedges1, inedges2)

outedges = struct;

outedges.source = [inedges1.source; inedges2.source];
outedges.sourcename = inedges1.sourcename;

outedges.target = [inedges1.target; inedges2.target];
outedges.targetname = inedges1.targetname;

if isfield(inedges1, 'sourcedesc') && isfield(inedges2, 'sourcedesc')
    
    outedges.sourcedesc = [inedges1.sourcedesc; inedges2.sourcedesc];
    outedges.sourcedescname = inedges1.sourcedescname;

end

if isfield(inedges1, 'targetdesc') && isfield(inedges2, 'targetdesc')
    
    outedges.targetdesc = [inedges1.targetdesc; inedges2.targetdesc];
    outedges.targetdescname = inedges1.targetdescname;
    
end

if isfield(inedges1, 'sourceid') && isfield(inedges2, 'sourceid')
    
    outedges.sourceid = [inedges1.sourceid; inedges2.sourceid];
    outedges.sourceidname = inedges1.sourceidname;

end

if isfield(inedges1, 'targetid') && isfield(inedges2, 'targetid')
    
    outedges.targetid = [inedges1.targetid; inedges2.targetid];
    outedges.targetidname = inedges1.targetidname;
    
end

if isfield(inedges1, 'weight') && isfield(inedges2, 'weight')
    
    outedges.weight = [inedges1.weight; inedges2.weight];

end

outedges.numedges = numel(outedges.source);


