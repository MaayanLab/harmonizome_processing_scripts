function outedges = sourcetargetswap(inedges)

outedges = struct;

outedges.source = inedges.target;
outedges.sourcename = inedges.targetname;

if isfield(inedges, 'targetdesc')
    outedges.sourcedesc = inedges.targetdesc;
    outedges.sourcedescname = inedges.targetdescname;
end

if isfield(inedges, 'targetid')
    outedges.sourceid = inedges.targetid;
    outedges.sourceidname = inedges.targetidname;
end

outedges.target = inedges.source;
outedges.targetname = inedges.sourcename;

if isfield(inedges, 'sourcedesc')
    outedges.targetdesc = inedges.sourcedesc;
    outedges.targetdescname = inedges.sourcedescname;
end

if isfield(inedges, 'sourceid')
    outedges.targetid = inedges.sourceid;
    outedges.targetidname = inedges.sourceidname;
end

if isfield(inedges, 'weight')
    outedges.weight = inedges.weight;
end

outedges.numedges = inedges.numedges;