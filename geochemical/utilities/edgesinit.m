function edges = edgesinit(numedges, source, sourcename, sourcedesc, sourcedescname, sourceid, sourceidname, target, targetname, targetdesc, targetdescname, targetid, targetidname, isweighted, weight)

edges = struct;

if ~isempty(source)
    edges.source = source;
else
    edges.source = repmat({'-666'}, numedges, 1);
end
edges.sourcename = sourcename;

if ~isempty(sourcedescname)
    if ~isempty(sourcedesc)
        edges.sourcedesc = sourcedesc;
    else
        edges.sourcedesc = repmat({'-666'}, numedges, 1);
    end
    edges.sourcedescname = sourcedescname;
end

if ~isempty(sourceidname)
    if ~isempty(sourceid)
        edges.sourceid = sourceid;
    else
        edges.sourceid = -666*ones([numedges 1]);
    end
    edges.sourceidname = sourceidname;
end

edges.numedges = numedges;

if ~isempty(target)
    edges.target = target;
else
    edges.target = repmat({'-666'}, numedges, 1);
end
edges.targetname = targetname;

if ~isempty(targetdescname)
    if ~isempty(targetdesc)
        edges.targetdesc = targetdesc;
    else
        edges.targetdesc = repmat({'-666'}, numedges, 1);
    end
    edges.targetdescname = targetdescname;
end

if ~isempty(targetidname)
    if ~isempty(targetid)
        edges.targetid = targetid;
    else
        edges.targetid = -666*ones([numedges 1]);
    end
    edges.targetidname = targetidname;
end

edges.numedges = numedges;

if isweighted
    if ~isempty(weight)
        edges.weight = weight;
    else
        edges.weight = zeros([numedges 1]);
    end
end


