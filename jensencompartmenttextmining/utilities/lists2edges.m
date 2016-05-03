function edges = lists2edges(lists)

edges = struct;

if isfield(lists, 'termdesc') && isfield(lists, 'termid')
    
    source = cell(lists.numterms, 1);
    sourcedesc = cell(lists.numterms, 1);
    sourceid = cell(lists.numterms, 1);
    
    for i = 1:1:lists.numterms

        source{i} = repmat(lists.term(i), lists.numentries(i), 1);
        sourcedesc{i} = repmat(lists.termdesc(i), lists.numentries(i), 1);
        sourceid{i} = lists.termid(i)*ones([lists.numentries(i) 1]);

    end

    edges.source = vertcat(source{:});
    edges.sourcename = lists.termname;
    edges.sourcedesc = vertcat(sourcedesc{:});
    edges.sourcedescname = lists.termdescname;
    edges.sourceid = vertcat(sourceid{:});
    edges.sourceidname = lists.termidname;
    
elseif isfield(lists, 'termdesc')
    
    source = cell(lists.numterms, 1);
    sourcedesc = cell(lists.numterms, 1);
    
    for i = 1:1:lists.numterms

        source{i} = repmat(lists.term(i), lists.numentries(i), 1);
        sourcedesc{i} = repmat(lists.termdesc(i), lists.numentries(i), 1);

    end

    edges.source = vertcat(source{:});
    edges.sourcename = lists.termname;
    edges.sourcedesc = vertcat(sourcedesc{:});
    edges.sourcedescname = lists.termdescname;
    
elseif isfield(lists, 'termid')
    
    source = cell(lists.numterms, 1);
    sourceid = cell(lists.numterms, 1);
    
    for i = 1:1:lists.numterms

        source{i} = repmat(lists.term(i), lists.numentries(i), 1);
        sourceid{i} = lists.termid(i)*ones([lists.numentries(i) 1]);

    end

    edges.source = vertcat(source{:});
    edges.sourcename = lists.termname;
    edges.sourceid = vertcat(sourceid{:});
    edges.sourceidname = lists.termidname;
    
else
    
    source = cell(lists.numterms, 1);
    
    for i = 1:1:lists.numterms

        source{i} = repmat(lists.term(i), lists.numentries(i), 1);

    end

    edges.source = vertcat(source{:});
    edges.sourcename = lists.termname;
    
end

edges.target = vertcat(lists.entries{:});
edges.targetname = lists.entryname;

if isfield(lists, 'entrydescs')
    
    edges.targetdesc = vertcat(lists.entrydescs{:});
    edges.targetdescname = lists.entrydescname;
    
end

if isfield(lists, 'entryids')
    
    edges.targetid = vertcat(lists.entryids{:});
    edges.targetidname = lists.entryidname;
    
end

if isfield(lists, 'weights')
    
    edges.weight = vertcat(lists.weights{:});
    
    [edges.weight, si] = sort(edges.weight, 'descend');
    
    edges.source = edges.source(si);
    
    if isfield(edges, 'sourcedesc')
        edges.sourcedesc = edges.sourcedesc(si);
    end
    
    if isfield(edges, 'sourceid')
        edges.sourceid = edges.sourceid(si);
    end
    
    edges.target = edges.target(si);
    
    if isfield(edges, 'targetdesc')
        edges.targetdesc = edges.targetdesc(si);
    end
    
    if isfield(edges, 'targetid')
        edges.targetid = edges.targetid(si);
    end
    
end

edges.numedges = numel(edges.source);


