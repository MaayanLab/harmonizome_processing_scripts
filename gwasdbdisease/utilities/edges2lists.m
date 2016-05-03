function lists = edges2lists(edges)

lists = struct;

[lists.term, ui, ri] = unique(edges.source);
lists.termname = edges.sourcename;

if isfield(edges, 'sourcedesc')
    lists.termdesc = edges.sourcedesc(ui);
    lists.termdescname = edges.sourcedescname;
end

if isfield(edges, 'sourceid')
    lists.termid = edges.sourceid(ui);
    lists.termidname = edges.sourceidname;
end

lists.numterms = numel(lists.term);

if isfield(edges, 'targetdesc') && isfield(edges, 'targetid') && isfield(edges, 'weight')
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = edges.targetdescname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = edges.targetidname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);
        entrydescs = edges.targetdesc(ri == i);
        entryids = edges.targetid(ri == i);
        weights = edges.weight(ri == i);

        [weights, si] = sort(weights, 'descend');
        entries = entries(si);
        entrydescs = entrydescs(si);
        entryids = entryids(si);

        [lists.entries{i}, ui, ~] = unique(entries, 'stable');
        lists.entrydescs{i} = entrydescs(ui);
        lists.entryids{i} = entryids(ui);
        lists.weights{i} = weights(ui);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
elseif isfield(edges, 'targetdesc') && isfield(edges, 'weight')
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = edges.targetdescname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);
        entrydescs = edges.targetdesc(ri == i);
        weights = edges.weight(ri == i);

        [weights, si] = sort(weights, 'descend');
        entries = entries(si);
        entrydescs = entrydescs(si);

        [lists.entries{i}, ui, ~] = unique(entries, 'stable');
        lists.entrydescs{i} = entrydescs(ui);
        lists.weights{i} = weights(ui);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
elseif isfield(edges, 'targetid') && isfield(edges, 'weight')
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = edges.targetidname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);
        entryids = edges.targetid(ri == i);
        weights = edges.weight(ri == i);

        [weights, si] = sort(weights, 'descend');
        entries = entries(si);
        entryids = entryids(si);

        [lists.entries{i}, ui, ~] = unique(entries, 'stable');
        lists.entryids{i} = entryids(ui);
        lists.weights{i} = weights(ui);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
elseif isfield(edges, 'targetdesc') && isfield(edges, 'targetid')
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = edges.targetdescname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = edges.targetidname;
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);
        entrydescs = edges.targetdesc(ri == i);
        entryids = edges.targetid(ri == i);

        [lists.entries{i}, ui, ~] = unique(entries);
        lists.entrydescs{i} = entrydescs(ui);
        lists.entryids{i} = entryids(ui);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
elseif isfield(edges, 'targetdesc')
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = edges.targetdescname;
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);
        entrydescs = edges.targetdesc(ri == i);

        [lists.entries{i}, ui, ~] = unique(entries);
        lists.entrydescs{i} = entrydescs(ui);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
elseif isfield(edges, 'targetid')
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = edges.targetidname;
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);
        entryids = edges.targetid(ri == i);

        [lists.entries{i}, ui, ~] = unique(entries);
        lists.entryids{i} = entryids(ui);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
elseif isfield(edges, 'weight')
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);
        weights = edges.weight(ri == i);

        [weights, si] = sort(weights, 'descend');
        entries = entries(si);

        [lists.entries{i}, ui, ~] = unique(entries, 'stable');
        lists.weights{i} = weights(ui);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
else
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = edges.targetname;
    lists.numentries = zeros([lists.numterms 1]);

    for i = 1:1:lists.numterms

        entries = edges.target(ri == i);

        lists.entries{i} = unique(entries);

        lists.numentries(i) = numel(lists.entries{i});

    end
    
end


