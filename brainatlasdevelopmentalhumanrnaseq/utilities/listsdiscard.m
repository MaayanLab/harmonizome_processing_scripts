function lists = listsdiscard(lists, discard)

lists.term(discard) = [];

if isfield(lists, 'termdesc')
    lists.termdesc(discard) = [];
end

if isfield(lists, 'termid')
    lists.termid(discard) = [];
end

lists.entries(discard) = [];

lists.numentries(discard) = [];

if isfield(lists, 'entrydescs')
    lists.entrydescs(discard) = [];
end

if isfield(lists, 'entryids')
    lists.entryids(discard) = [];
end

if isfield(lists, 'weights')
    lists.weights(discard) = [];
end

lists.numterms = numel(lists.term);


