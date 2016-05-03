function lists = listsinit(numterms, term, termname, termdesc, termdescname, termid, termidname, entries, entryname, entrydescs, entrydescname, entryids, entryidname, isweighted, weights)

lists = struct;

if ~isempty(term)
    lists.term = term;
else
    lists.term = repmat({'-666'}, numterms, 1);
end
lists.termname = termname;

if ~isempty(termdescname)
    if ~isempty(termdesc)
        lists.termdesc = termdesc;
    else
        lists.termdesc = repmat({'-666'}, numterms, 1);
    end
    lists.termdescname = termdescname;
end

if ~isempty(termidname)
    if ~isempty(termid)
        lists.termid = termid;
    else
        lists.termid = -666*ones([numterms 1]);
    end
    lists.termidname = termidname;
end

lists.numterms = numterms;
lists.numentries = zeros([numterms 1]);

if ~isempty(entries)
    lists.entries = entries;
    lists.numentries = cellfun(@numel, lists.entries);
else
    lists.entries = repmat({{'-666'}}, numterms, 1);
end
lists.entryname = entryname;

if ~isempty(entrydescname)
    if ~isempty(entrydescs)
        lists.entrydescs = entrydescs;
        lists.numentries = cellfun(@numel, lists.entrydescs);
    else
        lists.entrydescs = repmat({{'-666'}}, numterms, 1);
    end
    lists.entrydescname = entrydescname;
end

if ~isempty(entryidname)
    if ~isempty(entryids)
        lists.entryids = entryids;
        lists.numentries = cellfun(@numel, lists.entryids);
    else
        lists.entryids = repmat({-666}, numterms, 1);
    end
    lists.entryidname = entryidname;
end

if isweighted
    if ~isempty(weights)
        lists.weights = weights;
        lists.numentries = cellfun(@numel, lists.weights);
    else
        lists.weights = repmat({-666}, numterms, 1);
    end
end


