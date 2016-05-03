function cm = cminit(numterms, numentries, term, termname, termdesc, termdescname, termid, termidname, entry, entryname, entrydesc, entrydescname, entryid, entryidname, matrix)

cm = struct;

if ~isempty(term)
    cm.term = term;
else
    cm.term = repmat({'-666'}, numterms, 1);
end
cm.termname = termname;

if ~isempty(termdescname)
    if ~isempty(termdesc)
        cm.termdesc = termdesc;
    else
        cm.termdesc = repmat({'-666'}, numterms, 1);
    end
    cm.termdescname = termdescname;
end

if ~isempty(termidname)
    if ~isempty(termid)
        cm.termid = termid;
    else
        cm.termid = -666*ones([numterms 1]);
    end
    cm.termidname = termidname;
end

cm.numterms = numterms;

if ~isempty(entry)
    cm.entry = entry;
else
    cm.entry = repmat({'-666'}, numentries, 1);
end
cm.entryname = entryname;

if ~isempty(entrydescname)
    if ~isempty(entrydesc)
        cm.entrydesc = entrydesc;
    else
        cm.entrydesc = repmat({'-666'}, numentries, 1);
    end
    cm.entrydescname = entrydescname;
end

if ~isempty(entryidname)
    if ~isempty(entryid)
        cm.entryid = entryid;
    else
        cm.entryid = -666*ones([numentries 1]);
    end
    cm.entryidname = entryidname;
end

cm.numentries = numentries;

if ~isempty(matrix)
    cm.matrix = matrix;
else
    cm.matrix = zeros([numterms numentries]);
end


