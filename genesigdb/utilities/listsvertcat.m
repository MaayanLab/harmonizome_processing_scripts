function outlists = listsvertcat(inlists1, inlists2)

outlists = struct;

outlists.term = [inlists1.term; inlists2.term];
outlists.termname = inlists1.termname;

outlists.entries = [inlists1.entries; inlists2.entries];
outlists.entryname = inlists1.entryname;

outlists.numentries = [inlists1.numentries; inlists2.numentries];

if isfield(inlists1, 'termdesc') && isfield(inlists2, 'termdesc')
    
    outlists.termdesc = [inlists1.termdesc; inlists2.termdesc];
    outlists.termdescname = inlists1.termdescname;

end

if isfield(inlists1, 'entrydescs') && isfield(inlists2, 'entrydescs')
    
    outlists.entrydescs = [inlists1.entrydescs; inlists2.entrydescs];
    outlists.entrydescname = inlists1.entrydescname;
    
end

if isfield(inlists1, 'termid') && isfield(inlists2, 'termid')
    
    outlists.termid = [inlists1.termid; inlists2.termid];
    outlists.termidname = inlists1.termidname;

end

if isfield(inlists1, 'entryids') && isfield(inlists2, 'entryids')
    
    outlists.entryids = [inlists1.entryids; inlists2.entryids];
    outlists.entryidname = inlists1.entryidname;
    
end

if isfield(inlists1, 'weights') && isfield(inlists2, 'weights')
    
    outlists.weights = [inlists1.weights; inlists2.weights];

end

outlists.numterms = numel(outlists.term);


