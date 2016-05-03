function cm = edges2cm(edges)

cm = struct;


[cm.term, ui, ~] = unique(edges.source);
cm.termname = edges.sourcename;

if isfield(edges, 'sourcedesc')
    cm.termdesc = edges.sourcedesc(ui);
    cm.termdescname = edges.sourcedescname;
end

if isfield(edges, 'sourceid')
    cm.termid = edges.sourceid(ui);
    cm.termidname = edges.sourceidname;
end

cm.numterms = numel(cm.term);


[cm.entry, ui, ~] = unique(edges.target);
cm.entryname = edges.targetname;

if isfield(edges, 'targetdesc')
    cm.entrydesc = edges.targetdesc(ui);
    cm.entrydescname = edges.targetdescname;
end

if isfield(edges, 'targetid')
    cm.entryid = edges.targetid(ui);
    cm.entryidname = edges.targetidname;
end

cm.numentries = numel(cm.entry);


cm.matrix = zeros([cm.numterms cm.numentries]);

[~, i] = ismember(edges.source, cm.term);
[~, j] = ismember(edges.target, cm.entry);
k = sub2ind([cm.numterms cm.numentries], i, j);

if isfield(edges, 'weight')
    
    [~, si] = sort(edges.weight, 'ascend'); % this makes sure the max value gets placed into cm.matrix(k) if multiple edges map to cm.matrix(k)
    cm.matrix(k(si)) = edges.weight(si);
    
else
    
    cm.matrix(k) = 1;
    
end


