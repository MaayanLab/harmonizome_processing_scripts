function cm = cmrowdiscard(cm, discard)

cm.term(discard) = [];

if isfield(cm, 'termdesc')
    cm.termdesc(discard) = [];
end

if isfield(cm, 'termid')
    cm.termid(discard) = [];
end

cm.matrix(discard,:) = [];

if isfield(cm, 'stretchedmatrix')
    cm.stretchedmatrix(discard,:) = [];
end

if isfield(cm, 'weightmatrix')
    cm.weightmatrix(discard,:) = [];
end

if isfield(cm, 'chdist')
    cm = rmfield(cm, 'chdist');
end

cm.numterms = numel(cm.term);


