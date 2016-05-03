function cm = cmcoldiscard(cm, discard)

cm.entry(discard) = [];

if isfield(cm, 'entrydesc')
    cm.entrydesc(discard) = [];
end

if isfield(cm, 'entryid')
    cm.entryid(discard) = [];
end

cm.matrix(:,discard) = [];

if isfield(cm, 'stretchedmatrix')
    cm.stretchedmatrix(:,discard) = [];
end

if isfield(cm, 'weightmatrix')
    cm.weightmatrix(:,discard) = [];
end

if isfield(cm, 'chdist')
    cm = rmfield(cm, 'chdist');
end

cm.numentries = numel(cm.entry);


