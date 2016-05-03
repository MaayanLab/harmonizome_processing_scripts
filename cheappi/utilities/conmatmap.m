function outcm = conmatmap(incm, outterm, outentry)

if ~isempty(outterm) && ~isempty(outentry)
    
    intercm = struct;
    intercm.term = outterm;
    intercm.termname = incm.termname;
    intercm.numterms = numel(intercm.term);
    
    intercm.entry = incm.entry;
    intercm.entryname = incm.entryname;
    intercm.numentries = numel(intercm.entry);
    
    if isfield(incm, 'entrydesc')
        intercm.entrydesc = incm.entrydesc;
        intercm.entrydescname = incm.entrydescname;
    end
    
    if isfield(incm, 'entryid')
        intercm.entryid = incm.entryid;
        intercm.entryidname = incm.entryidname;
    end
    
    intercm.matrix = zeros([intercm.numterms intercm.numentries]);
    
    [o1, o2] = ismember(intercm.term, incm.term);
    o2(o2 == 0) = [];
    
    intercm.matrix(o1,:) = incm.matrix(o2,:);
    
    if isfield(incm, 'termdesc')
        intercm.termdesc = repmat({'-666'}, intercm.numterms, 1);
        intercm.termdescname = incm.termdescname;
        intercm.termdesc(o1) = incm.termdesc(o2);
    end
    
    if isfield(incm, 'termid')
        intercm.termid = -666*ones([intercm.numterms 1]);
        intercm.termidname = incm.termidname;
        intercm.termid(o1) = incm.termid(o2);
    end
    
    clear incm;
    
    outcm = struct;
    outcm.term = intercm.term;
    outcm.termname = intercm.termname;
    outcm.numterms = numel(outcm.term);
    
    outcm.entry = outentry;
    outcm.entryname = intercm.entryname;
    outcm.numentries = numel(outcm.entry);
    
    if isfield(intercm, 'termdesc')
        outcm.termdesc = intercm.termdesc;
        outcm.termdescname = intercm.termdescname;
    end
    
    if isfield(intercm, 'termid')
        outcm.termid = intercm.termid;
        outcm.termidname = intercm.termidname;
    end
    
    outcm.matrix = zeros([outcm.numterms outcm.numentries]);
    
    [o1, o2] = ismember(outcm.entry, intercm.entry);
    o2(o2 == 0) = [];
    
    outcm.matrix(:,o1) = intercm.matrix(:,o2);
    
    if isfield(intercm, 'entrydesc')
        outcm.entrydesc = repmat({'-666'}, outcm.numentries, 1);
        outcm.entrydescname = intercm.entrydescname;
        outcm.entrydesc(o1) = intercm.entrydesc(o2);
    end
    
    if isfield(intercm, 'entryid')
        outcm.entryid = -666*ones([outcm.numentries 1]);
        outcm.entryidname = intercm.entryidname;
        outcm.entryid(o1) = intercm.entryid(o2);
    end

elseif ~isempty(outterm)
    
    outcm = struct;
    outcm.term = outterm;
    outcm.termname = incm.termname;
    outcm.numterms = numel(outcm.term);
    
    outcm.entry = incm.entry;
    outcm.entryname = incm.entryname;
    outcm.numentries = numel(outcm.entry);
    
    if isfield(incm, 'entrydesc')
        outcm.entrydesc = incm.entrydesc;
        outcm.entrydescname = incm.entrydescname;
    end
    
    if isfield(incm, 'entryid')
        outcm.entryid = incm.entryid;
        outcm.entryidname = incm.entryidname;
    end
    
    outcm.matrix = zeros([outcm.numterms outcm.numentries]);
    
    [o1, o2] = ismember(outcm.term, incm.term);
    o2(o2 == 0) = [];
    
    outcm.matrix(o1,:) = incm.matrix(o2,:);
    
    if isfield(incm, 'termdesc')
        outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
        outcm.termdescname = incm.termdescname;
        outcm.termdesc(o1) = incm.termdesc(o2);
    end
    
    if isfield(incm, 'termid')
        outcm.termid = -666*ones([outcm.numterms 1]);
        outcm.termidname = incm.termidname;
        outcm.termid(o1) = incm.termid(o2);
    end

elseif ~isempty(outentry)
    
    outcm = struct;
    outcm.term = incm.term;
    outcm.termname = incm.termname;
    outcm.numterms = numel(outcm.term);
    
    outcm.entry = outentry;
    outcm.entryname = incm.entryname;
    outcm.numentries = numel(outcm.entry);
    
    if isfield(incm, 'termdesc')
        outcm.termdesc = incm.termdesc;
        outcm.termdescname = incm.termdescname;
    end
    
    if isfield(incm, 'termid')
        outcm.termid = incm.termid;
        outcm.termidname = incm.termidname;
    end
    
    outcm.matrix = zeros([outcm.numterms outcm.numentries]);
    
    [o1, o2] = ismember(outcm.entry, incm.entry);
    o2(o2 == 0) = [];
    
    outcm.matrix(:,o1) = incm.matrix(:,o2);
    
    if isfield(incm, 'entrydesc')
        outcm.entrydesc = repmat({'-666'}, outcm.numentries, 1);
        outcm.entrydescname = incm.entrydescname;
        outcm.entrydesc(o1) = incm.entrydesc(o2);
    end
    
    if isfield(incm, 'entryid')
        outcm.entryid = -666*ones([outcm.numentries 1]);
        outcm.entryidname = incm.entryidname;
        outcm.entryid(o1) = incm.entryid(o2);
    end
    
else
    
    warning('No terms or entries specified.');

end


