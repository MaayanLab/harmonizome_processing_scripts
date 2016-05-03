function lists = cm2lists(cm)

lists = struct;

lists.term = cm.term;
lists.termname = cm.termname;

if isfield(cm, 'termdesc')
    lists.termdesc = cm.termdesc;
    lists.termdescname = cm.termdescname;
end

if isfield(cm, 'termid')
    lists.termid = cm.termid;
    lists.termidname = cm.termidname;
end

lists.numterms = cm.numterms;

if isfield(cm, 'entrydesc') && isfield(cm, 'entryid') && sum(sum(cm.matrix ~= 0 & cm.matrix ~= 1)) > 0
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = cm.entrydescname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = cm.entryidname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        [weights, si] = sort(cm.matrix(i,:), 'descend');
        entries = cm.entry(si);
        entrydescs = cm.entrydesc(si);
        entryids = cm.entryid(si);
        
        keep = weights ~= 0;
        lists.entries{i} = entries(keep);
        lists.entrydescs{i} = entrydescs(keep);
        lists.entryids{i} = entryids(keep);
        lists.weights{i} = weights(keep)';
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
elseif isfield(cm, 'entryid') && sum(sum(cm.matrix ~= 0 & cm.matrix ~= 1)) > 0
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = cm.entryidname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        [weights, si] = sort(cm.matrix(i,:), 'descend');
        entries = cm.entry(si);
        entryids = cm.entryid(si);
        
        keep = weights ~= 0;
        lists.entries{i} = entries(keep);
        lists.entryids{i} = entryids(keep);
        lists.weights{i} = weights(keep)';
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
elseif isfield(cm, 'entrydesc') && sum(sum(cm.matrix ~= 0 & cm.matrix ~= 1)) > 0
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = cm.entrydescname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        [weights, si] = sort(cm.matrix(i,:), 'descend');
        entries = cm.entry(si);
        entrydescs = cm.entrydesc(si);
        
        keep = weights ~= 0;
        lists.entries{i} = entries(keep);
        lists.entrydescs{i} = entrydescs(keep);
        lists.weights{i} = weights(keep)';
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
elseif isfield(cm, 'entrydesc') && isfield(cm, 'entryid') && sum(sum(cm.matrix ~= 0 & cm.matrix ~= 1)) == 0
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = cm.entrydescname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = cm.entryidname;
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        keep = cm.matrix(i,:) ~= 0;
        lists.entries{i} = cm.entry(keep);
        lists.entrydescs{i} = cm.entrydesc(keep);
        lists.entryids{i} = cm.entryid(keep);
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
elseif sum(sum(cm.matrix ~= 0 & cm.matrix ~= 1)) > 0
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.weights = cell(lists.numterms, 1);
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        [weights, si] = sort(cm.matrix(i,:), 'descend');
        entries = cm.entry(si);
        
        keep = weights ~= 0;
        lists.entries{i} = entries(keep);
        lists.weights{i} = weights(keep)';
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
elseif isfield(cm, 'entrydesc') && sum(sum(cm.matrix ~= 0 & cm.matrix ~= 1)) == 0
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.entrydescs = cell(lists.numterms, 1);
    lists.entrydescname = cm.entrydescname;
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        keep = cm.matrix(i,:) ~= 0;
        lists.entries{i} = cm.entry(keep);
        lists.entrydescs{i} = cm.entrydesc(keep);
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
elseif isfield(cm, 'entryid') && sum(sum(cm.matrix ~= 0 & cm.matrix ~= 1)) == 0
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.entryids = cell(lists.numterms, 1);
    lists.entryidname = cm.entryidname;
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        keep = cm.matrix(i,:) ~= 0;
        lists.entries{i} = cm.entry(keep);
        lists.entryids{i} = cm.entryid(keep);
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
else
    
    lists.entries = cell(lists.numterms, 1);
    lists.entryname = cm.entryname;
    lists.numentries = zeros([lists.numterms 1]);
    
    for i = 1:1:lists.numterms
        
        keep = cm.matrix(i,:) ~= 0;
        lists.entries{i} = cm.entry(keep);
        
        lists.numentries(i) = numel(lists.entries{i});
        
    end
    
end

discard = lists.numentries == 0;

lists.term(discard) = [];

if isfield(lists, 'termdesc')
    lists.termdesc(discard) = [];
end

if isfield(lists, 'termid')
    lists.termid(discard) = [];
end

lists.entries(discard) = [];

if isfield(lists, 'entrydescs')
    lists.entrydescs(discard) = [];
end

if isfield(lists, 'entryids')
    lists.entryids(discard) = [];
end

if isfield(lists, 'weights')
    lists.weights(discard) = [];
end

lists.numentries(discard) = [];

lists.numterms = numel(lists.term);


