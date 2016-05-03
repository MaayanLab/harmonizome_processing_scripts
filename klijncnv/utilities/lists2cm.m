function cm = lists2cm(lists)

cm = struct;


cm.term = lists.term;
cm.termname = lists.termname;

if isfield(lists, 'termdesc')
    cm.termdesc = lists.termdesc;
    cm.termdescname = lists.termdescname;
end

if isfield(lists, 'termid')
    cm.termid = lists.termid;
    cm.termidname = lists.termidname;
end

cm.numterms = lists.numterms;


[cm.entry, ui, ~] = unique(vertcat(lists.entries{:}));
cm.entryname = lists.entryname;

if isfield(lists, 'entrydescs')
    entrydescs = vertcat(lists.entrydescs{:});
    cm.entrydesc = entrydescs(ui);
    cm.entrydescname = lists.entrydescname;
end

if isfield(lists, 'entryids')
    entryids = vertcat(lists.entryids{:});
    cm.entryid = entryids(ui);
    cm.entryidname = lists.entryidname;
end

cm.numentries = numel(cm.entry);


cm.matrix = zeros([cm.numterms cm.numentries]);

if isfield(lists, 'weights')
    
    for i = 1:1:lists.numterms
        
        [o1, o2] = ismember(cm.entry, lists.entries{i});
        o2(o2 == 0) = [];
        cm.matrix(i,o1) = lists.weights{i}(o2);
        
    end
    
else
    
    for i = 1:1:lists.numterms
        
        cm.matrix(i,:) = ismember(cm.entry, lists.entries{i});
        
    end
    
end


