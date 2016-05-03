function writecm(filename, cm)

rowlabelname1 = cm.termname;
rowlabel1 = cm.term;

if ~isfield(cm, 'termdesc')
    rowlabelname2 = [];
    rowlabel2 = [];
else
    rowlabelname2 = cm.termdescname;
    rowlabel2 = cm.termdesc;
end

if ~isfield(cm, 'termid')
    rowlabelname3 = [];
    rowlabel3 = [];
else
    rowlabelname3 = cm.termidname;
    rowlabel3 = cm.termid;
end

collabelname1 = cm.entryname;
collabel1 = cm.entry;

if ~isfield(cm, 'entrydesc')
    collabelname2 = [];
    collabel2 = [];
else
    collabelname2 = cm.entrydescname;
    collabel2 = cm.entrydesc;
end

if ~isfield(cm, 'entryid')
    collabelname3 = [];
    collabel3 = [];
else
    collabelname3 = cm.entryidname;
    collabel3 = cm.entryid;
end

writematrix(filename, cm.matrix, collabelname1, collabel1, collabelname2, collabel2, collabelname3, collabel3, rowlabelname1, rowlabel1, rowlabelname2, rowlabel2, rowlabelname3, rowlabel3);


