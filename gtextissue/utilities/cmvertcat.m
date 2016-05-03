function outcm = cmvertcat(incm1, incm2)

outentry = union(incm1.entry, incm2.entry);

incm1 = conmatmap(incm1, [], outentry);

incm2 = conmatmap(incm2, [], outentry);

outcm = struct;

outcm.term = [incm1.term; incm2.term];
outcm.termname = incm1.termname;

if isfield(incm1, 'termdesc') && isfield(incm2, 'termdesc')
    outcm.termdesc = [incm1.termdesc; incm2.termdesc];
    outcm.termdescname = incm1.termdescname;
end

if isfield(incm1, 'termid') && isfield(incm2, 'termid')
    outcm.termid = [incm1.termid; incm2.termid];
    outcm.termidname = incm1.termidname;
end

outcm.numterms = numel(outcm.term);

outcm.entry = incm1.entry;
outcm.entryname = incm1.entryname;

if isfield(incm1, 'entrydesc')
    outcm.entrydesc = incm1.entrydesc;
    outcm.entrydescname = incm1.entrydescname;
elseif isfield(incm2, 'entrydesc')
    outcm.entrydesc = incm2.entrydesc;
    outcm.entrydescname = incm2.entrydescname;
end

if isfield(incm1, 'entryid')
    outcm.entryid = incm1.entryid;
    outcm.entryidname = incm1.entryidname;
elseif isfield(incm2, 'entryid')
    outcm.entryid = incm2.entryid;
    outcm.entryidname = incm2.entryidname;
end

outcm.numentries = numel(outcm.entry);

outcm.matrix = [incm1.matrix; incm2.matrix];


