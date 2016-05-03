function outcm = cmtranspose(incm)

outcm.term = incm.entry;
outcm.termname = incm.entryname;

if isfield(incm, 'entrydesc')
    outcm.termdesc = incm.entrydesc;
    outcm.termdescname = incm.entrydescname;
end

if isfield(incm, 'entryid')
    outcm.termid = incm.entryid;
    outcm.termidname = incm.entryidname;
end

outcm.numterms = incm.numentries;

outcm.entry = incm.term;
outcm.entryname = incm.termname;

if isfield(incm, 'termdesc')
    outcm.entrydesc = incm.termdesc;
    outcm.entrydescname = incm.termdescname;
end

if isfield(incm, 'termid')
    outcm.entryid = incm.termid;
    outcm.entryidname = incm.termidname;
end

outcm.numentries = incm.numterms;

outcm.matrix = incm.matrix';


