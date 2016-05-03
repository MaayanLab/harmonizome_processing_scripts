function cm = cmtrim(cm, MinNumPerRow, MaxNumPerRow, MinNumPerCol, MaxNumPerCol)

numrows = cm.numterms;
numcols = cm.numentries;

discard = sum(cm.matrix ~= 0, 2) < MinNumPerRow | sum(cm.matrix ~= 0, 2) > MaxNumPerRow;
cm = cmrowdiscard(cm, discard);

if MaxNumPerCol == numrows - 1
    MaxNumPerCol = cm.numterms - 1;
end
numrows = cm.numterms;

discard = sum(cm.matrix ~= 0, 1) < MinNumPerCol | sum(cm.matrix ~= 0, 1) > MaxNumPerCol;
cm = cmcoldiscard(cm, discard);

if MaxNumPerRow == numcols - 1
    MaxNumPerRow = cm.numentries - 1;
end
numcols = cm.numentries;

discard = sum(cm.matrix ~= 0, 2) < MinNumPerRow | sum(cm.matrix ~= 0, 2) > MaxNumPerRow;
cm = cmrowdiscard(cm, discard);

if MaxNumPerCol == numrows - 1
    MaxNumPerCol = cm.numterms - 1;
end
numrows = cm.numterms;

discard = sum(cm.matrix ~= 0, 1) < MinNumPerCol | sum(cm.matrix ~= 0, 1) > MaxNumPerCol;
cm = cmcoldiscard(cm, discard);


