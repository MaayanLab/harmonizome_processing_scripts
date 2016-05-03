function cm = cmnantrim_frac(cm, MinFracPerRow, MaxFracPerRow, MinFracPerCol, MaxFracPerCol, FirstDim)

switch FirstDim
    
    case 'row'
        
        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 2) < MinFracPerRow*cm.numentries | sum(hitmat, 2) > MaxFracPerRow*cm.numentries;
        cm = cmrowdiscard(cm, discard);

        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 1) < MinFracPerCol*cm.numterms | sum(hitmat, 1) > MaxFracPerCol*cm.numterms;
        cm = cmcoldiscard(cm, discard);

        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 2) < MinFracPerRow*cm.numentries | sum(hitmat, 2) > MaxFracPerRow*cm.numentries;
        cm = cmrowdiscard(cm, discard);

        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 1) < MinFracPerCol*cm.numterms | sum(hitmat, 1) > MaxFracPerCol*cm.numterms;
        cm = cmcoldiscard(cm, discard);
        
    case 'column'

        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 1) < MinFracPerCol*cm.numterms | sum(hitmat, 1) > MaxFracPerCol*cm.numterms;
        cm = cmcoldiscard(cm, discard);
        
        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 2) < MinFracPerRow*cm.numentries | sum(hitmat, 2) > MaxFracPerRow*cm.numentries;
        cm = cmrowdiscard(cm, discard);

        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 1) < MinFracPerCol*cm.numterms | sum(hitmat, 1) > MaxFracPerCol*cm.numterms;
        cm = cmcoldiscard(cm, discard);
        
        hitmat = ~isnan(cm.matrix);
        discard = sum(hitmat, 2) < MinFracPerRow*cm.numentries | sum(hitmat, 2) > MaxFracPerRow*cm.numentries;
        cm = cmrowdiscard(cm, discard);
        
    otherwise
        
        disp('FirstDim must be row or column.');
        
end


