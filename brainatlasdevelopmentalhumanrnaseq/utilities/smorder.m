function sm = smorder(sm, idxorder)

sm.term = sm.term(idxorder);

if isfield(sm, 'termdesc')
    
    sm.termdesc = sm.termdesc(idxorder);
    
end

if isfield(sm, 'termid')
    
    sm.termid = sm.termid(idxorder);
    
end

sm.matrix = sm.matrix(idxorder,idxorder);

% sm.vector = sm.matrix(tril(ones(size(sm.matrix)), -1) == 1);

% sm.scaledmatrix = squarematrixscale(sm.matrix);


