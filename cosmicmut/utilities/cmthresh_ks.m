function [cm, nlpm] = cmthresh_ks(cm, FDR, type, method, discardemptyvectors, support)

if exist('support', 'var') == 0 || isempty(support)
    support = 'unbounded';
end

switch type
    
    case 'binary'
        
        pmatrix = 1 - kscumulativeprobability(cm.matrix, method, support);
        
        nlpm = cm;
        nlpm.matrix = -log10(pmatrix);
        
        is_significant = pmatrix <= FDR;
        
        cm.matrix = double(is_significant);
    
    case 'tertiary'
        
        pmatrix = kscumulativeprobability(cm.matrix, method, support);
        
        smatrix = sign(pmatrix - 0.5);
        
        pmatrix(pmatrix > 0.5) = 1 - pmatrix(pmatrix > 0.5);
        
        pmatrix = 2*pmatrix;
        
        nlpm = cm;
        nlpm.matrix = -log10(pmatrix).*smatrix;
        
        is_significant = pmatrix <= FDR;
        
        cm.matrix = is_significant.*smatrix;
        
    otherwise
        
        warning('Invalid type. Specify ''binary'' or ''tertiary''.');
        
end

if discardemptyvectors
    
    cm = cmtrim(cm, 1, Inf, 1, Inf);
    
end


