function cm = cm_tfidf_standardization(cm, tftype)

switch tftype
    case 'raw'
        coltfweights = 1;
        rowtfweights = 1;
    case 'normalized'
        coltfweights = 1./max(abs(cm.matrix), [], 1);
        rowtfweights = 1./max(abs(cm.matrix), [], 2);
    otherwise
        error('invalid tf type');
end

colhits = sum(cm.matrix~=0, 1);
colidfweights = log(1 + max(colhits)./colhits);
rowhits = sum(cm.matrix~=0, 2);
rowidfweights = log(1 + max(rowhits)./rowhits);

cm.matrix = bsxfun(@times, bsxfun(@times, cm.matrix, coltfweights.*colidfweights), rowtfweights.*rowidfweights);


