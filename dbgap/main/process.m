


clear all;

% get path to mapping files and add utilities to search path
currentpath = pwd;
parentpath = currentpath(1:find(currentpath=='\', 1, 'last')-1);
mappingfilespath = [parentpath '\mapping'];
utilitiespath = [parentpath '\utilities'];
searchpaths = strsplit(path, ';')';
if ~ismember(utilitiespath, searchpaths)
    addpath(utilitiespath, '-begin');
end
clear currentpath parentpath utilitiespath searchpaths;



% initialize edges structure
gene_atb.edges = edgesinit(44124, [], 'GeneSym', [], 'rsID_Chromosome_Position', [], 'GeneID', [], 'Trait', [], 'MESH Category', [], [], true, []);
context = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene2sym = repmat({'-666'}, gene_atb.edges.numedges, 1);


% read data
fid = fopen('input/RESULTS_all.TAB', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    context{i} = currcells{8};
    gene2sym{i} = currcells{19}; 
    
    gene_atb.edges.source{i} = currcells{15};
    
    gene_atb.edges.sourcedesc{i} = [currcells{3} '_' currcells{5} '_' currcells{6}];
    
    gene_atb.edges.sourceid(i) = str2double(currcells{14});
    
    gene_atb.edges.target{i} = currcells{2};
    
    gene_atb.edges.targetdesc{i} = currcells{11};
    
    gene_atb.edges.weight(i) = str2double(currcells{4});  % p-value
    
end

fclose(fid);


% discard intergenic SNPs
discard = strcmp(context, 'Intergenic'); % same as discard = ~strcmp(gene_atb.edges.source, gene2sym);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% convert to matrix format (attribute table)
gene_atb.edges.weight = -log10(gene_atb.edges.weight);
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);


% save result
if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');
if issparse(gene_atb.cm.matrix)
    gene_atb.cm.matrix = full(gene_atb.cm.matrix);
end
savepath = 'output/';
mkdir(savepath);
writecm([savepath 'gene_attribute_matrix_cleaned'], gene_atb.cm);
gzip([savepath 'gene_attribute_matrix_cleaned.txt']);
delete([savepath 'gene_attribute_matrix_cleaned.txt']);






% get standardized matrix
clearvars -except savepath gene_atb;

% threshfrac = 0.1;
% type = 'binary';
% method = 'matrix';
% discardemptyvectors = true;
% [~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
gene_atb.cm = cm_standardize_ignorezeros(gene_atb.cm);

% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);

% save result
if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
if issparse(gene_atb.cm.matrix)
    gene_atb.cm.matrix = full(gene_atb.cm.matrix);
end
writecm([savepath 'gene_attribute_matrix_standardized'], gene_atb.cm);
gzip([savepath 'gene_attribute_matrix_standardized.txt']);
delete([savepath 'gene_attribute_matrix_standardized.txt']);


