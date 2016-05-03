


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

%{

% initialize edges structure
gene_atb.edges = edgesinit(10000000, [], 'GeneSym', [], 'EvidenceType', [], 'GeneID', [], 'Disease', [], 'Mesh or Omim ID', [], [], true, []);
gene_atb.edges.sourcedesc = repmat({'inferred'}, gene_atb.edges.numedges, 1);


% read data
fid = fopen('input/CTD_genes_diseases_filtered2.tsv', 'r');

i = 0;

while ~feof(fid)
    
    i = i + 1;
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if currcells{1} > 0
        gene_atb.edges.source{i} = currcells{1};
    end
    
    if currcells{2} > 0
        gene_atb.edges.sourceid(i) = str2double(currcells{2});
    end
    
    if currcells{3} > 0
        gene_atb.edges.target{i} = currcells{3};
    end
    
    if currcells{4} > 0
        gene_atb.edges.targetdesc{i} = currcells{4};
    end
    
    if numel(currcells{5}) > 0
        gene_atb.edges.sourcedesc{i} = currcells{5};
    end
    
    if currcells{7} > 0
        gene_atb.edges.weight(i) = str2double(currcells{7});
    end
    
end

fclose(fid);

if i < gene_atb.edges.numedges
    discard = false([gene_atb.edges.numedges 1]);
    discard(i+1:end) = true;
    gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
end

discard = strcmp(gene_atb.edges.source, '-666') | gene_atb.edges.sourceid == -666 | strcmp(gene_atb.edges.target, '-666') | strcmp(gene_atb.edges.targetdesc, '-666');
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% normalize weights
gene_atb.edges.weight(~strcmp(gene_atb.edges.sourcedesc, 'inferred')) = 1000;
gene_atb.edges.weight = gene_atb.edges.weight/1000;


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


edges = gene_atb.edges;
clear gene_atb;
save('input/gene_disease_ctd_filtered_20150415.mat', '-struct', 'edges');
%}


gene_atb.edges = load('input/gene_disease_ctd_filtered_20150415.mat', '-mat');


% convert to matrix format (attribute table)
gene_atb.edges = rmfield(gene_atb.edges, {'sourcedesc' 'sourcedescname'});  %  'weight'});
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'both', true, true, mappingfilespath);
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


