


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
gene_atb.edges = edgesinit(1041228, [], 'GeneSym', [], 'EvidenceType', [], 'GeneID', [], 'Chemical', [], 'PubchemCID', [], [], false, []);
taxid = -666*ones([gene_atb.edges.numedges 1]);


% read data
fid = fopen('input/CTD_chem_gene_ixns.tsv', 'r');

for i = 1:1:28
    currline = fgetl(fid);
end

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb.edges.source{i} = currcells{4};
    
    gene_atb.edges.sourceid(i) = str2double(currcells{5});
    
    gene_atb.edges.target{i} = currcells{1};
    
    gene_atb.edges.targetdesc{i} = currcells{2};
    
    taxid(i) = str2double(currcells{8});
    
    gene_atb.edges.sourcedesc{i} = currcells{6};
    
end

fclose(fid);


% discard interactions that are not from human, mouse, or rat data
keep = cellfun(@numel, gene_atb.edges.source) > 0 & cellfun(@numel, gene_atb.edges.target) > 0 & ismember(taxid, [10090 9606 10116]);
gene_atb.edges = edgesdiscard(gene_atb.edges, ~keep);


save('input/gene_chemical_ctd_filtered_20150415.mat', '-struct', 'gene_atb');
%}


gene_atb = load('input/gene_chemical_ctd_filtered_20150415.mat', '-mat', 'edges');


% discard interactions that are not gene or protein type (alternative is
% mRNA, i.e. chemical regulates expression)
keep = ismember(gene_atb.edges.sourcedesc, {'protein' 'gene'});
gene_atb.edges = edgesdiscard(gene_atb.edges, ~keep);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.edges = rmfield(gene_atb.edges, {'sourcedesc' 'sourcedescname'});
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


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;


% convert to full matrix
gene_atb.cm.matrix = full(gene_atb.cm.matrix);


% remove rows and cols with too many connections
threshfrac = 1/2;
gene_atb.cm = cmtrim_frac(gene_atb.cm, 0, threshfrac, 0, threshfrac, 'column');
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


