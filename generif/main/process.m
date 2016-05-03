


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
gene_atb.edges = edgesinit(2551337, [], 'GeneID', [], [], [], 'GeneID', [], 'Word', [], [], [], [], false, []);


% read data
fid = fopen('input/generifs_basic_20150416_human_filtered_uniqueedges_selectedwords3.txt', 'r');

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t');
    
    gene_atb.edges.source{i} = currcells{1};
    
    gene_atb.edges.sourceid(i) = str2double(currcells{1});
    
    gene_atb.edges.target{i} = currcells{2};
    
end

fclose(fid);

edges = gene_atb.edges;
clear gene_atb;
save('input/gene_word_generif_20150415.mat', '-struct', 'edges');
%}


gene_atb.edges = load('input/gene_word_generif_20150415.mat', '-mat');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup([], gene_atb.edges.sourceid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
% close force all;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;


% convert to full matrix
gene_atb.cm.matrix = full(gene_atb.cm.matrix);


% remove rows and cols with too many connections
% threshfrac = 1/2;
% gene_atb.cm = cmtrim_frac(gene_atb.cm, 0, threshfrac, 0, threshfrac, 'column');
% gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);
% ALREADY BELOW CUTOFF


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


