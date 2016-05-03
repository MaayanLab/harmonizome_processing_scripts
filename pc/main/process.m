


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
gene_atb.edges = edgesinit(4901543, [], 'GeneSym', [], [], [], 'GeneID', [], 'GeneSym', [], [], [], 'GeneID', false, []);


% read data
fid = fopen('input/Pathway Commons.7.All.BINARY_SIF.hgnc.sif', 'r');

discard = true([gene_atb.edges.numedges 1]);

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if ismember(currcells{2}, {'controls-state-change-of'; 'controls-transport-of'; 'controls-phosphorylation-of'; 'interacts-with'; 'in-complex-with'})
        
        discard(i) = false;
    
        gene_atb.edges.source{i} = currcells{1};
    
        gene_atb.edges.target{i} = currcells{3};
        
    end
    
end

fclose(fid);

gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map attribute identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.target, gene_atb.edges.targetname, gene_atb.edges.targetid, gene_atb.edges.targetidname, discard] = genesymlookup(gene_atb.edges.target, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% get edges from source to target and from target to source
gene_atb.edges = edges2symedges(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


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


