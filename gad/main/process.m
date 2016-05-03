


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
gene_atb.edges = edgesinit(178795, [], 'GeneSym', [], 'Context', [], [], [], 'Disease', [], 'DiseaseClass', [], [], true, []);
evidence = repmat({'-666'}, gene_atb.edges.numedges, 1);
discard = false([gene_atb.edges.numedges 1]);

% read data
fid = fopen('input/associations.tsv', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.edges.numedges
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells) < 89
        currcells(end+1:89) = {''};
        discard(i) = true;
    end
    
    gene_atb.edges.source{i} = currcells{6};
    
    gene_atb.edges.sourcedesc{i} = currcells{14};
    
    gene_atb.edges.target{i} = strtrim(strrep(lower(currcells{3}), '|', '; '));
    
    gene_atb.edges.targetdesc{i} = lower(currcells{4});
    
    gene_atb.edges.weight(i) = str2double(currcells{8});  % p-value
    
    evidence{i} = currcells{39};
    
end

fclose(fid);


% discard lines improperly delimited, lines missing gene symbol, lines missing disease, intergenic variants, or variants where evidence is for NO association 
discard = discard | strcmp(gene_atb.edges.sourcedesc, 'intergenic') | strcmp(evidence, 'N') | strcmp(gene_atb.edges.source, '') | strcmp(gene_atb.edges.target, '');
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% convert to matrix format (attribute table)
gene_atb.edges = rmfield(gene_atb.edges, {'weight' 'sourcedesc' 'sourcedescname'});
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


