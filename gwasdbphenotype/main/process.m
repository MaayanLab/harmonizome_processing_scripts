


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
gene_atb.edges = edgesinit(500000, [], 'GeneSym', [], 'rsID_Chromosome_Position_Context', [], [], [], 'Phenotype', [], 'HPOID', [], [], true, []);
% gene_atb.edges = edgesinit(500000, [], 'GeneSym', [], 'rsID_Chromosome_Position_Context', [], [], [], 'Disease', [], 'DOID', [], [], true, []);
context = repmat({'-666'}, gene_atb.edges.numedges, 1);

% read data
fid = fopen('input/gwasdb_20140908_snp_trait', 'r');

currline = fgetl(fid);

i = 0;

while ~feof(fid)
    
    currline = fgetl(fid);
    
    currcells = strsplitbyadr(currline, '\t');

    target = strsplitbyadr(currcells{17}, '|'); % HPO
%     target = strsplitbyadr(currcells{19}, '|'); % DO
    
    targetdesc = strsplitbyadr(currcells{16}, '|'); % HPOID
%     targetdesc = strsplitbyadr(currcells{18}, '|'); % DOID
    
    for j = 1:1:numel(target)
        
        i = i + 1;
        
        gene_atb.edges.source{i} = currcells{29};
        
        gene_atb.edges.sourcedesc{i} = [currcells{3} '_' currcells{1} '_' currcells{2} '_' currcells{30}];
        
        gene_atb.edges.weight(i) = str2double(currcells{8});  % p-value
        
        context{i} = currcells{30};
        
        gene_atb.edges.target{i} = target{j};

        gene_atb.edges.targetdesc{i} = targetdesc{j};
        
    end
    
end

fclose(fid);

if i < gene_atb.edges.numedges
    discard = false([gene_atb.edges.numedges 1]);
    discard(i+1:end) = true;
    gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
    context(discard) = [];
end


% discard incomplete information
discard = strcmp(context, 'NA') | strcmp(gene_atb.edges.source, 'NA') | strcmp(gene_atb.edges.target, 'NA') | strcmp(gene_atb.edges.targetdesc, 'NA') | isnan(gene_atb.edges.weight);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map attributes to ontology and discard un-mapped attributes
atb_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');
id_alt = load('output/gene_attribute_matrix_imported.mat', '-mat');
% atb_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');
% id_alt = load('output/gene_attribute_matrix_imported.mat', '-mat');

%%%%%%
gene_atb.edges.targetdesc = cellfun(@(x) strrep(x, 'HPOID:', 'HP:'), gene_atb.edges.targetdesc, 'UniformOutput', false);
%%%%%%

[o1, o2] = ismember(gene_atb.edges.targetdesc, id_alt.edges.targetdesc);
o2(o2 == 0) = [];

[o3, o4] = ismember(gene_atb.edges.targetdesc, atb_atb.cm.entrydesc);
o4(o4 == 0) = [];

gene_atb.edges.target(o1) = id_alt.edges.source(o2);
gene_atb.edges.targetdesc(o1) = id_alt.edges.sourcedesc(o2);
gene_atb.edges.target(o3) = atb_atb.cm.entry(o4);
gene_atb.edges.targetdesc(o3) = atb_atb.cm.entrydesc(o4);

discard = ~o1 & ~o3;
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% remove duplicate edges
gene_atb.edges.weight = -log10(gene_atb.edges.weight);
gene_atb.edges = edgesunique(gene_atb.edges);
maxrealval = -log10(realmin);
gene_atb.edges.weight(gene_atb.edges.weight > maxrealval) = maxrealval;


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% make sure all appropriate connections to parent attributes are made
% gene_atb.em = cm2em_qryavg(gene_atb.cm, atb_atb.cm, false, []);
% discard = ismember(gene_atb.em.entry, gene_atb.cm.entry);
% gene_atb.em = cmcoldiscard(gene_atb.em, discard);
% if gene_atb.em.numentries > 0 && gene_atb.em.numterms == gene_atb.cm.numterms && all(strcmp(gene_atb.em.term, gene_atb.cm.term))
%     gene_atb.cm = cmtranspose(cmvertcat_prealigned(cmtranspose(gene_atb.cm), cmtranspose(gene_atb.em)));
%     gene_atb = rmfield(gene_atb, 'em');
% else
%     error('problem!');
% end
% gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);
gene_atb.em = cm2em_qryavg(gene_atb.cm, atb_atb.cm, false, []);
gene_atb.cm = conmatmap(gene_atb.cm, gene_atb.em.term, gene_atb.em.entry);
gene_atb.em.matrix(gene_atb.cm.matrix~=0) = gene_atb.cm.matrix(gene_atb.cm.matrix~=0);
gene_atb.cm = gene_atb.em;
gene_atb = rmfield(gene_atb, 'em');
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
gene_atb.cm = cmcluster_nogram(gene_atb.cm, true);


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
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
gene_atb.cm = cmcluster_nogram(gene_atb.cm, true);

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


