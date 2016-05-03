


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
gene_atb.edges = edgesinit(319907, [], 'GeneID', [], 'BioGrid ID', [], 'GeneID', [], 'GeneID', [], 'BioGrid ID', [], 'GeneID', false, []);
detectionmethod = cell(gene_atb.edges.numedges, 1);
interactiontype = cell(gene_atb.edges.numedges, 1);
sourcetaxid = cell(gene_atb.edges.numedges, 1);
targettaxid = cell(gene_atb.edges.numedges, 1);

% read human data
fid = fopen('input/dataset_20160407_original.mitab', 'r');
currline = fgetl(fid);
for i = 1:1:gene_atb.edges.numedges
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    
    match = regexp(currcells{10}, 'taxid:(?<taxid>\d+)', 'names');
    sourcetaxid{i} = match.taxid;
    match = regexp(currcells{11}, 'taxid:(?<taxid>\d+)', 'names');
    targettaxid{i} = match.taxid;
    
    if strcmp(sourcetaxid{i}, '9606') && strcmp(targettaxid{i}, '9606')
        match = regexp(currcells{1}, 'entrez gene/locuslink:(?<geneid>\d+)', 'names');
        gene_atb.edges.source{i} = match.geneid;
        gene_atb.edges.sourceid(i) = str2double(gene_atb.edges.source{i});
        match = regexp(currcells{3}, 'biogrid:(?<biogridid>\d+)', 'names');
        gene_atb.edges.sourcedesc{i} = match.biogridid;

        match = regexp(currcells{2}, 'entrez gene/locuslink:(?<geneid>\d+)', 'names');
        gene_atb.edges.target{i} = match.geneid;
        gene_atb.edges.targetid(i) = str2double(gene_atb.edges.target{i});
        match = regexp(currcells{4}, 'biogrid:(?<biogridid>\d+)', 'names');
        gene_atb.edges.targetdesc{i} = match.biogridid;

        match = regexp(currcells{7}, '\((?<detectionmethod>.+)\)', 'names');
        detectionmethod{i} = lower(strtrim(match.detectionmethod));

        match = regexp(currcells{12}, '\((?<interactiontype>.+)\)', 'names');
        interactiontype{i} = lower(strtrim(match.interactiontype));
    end
end
fclose(fid);
discard = cellfun(@isempty, gene_atb.edges.source) | ismember(gene_atb.edges.source, {'-666' '' '-'}) | cellfun(@isempty, gene_atb.edges.target) | ismember(gene_atb.edges.target, {'-666' '' '-'});
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
detectionmethod(discard) = [];
interactiontype(discard) = [];
clear fid ans currcells currline match discard i sourcetaxid targettaxid;

% save unmodified data
save('output/gene_attribute_edges_unmodified.mat', '-struct', 'gene_atb');
save('output/extra_vars_unmodified.mat', '-mat', 'interactiontype', 'detectionmethod');

% discard interactions that are not confirmed protein-protein
[u_interactiontype, ~, ri] = unique(interactiontype);
c_interactiontype = zeros(size(u_interactiontype));
for i = 1:1:numel(u_interactiontype)
c_interactiontype(i) = sum(ri==i);
end
[c_interactiontype, si] = sort(c_interactiontype, 'descend');
u_interactiontype = u_interactiontype(si);
% discard = ~ismember(interactiontype, {'physical association' 'direct interaction' 'colocalization'});
discard = ~strcmp(interactiontype, 'direct interaction');
detectionmethod(discard) = [];
interactiontype(discard) = [];
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
clear interactiontype detectionmethod u_interactiontype c_interactiontype ri i si discard;

% convert to matrix format
gene_atb.cm = edges2cm(edges2symedges(edgesunique(gene_atb.edges)));
gene_atb = rmfield(gene_atb, 'edges');

% map gene identifiers to entrez gene symbols and discard rows corresponding to un-mapped symbols
% 1 unmapped out of 12497 (0%)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
disp([num2str(sum(discard)) ' unmapped out of ' num2str(numel(discard)) ' (' num2str(round(100*sum(discard)/numel(discard))) '%)']);
clear discard;

% map attribute identifiers to entrez gene symbols and discard rows corresponding to un-mapped symbols
% 1 unmapped out of 12497 (0%)
[gene_atb.cm.entry, gene_atb.cm.entryname, gene_atb.cm.entryid, gene_atb.cm.entryidname, discard] = genesymlookup([], gene_atb.cm.entryid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
disp([num2str(sum(discard)) ' unmapped out of ' num2str(numel(discard)) ' (' num2str(round(100*sum(discard)/numel(discard))) '%)']);
clear discard;

% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end

% merge cols corresponding to the same gene
if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'max');
end

% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);

% confirm matrix is symmetric
if gene_atb.cm.numterms ~= gene_atb.cm.numentries || ~all(strcmp(gene_atb.cm.term, gene_atb.cm.entry)) || ~all(all(gene_atb.cm.matrix == gene_atb.cm.matrix'))
    error('matrix is not symmetric');
else
    disp('matrix is symmetric');
end

% save cleaned data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);

% standardize
tftype = 'raw'; % doesn't matter for binary data
gene_atb.cm = cm_tfidf_standardization(gene_atb.cm, tftype);
figure(1); clf; hist(gene_atb.cm.matrix(gene_atb.cm.matrix~=0), 100);
support = 'unbounded';
fignum = 2;
gene_atb.cm = cm_ksdensity_standardization_sparsematrix(gene_atb.cm, support, fignum);
clear tftype support fignum;

% cluster and view standardized matrix
gene_atb.cm = cmcluster_nogram(gene_atb.cm, true);

% save standardized data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);
writecm('output/gene_attribute_matrix_standardized', gene_atb.cm);
gzip('output/gene_attribute_matrix_standardized.txt');
delete('output/gene_attribute_matrix_standardized.txt');
genes = gene_atb.cm.term;
atbs = gene_atb.cm.entry;

% threshold
gene_atb.cm.matrix = (gene_atb.cm.matrix > 0) - (gene_atb.cm.matrix < 0);

% save thresholded data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_thresholded.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);
writecm('output/gene_attribute_matrix_thresholded', gene_atb.cm);
gzip('output/gene_attribute_matrix_thresholded.txt');
delete('output/gene_attribute_matrix_thresholded.txt');
clear gene_atb;

% align cleaned matrix with standardized matrix and save
load('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
cm = conmatmap(cm, genes, atbs);
save('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
writecm('output/gene_attribute_matrix_cleaned', cm);
gzip('output/gene_attribute_matrix_cleaned.txt');
delete('output/gene_attribute_matrix_cleaned.txt');
clear cm genes atbs mappingfilespath;


