


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

% initialize connectivity matrix structure
gene_atb.cm = cminit(5711, 216, [], 'GeneSym', [], [], [], [], [], 'CellLine', [], 'Tissue', [], [], []);

% read data into connectivity matrix
fid = fopen('input/dataset_20160407_original.gct', 'r');
currline = fgetl(fid);
currline = fgetl(fid);
currline = fgetl(fid);
currcells = strsplitbyadr(currline, '\t');
currcells(1:2) = [];
for i = 1:1:gene_atb.cm.numentries
    subcells = strsplitbyadr(currcells{i}, '_');
    gene_atb.cm.entry{i} = subcells{1};
    gene_atb.cm.entrydesc{i} = lower(strjoin(subcells(2:end), ' '));
end
for i = 1:1:gene_atb.cm.numterms
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    gene_atb.cm.term{i} = currcells{2};
    gene_atb.cm.matrix(i,:) = str2double(currcells(3:end));
end
fclose(fid);
clear ans currcells currline i subcells fid;

% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
% 226 unmapped out of 5711 (4%)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
disp([num2str(sum(discard)) ' unmapped out of ' num2str(numel(discard)) ' (' num2str(round(100*sum(discard)/numel(discard))) '%)']);
clear discard;

% remove rows and columns with missing values
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 1, 1, 1, 1, 'row');
end

% impute remaining missing values
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % FALSE
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'row');
end


% merge measurements corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % TRUE
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end

% remove suspect rows that have identical values to at least one other row
sm = cm2sm_cosine_nocluster(gene_atb.cm);
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectgenes = sm.term(numidentical > 0);
discard = ismember(gene_atb.cm.term, suspectgenes);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
clear sm numidentical suspectgenes discard;

% save cleaned data
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');

% standardize and threshold
threshfrac = 0.1;
type = 'tertiary';
discardemptyvectors = false;
method = 'rows';
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
method = 'matrix';
[cm, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
cm.matrix = sparse(cm.matrix);
clear threshfrac type discardemptyvectors method;

% cluster and view standardized matrix
gene_atb.cm = cmcluster_nogram(gene_atb.cm, true);

% save standardized data
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
writecm('output/gene_attribute_matrix_standardized', gene_atb.cm);
gzip('output/gene_attribute_matrix_standardized.txt');
delete('output/gene_attribute_matrix_standardized.txt');
genes = gene_atb.cm.term;
atbs = gene_atb.cm.entry;
clear gene_atb;

% align thresholded matrix with standardized matrix and save
cm = conmatmap(cm, genes, atbs);
save('output/gene_attribute_matrix_thresholded.mat', '-mat', 'cm');
cm.matrix = full(cm.matrix);
writecm('output/gene_attribute_matrix_thresholded', cm);
gzip('output/gene_attribute_matrix_thresholded.txt');
delete('output/gene_attribute_matrix_thresholded.txt');
clear cm;

% align cleaned matrix with standardized matrix and save
load('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
cm = conmatmap(cm, genes, atbs);
save('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
writecm('output/gene_attribute_matrix_cleaned', cm);
gzip('output/gene_attribute_matrix_cleaned.txt');
delete('output/gene_attribute_matrix_cleaned.txt');
clear cm genes atbs mappingfilespath;


