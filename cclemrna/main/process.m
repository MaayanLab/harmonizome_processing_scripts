


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
gene_cl.cm = cminit(18988, 1037, [], 'GeneSym', [], [], [], [], [], 'CellLine', [], 'Tissue', [], [], []);

% read data into connectivity matrix
fid = fopen('input/CCLE_Expression_Entrez_2012-09-29.gct', 'r');

currline = fgetl(fid);
currline = fgetl(fid);
currline = fgetl(fid);

currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
currcells(1:2) = [];

for i = 1:1:gene_cl.cm.numentries
    
    subcells = strsplit(currcells{i}, '_');
    
    gene_cl.cm.entry{i} = subcells{1};
    gene_cl.cm.entrydesc{i} = lower(strjoin(subcells(2:end), ' '));
    
end

for i = 1:1:gene_cl.cm.numterms
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_cl.cm.term{i} = currcells{2};
    
    gene_cl.cm.matrix(i,:) = str2double(currcells(3:end));
    
end

fclose(fid);

% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
[gene_cl.cm.term, gene_cl.cm.termname, gene_cl.cm.termid, gene_cl.cm.termidname, discard] = genesymlookup(gene_cl.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_cl.cm = cmrowdiscard(gene_cl.cm, discard);

if sum(isnan(gene_cl.cm.matrix(:))) > 0
    % remove rows and columns with more than 5% values missing
    
    gene_cl.cm = cmnantrim_frac(gene_cl.cm, 0.95, Inf, 0.95, Inf, 'row');
    
end

if sum(isnan(gene_cl.cm.matrix(:))) > 0
    % impute remaining missing values
    
    gene_cl.cm = cmnanimpute(gene_cl.cm, 'row');
    
end

% view data distributions
figure(1); clf; hist(gene_cl.cm.matrix, 100);

if numel(unique(gene_cl.cm.term)) < gene_cl.cm.numterms
    % merge measurements corresponding to the same gene
    
    gene_cl.cm = cmrowmerge(gene_cl.cm, 'mean');
    
end

% quantile normalize
% gene_cl.cm.matrix = quantilenormalization(gene_cl.cm.matrix);

% cluster and view data
cgo = cm2clustergram(gene_cl.cm, 'none', 'all', 'cosine', 'average');

close force all;

gene_cl.cm = conmatmap(gene_cl.cm, cgo.RowLabels, cgo.ColumnLabels');

HeatMap(gene_cl.cm.matrix, 'Colormap', redbluecmap);

% view data distributions
figure(2); clf; hist(gene_cl.cm.matrix, 100);

% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_cl');

% make up down lists for qiaonan
gene_atb.zm = gene_atb.cm;
gene_atb.zm.matrix = zscore(gene_atb.zm.matrix, 0, 2);
atb_gene.lists = cm2lists(cmtranspose(gene_atb.zm));
for i = 1:1:atb_gene.lists.numterms
c1 = max(0, atb_gene.lists.weights{i}(1000));
c2 = min(0, atb_gene.lists.weights{i}(end-999));
discard = atb_gene.lists.weights{i} < c1 & atb_gene.lists.weights{i} > c2;
atb_gene.lists.weights{i}(discard) = [];
atb_gene.lists.entries{i}(discard) = [];
atb_gene.lists.entryids{i}(discard) = [];
atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
end
writelists(atb_gene.lists, 'input/cellline_gene_zscore_1000up1000dn', 'fuzzy');






% get standardized matrix

gene_atb = load('output/gene_attribute_matrix_imported.mat', '-mat', 'cm');

threshfrac = 0.05;
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);

threshfrac = 0.1;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);

% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);

% save result
if exist('output/gene_attribute_matrix_standardized.mat', 'file') == 0
    save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
else
    error('file already exists?!');
end


