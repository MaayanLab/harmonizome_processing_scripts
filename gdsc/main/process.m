


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
gene_cl.cm = cminit(13321, 624, [], 'GeneSym', [], [], [], [], [], 'CellLine', [], [], [], [], []);

% read data into connectivity matrix
fid = fopen('input/en_input_w5.csv', 'r');

currline = fgetl(fid);

gene_cl.cm.entry = strsplit(currline, ',', 'CollapseDelimiters', false)';

for i = 1:1:gene_cl.cm.numterms
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, ',', 'CollapseDelimiters', false);
    
    gene_cl.cm.term{i} = currcells{1};
    
    gene_cl.cm.matrix(i,:) = str2double(currcells(2:end));
    
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

if (max(gene_cl.cm.matrix(:)) - min(gene_cl.cm.matrix(:))) > 30
    % log2 transformation

    gene_cl.cm.matrix = log2(gene_cl.cm.matrix);

end

if numel(unique(gene_cl.cm.term)) < gene_cl.cm.numterms
    % merge measurements corresponding to the same gene
    
    gene_cl.cm = cmrowmerge(gene_cl.cm, 'mean');
    
end

% view data distributions
figure(2); clf; hist(gene_cl.cm.matrix, 100);

% quantile normalize
% gene_cl.cm.matrix = quantilenormalization(gene_cl.cm.matrix);

% view data distributions
% figure(3); clf; hist(gene_cl.cm.matrix, 100);

% cluster and view data
cgo = cm2clustergram(gene_cl.cm, 'none', 'all', 'cosine', 'average');

close force all;

gene_cl.cm = conmatmap(gene_cl.cm, cgo.RowLabels, cgo.ColumnLabels');

HeatMap(gene_cl.cm.matrix, 'Colormap', redbluecmap);

% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_cl');






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


