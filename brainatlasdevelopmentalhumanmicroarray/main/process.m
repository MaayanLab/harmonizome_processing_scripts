


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


% get row meta data
row = struct;
row.numterms = 17604;
row.term = repmat({'-666'}, row.numterms, 1);
row.termname = 'GeneSym';
row.termdesc = repmat({'-666'}, row.numterms, 1);
row.termdescname = 'EnsembleAcc';
row.termid = -666*ones([row.numterms 1]);
row.termidname = 'GeneID';

fid = fopen('input/HumanDevelopmentalTranscriptomeMicroarray/rows_metadata_fixed.txt', 'r');

currline = fgetl(fid);

for i = 1:1:row.numterms
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    row.term{i} = currcells{4};

    if numel(currcells{3}) > 0
        
        row.termdesc{i} = currcells{3};
        
    end
    
    if numel(currcells{5}) > 0
        
        row.termid(i) = str2double(currcells{5});
        
    end
    
end

fclose(fid);


% get col meta data
col = struct;
col.numterms = 492;
col.term = repmat({'-666'}, col.numterms, 1);
col.termname = 'Structure_Age_Gender_DonorID';
col.termdesc = repmat({'-666'}, col.numterms, 1);
col.termdescname = 'Age_Gender';
col.termid = -666*ones([col.numterms 1]);
col.termidname = 'StructureID';

fid = fopen('input/HumanDevelopmentalTranscriptomeMicroarray/columns_metadata.txt', 'r');

currline = fgetl(fid);

for i = 1:1:col.numterms
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    col.term{i} = [currcells{8} '_' currcells{4} '_' currcells{5} '_' currcells{2}];
    
    col.termdesc{i} = [currcells{4} '_' currcells{5}];
    
    col.termid(i) = str2double(currcells{6});
    
end

fclose(fid);


% initialize cm structure
gene_atb.cm = cminit(row.numterms, col.numterms, row.term, row.termname, row.termdesc, row.termdescname, row.termid, row.termidname, col.term, col.termname, col.termdesc, col.termdescname, col.termid, col.termidname, []);


% read data
gene_atb.cm.matrix = csvread('input/HumanDevelopmentalTranscriptomeMicroarray/expression_matrix.csv', 0, 1);


save('input/gene_tissuesample_brainatlas_humandevelopmentarray_unmodified.mat', '-struct', 'gene_atb');
%}


gene_atb = load('input/gene_tissuesample_brainatlas_humandevelopmentarray_unmodified.mat', '-mat', 'cm');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% view column distributions
figure(1);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% quantile normalize
gene_atb.cm.matrix = quantilenormalization(gene_atb.cm.matrix);


% view column distributions
figure(2);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end


% view column distributions
figure(3);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 4);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors); % 0.2671 by Age_Gender
% gene_atb.cm.matrix = robustzscore(gene_atb.cm.matrix, 2); % 0? by Age_Gender
% gene_atb.cm.matrix = zscore(gene_atb.cm.matrix, 0, 2); % 0.2551 by Age_Gender
daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm)); % 0.0490 by Age_Gender


% search for best threshold (lowest davies-bouldin index)
threshfrac = [(0.01:0.01:0.04) (0.05:0.05:0.95) (0.96:0.01:0.99)]';
type = 'tertiary';
method = 'matrix';
discardemptyvectors = false;

daviesbouldin = zeros([numel(threshfrac) 1]);

for i = 1:1:numel(threshfrac)
    
    tm = cmthresh(gene_atb.cm, threshfrac(i), type, method, discardemptyvectors);

    daviesbouldin(i) = cm2daviesbouldin(cmtranspose(tm));
    
end

figure(5); clf; plot(threshfrac, daviesbouldin, '-ok'); ylabel('davies-bouldin index'); xlabel('threshold');
%}


% threshold.  note, we could equivalently just use cmthresh to get the
% thresholded matrix, but using cmthresh_ks here in case we ever want the
% final normalized matrix
threshfrac = 0.15;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[gene_atb.cm, nlpm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm)); % 0.3278 by Age_Gender


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 10);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');






clearvars -except gene_atb nlpm;

gene_atb.cm = nlpm;
daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm)); % 0.2694 by Age_Gender
clear nlpm;


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 11);


% save result
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');


