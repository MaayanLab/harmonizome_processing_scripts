


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





% get sample meta data
numsamples = 100000;
sample_id = cell(numsamples, 1);
tissue_general = cell(numsamples, 1);
tissue_specific = cell(numsamples, 1);

fid = fopen('input/GTEx_Data_V6_Annotations_SampleAttributesDS_released20151019.txt', 'r');

currline = fgetl(fid);
currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
sample_id_hit = strcmp(currcells, 'SAMPID');
tissue_general_hit = strcmp(currcells, 'SMTS');
tissue_specific_hit = strcmp(currcells, 'SMTSD');

i = 0;
while ~feof(fid)
    i = i + 1;
    currline = fgetl(fid);
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    sample_id{i} = currcells{sample_id_hit};
    tissue_general{i} = currcells{tissue_general_hit};
    tissue_specific{i} = currcells{tissue_specific_hit};
end
fclose(fid);

if i < numsamples
    sample_id(i+1:end) = [];
    tissue_general(i+1:end) = [];
    tissue_specific(i+1:end) = [];
end
numsamples = numel(sample_id);

tissue_general = lower(tissue_general);
tissue_specific = lower(tissue_specific);


% get subject meta data
numsubjects = 10000;
subject_id = cell(numsubjects, 1);
gender = cell(numsubjects, 1);
age = cell(numsubjects, 1);

fid = fopen('input/GTEx_Data_V6_Annotations_SubjectPhenotypesDS_released20151019.txt', 'r');

currline = fgetl(fid);

i = 0;
while ~feof(fid)
    i = i + 1;
    currline = fgetl(fid);
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    subject_id{i} = currcells{1};
    if currcells{2} == '1'
        gender{i} = 'male';
    elseif currcells{2} == '2'
        gender{i} = 'female';
    end
    age{i} = currcells{3};
end
fclose(fid);

if i < numsubjects
    subject_id(i+1:end) = [];
    gender(i+1:end) = [];
    age(i+1:end) = [];
end
numsubjects = numel(subject_id);


% map subject meta data to samples
sample.numsamples = numsamples;
sample.sample_id = sample_id;
sample.tissue_general = tissue_general;
sample.tissue_specific = tissue_specific;
sample.subject_id = cell(sample.numsamples, 1);
for i = 1:1:sample.numsamples
    di = find(sample.sample_id{i} == '-');
    sample.subject_id{i} = sample.sample_id{i}(1:di(2)-1);
end
[o1, o2] = ismember(sample.subject_id, subject_id);
o2(o2 == 0) = [];
sample.gender = cell(sample.numsamples, 1);
sample.age = cell(sample.numsamples, 1);
sample.gender(o1) = gender(o2);
sample.age(o1) = age(o2);

discard = ~o1 | cellfun(@isempty, sample.tissue_general);
sample.sample_id(discard) = [];
sample.tissue_general(discard) = [];
sample.tissue_specific(discard) = [];
sample.subject_id(discard) = [];
sample.gender(discard) = [];
sample.age(discard) = [];
sample.numsamples = numel(sample.sample_id);

clearvars -except sample;


% get sample and subject counts
tissue_general = unique(sample.tissue_general);
numtissues = numel(tissue_general);
age = unique(sample.age);
numages = numel(age);
sampleids = cell(numtissues, numages);
numsamples = zeros(numtissues, numages);
subjectids = cell(numtissues, numages);
numsubjects = zeros(numtissues, numages);
for i = 1:1:numtissues
    for j = 1:1:numages
        hit = strcmp(sample.tissue_general, tissue_general{i}) & strcmp(sample.age, age{j});
        sampleids{i,j} = sample.sample_id(hit);
        numsamples(i,j) = numel(sampleids{i,j});
        subjectids{i,j} = unique(sample.subject_id(hit));
        numsubjects(i,j) = numel(subjectids{i,j});
    end
end


% % discard samples if no replicate exists for that tissue, age, subject triplet
% subject_id = unique(sample.subject_id);
% numsubjects = numel(subject_id);
% discard = false(sample.numsamples, 1);
% for i = 1:1:numtissues
%     for j = 1:1:numages
%         for k = 1:1:numsubjects
%             hit = strcmp(sample.tissue_general, tissue_general{i}) & strcmp(sample.age, age{j}) & strcmp(sample.subject_id, subject_id{k});
%             if sum(hit) == 1
%                 discard(hit) = true;
%             end
%         end
%     end
% end
% sample.sample_id(discard) = [];
% sample.tissue_general(discard) = [];
% sample.tissue_specific(discard) = [];
% sample.subject_id(discard) = [];
% sample.gender(discard) = [];
% sample.age(discard) = [];
% sample.numsamples = numel(sample.sample_id);
% 
% 
% % get sample and subject counts
% tissue_general = unique(sample.tissue_general);
% numtissues = numel(tissue_general);
% age = unique(sample.age);
% numages = numel(age);
% sampleids = cell(numtissues, numages);
% numsamples = zeros(numtissues, numages);
% subjectids = cell(numtissues, numages);
% numsubjects = zeros(numtissues, numages);
% for i = 1:1:numtissues
%     for j = 1:1:numages
%         hit = strcmp(sample.tissue_general, tissue_general{i}) & strcmp(sample.age, age{j});
%         sampleids{i,j} = sample.sample_id(hit);
%         numsamples(i,j) = numel(sampleids{i,j});
%         subjectids{i,j} = unique(sample.subject_id(hit));
%         numsubjects(i,j) = numel(subjectids{i,j});
%     end
% end


clearvars -except sample;


% initialize cm structure
fid = fopen('input/All_Tissue_Site_Details_Analysis_released20151019.combined.rpkm.gct', 'r');

currline = fgetl(fid);
currline = fgetl(fid);
currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
numrows = str2double(currcells{1});
numcols = str2double(currcells{2});

gene_atb.cm = cminit(numrows, numcols, [], 'GeneSym', [], 'Ensemble Acc', [], [], [], 'Sample ID', [], 'Subject ID_Tissue_Age_Gender', [], [], []);


% read column labels
currline = fgetl(fid);
currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
gene_atb.cm.entry = currcells(3:end)';
for i = 1:1:gene_atb.cm.numentries
    hit = strcmp(sample.sample_id, gene_atb.cm.entry{i});
    if sum(hit) == 1
        gene_atb.cm.entrydesc{i} = [sample.subject_id{hit} '_' sample.tissue_general{hit} '_' sample.age{hit} '_' sample.gender{hit}];
    else
        gene_atb.cm.entrydesc{i} = '-666';
    end
end


% read row labels
for i = 1:1:gene_atb.cm.numterms
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    gene_atb.cm.term{i} = currcells{2};
    gene_atb.cm.termdesc{i} = currcells{1};
end

fclose(fid);


% read matrix
gene_atb.cm.matrix = dlmread('input/All_Tissue_Site_Details_Analysis_released20151019.combined.rpkm.gct', '\t', 3, 2);

save('input/gene_tissuesample_gtex_unmodified20160103.mat', '-struct', 'gene_atb'); % too big to save


% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
missed = gene_atb.cm.term(discard);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);

save('input/gene_tissuesample_gtex_genesmapped20160103.mat', '-struct', 'gene_atb');


% remove columns missing metadata
discard = strcmp(gene_atb.cm.entrydesc, '-666');
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
clear discard;


% handle NaNs (there are none)
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % FALSE
    % remove rows and columns with more than 5% missing values
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.95, Inf, 0.95, Inf, 'column');
end

rowhasnansorzeros = sum(isnan(gene_atb.cm.matrix) | (gene_atb.cm.matrix==0), 2) > 0;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % FALSE
    % impute remaining missing values
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end


% % view column distributions
% figure(1);
% clf;
% subplot(2, 2, 1);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,:)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cellfun(@(x) x(1:find(x=='_',1,'first')-1), atb_gene.sm.termdesc, 'uniformoutput', false);
atb_gene.sm.termdescname = 'Tissue_Age_Gender';
atb_gene.sm.termdesc = cellfun(@(x) x(find(x=='_',1,'first')+1:end), atb_gene.sm.termdesc, 'uniformoutput', false);
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_genesmapped, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_genesmapped = interp1(fpr, tpr, fpr_x);
figure(11);
clf;
subplot(2,3,1);
plot(fpr_x, tpr_genesmapped, '-k');
title('genesmapped');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% find lower detection limit
detlim = zeros([1000 gene_atb.cm.numentries]);
posmin = zeros([1 gene_atb.cm.numentries]);
pos25 = zeros([1 gene_atb.cm.numentries]);
pos50 = zeros([1 gene_atb.cm.numentries]);
pos75 = zeros([1 gene_atb.cm.numentries]);
posmax = zeros([1 gene_atb.cm.numentries]);

for i = 1:1:gene_atb.cm.numentries
    x = gene_atb.cm.matrix(gene_atb.cm.matrix(:,i) > 0,i);
    x = sort(x);
    detlim(:,i) = x(1:1000);
    posmin(i) = x(1);
    pos25(i) = x(round(0.25*numel(x)));
    pos50(i) = x(round(0.50*numel(x)));
    pos75(i) = x(round(0.75*numel(x)));
    posmax(i) = x(end);
end

[posmin, si] = sort(posmin);
pos25 = pos25(si);
pos50 = pos50(si);
pos75 = pos75(si);
posmax = posmax(si);

% figure(2);
% clf;
% subplot(1,2,1);
% hist(log2(detlim), 100);
% subplot(1,2,2);
% plot(1:1:gene_atb.cm.numentries, log2(posmin), '-ok', 1:1:gene_atb.cm.numentries, log2(pos25), '-ok', 1:1:gene_atb.cm.numentries, log2(pos50), '-ok', 1:1:gene_atb.cm.numentries, log2(pos75), '-ok', 1:1:gene_atb.cm.numentries, log2(posmax), '-ok');


% handle zeros (expression level below detection)
gene_atb.cm.matrix(gene_atb.cm.matrix == 0) = NaN;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % remove rows and columns with more than 5% undetected values
    frac = (0.05:0.05:0.95)';
    for i = 1:1:numel(frac)
        gene_atb.cm = cmnantrim_frac(gene_atb.cm, frac(i), Inf, frac(i), Inf, 'column');
    end
end

rowhaszeros = sum(isnan(gene_atb.cm.matrix), 2) > 0;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % impute remaining undetected values
%     gene_atb.cm.matrix(isnan(gene_atb.cm.matrix)) = mean(posmin)/2;
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end
clear detlim frac pos25 pos50 pos75 posmax posmin si x;


% % view column distributions
% figure(3);
% clf;
% subplot(2, 2, 1);
% hist(log2(gene_atb.cm.matrix), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(log2(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(log2(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(log2(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cellfun(@(x) x(1:find(x=='_',1,'first')-1), atb_gene.sm.termdesc, 'uniformoutput', false);
atb_gene.sm.termdescname = 'Tissue_Age_Gender';
atb_gene.sm.termdesc = cellfun(@(x) x(find(x=='_',1,'first')+1:end), atb_gene.sm.termdesc, 'uniformoutput', false);
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_zeroshandled, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_zeroshandled = interp1(fpr, tpr, fpr_x);
figure(11);
subplot(2,3,2);
plot(fpr_x, tpr_zeroshandled, '-k');
title('zeroshandled');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;

% quantify extent to which gene similarity predicts ppis
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition = {'zeroshandled_angularcosine'};
clear X Y stats;
sm = cm2sm_angularpearson_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,2) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(2) = {'zeroshandled_angularpearson'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2pm(sm);
sm = pm2clrm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,3) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(3) = {'zeroshandled_angularcosine_clr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2rm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,4) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(4) = {'zeroshandled_angularcosine_rwr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm.matrix(1:sm.numterms+1:sm.numterms^2) = 0;
d = sum(sm.matrix, 2);
sm.matrix = bsxfun(@rdivide, bsxfun(@rdivide, sm.matrix, sqrt(d)), sqrt(d)');
clear d;
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,5) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(5) = {'zeroshandled_angularcosine_ngl'};
clear X Y stats;
sm = cm2sm_angularspearman_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,6) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(6) = {'zeroshandled_angularspearman'};
clear X Y stats;

% log2 transform
gene_atb.cm.matrix = log2(gene_atb.cm.matrix);


% % view column distributions
% figure(4);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cellfun(@(x) x(1:find(x=='_',1,'first')-1), atb_gene.sm.termdesc, 'uniformoutput', false);
atb_gene.sm.termdescname = 'Tissue_Age_Gender';
atb_gene.sm.termdesc = cellfun(@(x) x(find(x=='_',1,'first')+1:end), atb_gene.sm.termdesc, 'uniformoutput', false);
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_log2, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_log2 = interp1(fpr, tpr, fpr_x);
figure(11);
subplot(2,3,3);
plot(fpr_x, tpr_log2, '-k');
title('log2');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;

% quantify extent to which gene similarity predicts ppis
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,7) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(7) = {'zeroshandled_log2_angularcosine'};
clear X Y stats;
sm = cm2sm_angularpearson_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,8) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(8) = {'zeroshandled_log2_angularpearson'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2pm(sm);
sm = pm2clrm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,9) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(9) = {'zeroshandled_log2_angularcosine_clr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2rm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,10) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(10) = {'zeroshandled_log2_angularcosine_rwr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm.matrix(1:sm.numterms+1:sm.numterms^2) = 0;
d = sum(sm.matrix, 2);
sm.matrix = bsxfun(@rdivide, bsxfun(@rdivide, sm.matrix, sqrt(d)), sqrt(d)');
clear d;
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,11) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(11) = {'zeroshandled_log2_angularcosine_ngl'};
clear X Y stats;
sm = cm2sm_angularspearman_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,12) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(12) = {'zeroshandled_log2_angularspearman'};
clear X Y stats;


% quantile normalization
gene_atb.cm.matrix = quantilenormalization(gene_atb.cm.matrix);


% % view column distributions
% figure(5);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cellfun(@(x) x(1:find(x=='_',1,'first')-1), atb_gene.sm.termdesc, 'uniformoutput', false);
atb_gene.sm.termdescname = 'Tissue_Age_Gender';
atb_gene.sm.termdesc = cellfun(@(x) x(find(x=='_',1,'first')+1:end), atb_gene.sm.termdesc, 'uniformoutput', false);
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_qn, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_qn = interp1(fpr, tpr, fpr_x);
figure(11);
subplot(2,3,4);
plot(fpr_x, tpr_qn, '-k');
title('qn');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;

% quantify extent to which gene similarity predicts ppis
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,13) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(13) = {'zeroshandled_log2_qn_angularcosine'};
clear X Y stats;
sm = cm2sm_angularpearson_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,14) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(14) = {'zeroshandled_log2_qn_angularpearson'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2pm(sm);
sm = pm2clrm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,15) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(15) = {'zeroshandled_log2_qn_angularcosine_clr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2rm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,16) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(16) = {'zeroshandled_log2_qn_angularcosine_rwr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm.matrix(1:sm.numterms+1:sm.numterms^2) = 0;
d = sum(sm.matrix, 2);
sm.matrix = bsxfun(@rdivide, bsxfun(@rdivide, sm.matrix, sqrt(d)), sqrt(d)');
clear d;
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,17) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(17) = {'zeroshandled_log2_qn_angularcosine_ngl'};
clear X Y stats;
sm = cm2sm_angularspearman_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,18) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(18) = {'zeroshandled_log2_qn_angularspearman'};
clear X Y stats;


% merge measurements corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % TRUE
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end


% % view column distributions
% figure(6);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cellfun(@(x) x(1:find(x=='_',1,'first')-1), atb_gene.sm.termdesc, 'uniformoutput', false);
atb_gene.sm.termdescname = 'Tissue_Age_Gender';
atb_gene.sm.termdesc = cellfun(@(x) x(find(x=='_',1,'first')+1:end), atb_gene.sm.termdesc, 'uniformoutput', false);
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_mergegenes, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_mergegenes = interp1(fpr, tpr, fpr_x);
figure(11);
subplot(2,3,5);
plot(fpr_x, tpr_mergegenes, '-k');
title('mergegenes');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;

% quantify extent to which gene similarity predicts ppis
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,19) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(19) = {'zeroshandled_log2_qn_mergerows_angularcosine'};
clear X Y stats;
sm = cm2sm_angularpearson_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,20) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(20) = {'zeroshandled_log2_qn_mergerows_angularpearson'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2pm(sm);
sm = pm2clrm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,21) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(21) = {'zeroshandled_log2_qn_mergerows_angularcosine_clr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2rm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,22) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(22) = {'zeroshandled_log2_qn_mergerows_angularcosine_rwr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm.matrix(1:sm.numterms+1:sm.numterms^2) = 0;
d = sum(sm.matrix, 2);
sm.matrix = bsxfun(@rdivide, bsxfun(@rdivide, sm.matrix, sqrt(d)), sqrt(d)');
clear d;
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,23) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(23) = {'zeroshandled_log2_qn_mergerows_angularcosine_ngl'};
clear X Y stats;
sm = cm2sm_angularspearman_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,24) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(24) = {'zeroshandled_log2_qn_mergerows_angularspearman'};
clear X Y stats;


% cluster and view matrix
gene_atb.cm = cmcluster(gene_atb.cm, false);


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');




clearvars -except gene_atb;

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);


% standardize and threshold
threshfrac = 0.1;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[tm, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);


% cluster and view matrix
gene_atb.cm = cmcluster(gene_atb.cm, false);


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cellfun(@(x) x(1:find(x=='_',1,'first')-1), atb_gene.sm.termdesc, 'uniformoutput', false);
atb_gene.sm.termdescname = 'Tissue_Age_Gender';
atb_gene.sm.termdesc = cellfun(@(x) x(find(x=='_',1,'first')+1:end), atb_gene.sm.termdesc, 'uniformoutput', false);
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_std, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_std = interp1(fpr, tpr, fpr_x);
figure(12);
clf;
subplot(1,2,1);
plot(fpr_x, tpr_std, '-k');
title('std');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;

% quantify extent to which gene similarity predicts ppis
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,25) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(25) = {'zeroshandled_log2_qn_mergerows_standardized_angularcosine'};
clear X Y stats;
sm = cm2sm_angularpearson_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,26) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(26) = {'zeroshandled_log2_qn_mergerows_standardized_angularpearson'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2pm(sm);
sm = pm2clrm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,27) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(27) = {'zeroshandled_log2_qn_mergerows_standardized_angularcosine_clr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm = sm2rm(sm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,28) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(28) = {'zeroshandled_log2_qn_mergerows_standardized_angularcosine_rwr'};
clear X Y stats;
sm = cm2sm_angularcosine_nocluster(gene_atb.cm);
sm.matrix(1:sm.numterms+1:sm.numterms^2) = 0;
d = sum(sm.matrix, 2);
sm.matrix = bsxfun(@rdivide, bsxfun(@rdivide, sm.matrix, sqrt(d)), sqrt(d)');
clear d;
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,29) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(29) = {'zeroshandled_log2_qn_mergerows_standardized_angularcosine_ngl'};
clear X Y stats;
sm = cm2sm_angularspearman_nocluster(gene_atb.cm);
gene_gene = load('C:/Users/Andrew/Dropbox/Ioutput/gene_attribute_matrix_imported.mat', '-mat');
gene_gene.cm.matrix = full(gene_gene.cm.matrix);
commongenes = intersect(sm.term, gene_gene.cm.term);
sm = simmatmap(sm, commongenes);
X = sm.matrix(tril(true(size(sm.matrix)),-1));
clear sm;
gene_gene.cm = conmatmap(gene_gene.cm, commongenes, commongenes);
Y = gene_gene.cm.matrix(tril(true(size(gene_gene.cm.matrix)),-1));
clear gene_gene commongenes;
stats = classifierperformance(Y, X, [], 5);
summary(:,30) = [stats.idx stats.val stats.f1s stats.mcc stats.fdr stats.mcr stats.tpr stats.fpr stats.fdrmin stats.idxfdrmin stats.valfdrmin stats.auc]';
condition(30) = {'zeroshandled_log2_qn_mergerows_standardized_angularspearman'};
clear X Y stats;


% save standardized matrix
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');


% cluster thresholded matrix
gene_atb.cm = conmatmap(tm, gene_atb.cm.term, gene_atb.cm.entry);
HeatMap(gene_atb.cm.matrix, 'colormap', redbluecmap);
clear tm;


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cellfun(@(x) x(1:find(x=='_',1,'first')-1), atb_gene.sm.termdesc, 'uniformoutput', false);
atb_gene.sm.termdescname = 'Tissue_Age_Gender';
atb_gene.sm.termdesc = cellfun(@(x) x(find(x=='_',1,'first')+1:end), atb_gene.sm.termdesc, 'uniformoutput', false);
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_thr, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_thr = interp1(fpr, tpr, fpr_x);
figure(12);
subplot(1,2,2);
plot(fpr_x, tpr_thr, '-k');
title('thr');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save thresholded matrix
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');











%{
tissue = cell(gene_atb.cm.numentries, 1);
for i = 1:1:gene_atb.cm.numentries
ui = find(gene_atb.cm.entrydesc{i} == '_');
tissue{i} = gene_atb.cm.entrydesc{i}(ui(1)+1:ui(2)-1);
end
[ut, ui, ri] = unique(tissue);
nt = zeros(numel(ut), 1);
for i = 1:1:numel(ut)
    nt(i) = sum(ri==i);
end
isbrain = strcmp(tissue, 'brain');
temp = cmcoldiscard(gene_atb.cm, ~isbrain);
[~, X] = pca(temp.matrix');

age = cell(temp.numentries, 1);
for i = 1:1:temp.numentries
ui = find(temp.entrydesc{i} == '_');
age{i} = temp.entrydesc{i}(ui(2)+1:ui(3)-1);
end
ua = unique(age);
figure(100);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,1), X(hit,2), 'or');
hit = ismember(age, ua(5:6));
plot(X(hit,1), X(hit,2), 'sb');
hold off;

gender = cell(temp.numentries, 1);
for i = 1:1:temp.numentries
ui = find(temp.entrydesc{i} == '_');
gender{i} = temp.entrydesc{i}(ui(3)+1:end);
end
ug = unique(gender);
figure(200);
clf;
hold on;
hit = strcmp(gender, ug{1});
plot(X(hit,1), X(hit,2), 'or');
hit = strcmp(gender, ug{2});
plot(X(hit,1), X(hit,2), 'sb');
hold off;

tissuesubtype = cell(temp.numentries, 1);
[o1, o2] = ismember(temp.entry, sample.sample_id);
o2(o2 == 0) = [];
tissuesubtype(o1) = sample.tissue_specific(o2);
[uts, ui, ri] = unique(tissuesubtype);
nts = zeros(numel(uts), 1);
for i = 1:1:numel(uts)
    nts(i) = sum(ri==i);
end
figure(300);
clf;
hold on;
hit = strcmp(tissuesubtype, uts{1});
plot(X(hit,1), X(hit,2), 'or');
hit = strcmp(tissuesubtype, uts{2});
plot(X(hit,1), X(hit,2), 'og');
hit = strcmp(tissuesubtype, uts{3});
plot(X(hit,1), X(hit,2), 'ob');
hit = strcmp(tissuesubtype, uts{4});
plot(X(hit,1), X(hit,2), 'oc');
hit = strcmp(tissuesubtype, uts{5});
plot(X(hit,1), X(hit,2), 'om');
hit = strcmp(tissuesubtype, uts{6});
plot(X(hit,1), X(hit,2), 'ok');
hit = strcmp(tissuesubtype, uts{7});
plot(X(hit,1), X(hit,2), 'oy');
hit = strcmp(tissuesubtype, uts{8});
plot(X(hit,1), X(hit,2), 'xr');
hit = strcmp(tissuesubtype, uts{9});
plot(X(hit,1), X(hit,2), 'xg');
hit = strcmp(tissuesubtype, uts{10});
plot(X(hit,1), X(hit,2), 'xb');
hit = strcmp(tissuesubtype, uts{11});
plot(X(hit,1), X(hit,2), 'xc');
hit = strcmp(tissuesubtype, uts{12});
plot(X(hit,1), X(hit,2), 'xm');
hit = strcmp(tissuesubtype, uts{13});
plot(X(hit,1), X(hit,2), 'xk');
hold off;
legend(uts, 'location', 'southeast');

figure(300);
clf;
hold on;
hit = strcmp(tissuesubtype, uts{1});
plot(X(hit,1), X(hit,2), 'or');
hit = strcmp(tissuesubtype, uts{2});
plot(X(hit,1), X(hit,2), 'og');
hit = strcmp(tissuesubtype, uts{3});
plot(X(hit,1), X(hit,2), 'ob');
hit = strcmp(tissuesubtype, uts{4});
plot(X(hit,1), X(hit,2), 'oc');
hit = strcmp(tissuesubtype, uts{5});
plot(X(hit,1), X(hit,2), 'om');
hit = strcmp(tissuesubtype, uts{6});
plot(X(hit,1), X(hit,2), 'ok');
hit = strcmp(tissuesubtype, uts{7});
plot(X(hit,1), X(hit,2), 'oy');
hit = strcmp(tissuesubtype, uts{8});
plot(X(hit,1), X(hit,2), 'xr');
hit = strcmp(tissuesubtype, uts{9});
plot(X(hit,1), X(hit,2), 'xg');
hit = strcmp(tissuesubtype, uts{10});
plot(X(hit,1), X(hit,2), 'xb');
hit = strcmp(tissuesubtype, uts{11});
plot(X(hit,1), X(hit,2), 'xc');
hit = strcmp(tissuesubtype, uts{12});
plot(X(hit,1), X(hit,2), 'xm');
hit = strcmp(tissuesubtype, uts{13});
plot(X(hit,1), X(hit,2), 'xk');
hold off;
legend(uts, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-2');
figure(301);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,1), X(hit,2), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,1), X(hit,2), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,1), X(hit,2), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,1), X(hit,2), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,1), X(hit,2), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,1), X(hit,2), 'ok');
hold off;
legend(ua, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-2');
figure(302);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,1), X(hit,2), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,1), X(hit,2), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,1), X(hit,2), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-2');

figure(400);
clf;
hold on;
hit = strcmp(tissuesubtype, uts{1});
plot(X(hit,1), X(hit,3), 'or');
hit = strcmp(tissuesubtype, uts{2});
plot(X(hit,1), X(hit,3), 'og');
hit = strcmp(tissuesubtype, uts{3});
plot(X(hit,1), X(hit,3), 'ob');
hit = strcmp(tissuesubtype, uts{4});
plot(X(hit,1), X(hit,3), 'oc');
hit = strcmp(tissuesubtype, uts{5});
plot(X(hit,1), X(hit,3), 'om');
hit = strcmp(tissuesubtype, uts{6});
plot(X(hit,1), X(hit,3), 'ok');
hit = strcmp(tissuesubtype, uts{7});
plot(X(hit,1), X(hit,3), 'oy');
hit = strcmp(tissuesubtype, uts{8});
plot(X(hit,1), X(hit,3), 'xr');
hit = strcmp(tissuesubtype, uts{9});
plot(X(hit,1), X(hit,3), 'xg');
hit = strcmp(tissuesubtype, uts{10});
plot(X(hit,1), X(hit,3), 'xb');
hit = strcmp(tissuesubtype, uts{11});
plot(X(hit,1), X(hit,3), 'xc');
hit = strcmp(tissuesubtype, uts{12});
plot(X(hit,1), X(hit,3), 'xm');
hit = strcmp(tissuesubtype, uts{13});
plot(X(hit,1), X(hit,3), 'xk');
hold off;
legend(uts, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-3');
figure(401);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,1), X(hit,3), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,1), X(hit,3), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,1), X(hit,3), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,1), X(hit,3), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,1), X(hit,3), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,1), X(hit,3), 'ok');
hold off;
legend(ua, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-3');
figure(402);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,1), X(hit,3), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,1), X(hit,3), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,1), X(hit,3), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-3');

figure(500);
clf;
hold on;
hit = strcmp(tissuesubtype, uts{1});
plot(X(hit,2), X(hit,3), 'or');
hit = strcmp(tissuesubtype, uts{2});
plot(X(hit,2), X(hit,3), 'og');
hit = strcmp(tissuesubtype, uts{3});
plot(X(hit,2), X(hit,3), 'ob');
hit = strcmp(tissuesubtype, uts{4});
plot(X(hit,2), X(hit,3), 'oc');
hit = strcmp(tissuesubtype, uts{5});
plot(X(hit,2), X(hit,3), 'om');
hit = strcmp(tissuesubtype, uts{6});
plot(X(hit,2), X(hit,3), 'ok');
hit = strcmp(tissuesubtype, uts{7});
plot(X(hit,2), X(hit,3), 'oy');
hit = strcmp(tissuesubtype, uts{8});
plot(X(hit,2), X(hit,3), 'xr');
hit = strcmp(tissuesubtype, uts{9});
plot(X(hit,2), X(hit,3), 'xg');
hit = strcmp(tissuesubtype, uts{10});
plot(X(hit,2), X(hit,3), 'xb');
hit = strcmp(tissuesubtype, uts{11});
plot(X(hit,2), X(hit,3), 'xc');
hit = strcmp(tissuesubtype, uts{12});
plot(X(hit,2), X(hit,3), 'xm');
hit = strcmp(tissuesubtype, uts{13});
plot(X(hit,2), X(hit,3), 'xk');
hold off;
legend(uts, 'location', 'southeast');
xlabel('PC-2');
ylabel('PC-3');
figure(501);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,2), X(hit,3), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,2), X(hit,3), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,2), X(hit,3), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,2), X(hit,3), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,2), X(hit,3), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,2), X(hit,3), 'ok');
hold off;
legend(ua, 'location', 'southeast');
xlabel('PC-2');
ylabel('PC-3');
figure(502);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,2), X(hit,3), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,2), X(hit,3), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,2), X(hit,3), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-2');
ylabel('PC-3');

figure(600);
clf;
hold on;
hit = strcmp(tissuesubtype, uts{1});
plot(X(hit,1), X(hit,4), 'or');
hit = strcmp(tissuesubtype, uts{2});
plot(X(hit,1), X(hit,4), 'og');
hit = strcmp(tissuesubtype, uts{3});
plot(X(hit,1), X(hit,4), 'ob');
hit = strcmp(tissuesubtype, uts{4});
plot(X(hit,1), X(hit,4), 'oc');
hit = strcmp(tissuesubtype, uts{5});
plot(X(hit,1), X(hit,4), 'om');
hit = strcmp(tissuesubtype, uts{6});
plot(X(hit,1), X(hit,4), 'ok');
hit = strcmp(tissuesubtype, uts{7});
plot(X(hit,1), X(hit,4), 'oy');
hit = strcmp(tissuesubtype, uts{8});
plot(X(hit,1), X(hit,4), 'xr');
hit = strcmp(tissuesubtype, uts{9});
plot(X(hit,1), X(hit,4), 'xg');
hit = strcmp(tissuesubtype, uts{10});
plot(X(hit,1), X(hit,4), 'xb');
hit = strcmp(tissuesubtype, uts{11});
plot(X(hit,1), X(hit,4), 'xc');
hit = strcmp(tissuesubtype, uts{12});
plot(X(hit,1), X(hit,4), 'xm');
hit = strcmp(tissuesubtype, uts{13});
plot(X(hit,1), X(hit,4), 'xk');
hold off;
legend(uts, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-4');
figure(601);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,1), X(hit,4), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,1), X(hit,4), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,1), X(hit,4), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,1), X(hit,4), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,1), X(hit,4), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,1), X(hit,4), 'ok');
hold off;
legend(ua, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-4');
figure(602);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,1), X(hit,4), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,1), X(hit,4), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,1), X(hit,4), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-1');
ylabel('PC-4');

figure(700);
clf;
hold on;
hit = strcmp(tissuesubtype, uts{1});
plot(X(hit,3), X(hit,4), 'or');
hit = strcmp(tissuesubtype, uts{2});
plot(X(hit,3), X(hit,4), 'og');
hit = strcmp(tissuesubtype, uts{3});
plot(X(hit,3), X(hit,4), 'ob');
hit = strcmp(tissuesubtype, uts{4});
plot(X(hit,3), X(hit,4), 'oc');
hit = strcmp(tissuesubtype, uts{5});
plot(X(hit,3), X(hit,4), 'om');
hit = strcmp(tissuesubtype, uts{6});
plot(X(hit,3), X(hit,4), 'ok');
hit = strcmp(tissuesubtype, uts{7});
plot(X(hit,3), X(hit,4), 'oy');
hit = strcmp(tissuesubtype, uts{8});
plot(X(hit,3), X(hit,4), 'xr');
hit = strcmp(tissuesubtype, uts{9});
plot(X(hit,3), X(hit,4), 'xg');
hit = strcmp(tissuesubtype, uts{10});
plot(X(hit,3), X(hit,4), 'xb');
hit = strcmp(tissuesubtype, uts{11});
plot(X(hit,3), X(hit,4), 'xc');
hit = strcmp(tissuesubtype, uts{12});
plot(X(hit,3), X(hit,4), 'xm');
hit = strcmp(tissuesubtype, uts{13});
plot(X(hit,3), X(hit,4), 'xk');
hold off;
legend(uts, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-4');
figure(701);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,3), X(hit,4), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,3), X(hit,4), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,3), X(hit,4), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,3), X(hit,4), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,3), X(hit,4), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,3), X(hit,4), 'ok');
hold off;
legend(ua, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-4');
figure(702);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,4), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,4), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,4), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-4');

figure(801);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,1), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-1');
figure(802);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,2), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,2), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,2), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-2');
figure(803);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,3), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,3), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,3), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-3');
figure(804);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,4), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,4), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,4), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-4');
figure(805);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,5), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,5), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,5), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-5');
figure(806);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,6), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,6), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,6), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-6');
figure(807);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,7), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,7), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,7), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-7');
figure(808);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,8), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,8), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,8), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-8');
figure(809);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,9), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,9), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,9), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-9');
figure(810);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), X(hit,10), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), X(hit,10), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), X(hit,10), 'ob');
hold off;
legend({'young';'middle';'old'}, 'location', 'southeast');
xlabel('PC-3');
ylabel('PC-10');

figure(911);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,1), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,1), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,1), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-1');
ylabel('age rank');
ylim([0 4]);
figure(912);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,1), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,1), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,1), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,1), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,1), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,1), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-1');
ylabel('age rank');
ylim([0 7]);
figure(913);
clf;
boxplot(X(:,1), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-1');


figure(921);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,2), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,2), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,2), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-2');
ylabel('age rank');
ylim([0 4]);
figure(922);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,2), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,2), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,2), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,2), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,2), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,2), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-2');
ylabel('age rank');
ylim([0 7]);
figure(923);
clf;
boxplot(X(:,2), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-2');


figure(931);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,3), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,3), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,3), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-3');
ylabel('age rank');
ylim([0 4]);
figure(932);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,3), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,3), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,3), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,3), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,3), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,3), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-3');
ylabel('age rank');
ylim([0 7]);
figure(933);
clf;
boxplot(X(:,3), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-1');


figure(941);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,4), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,4), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,4), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-4');
ylabel('age rank');
ylim([0 4]);
figure(942);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,4), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,4), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,4), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,4), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,4), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,4), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-4');
ylabel('age rank');
ylim([0 7]);
figure(943);
clf;
boxplot(X(:,4), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-4');


figure(951);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,5), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,5), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,5), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-5');
ylabel('age rank');
ylim([0 4]);
figure(952);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,5), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,5), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,5), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,5), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,5), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,5), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-5');
ylabel('age rank');
ylim([0 7]);
figure(953);
clf;
boxplot(X(:,5), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-5');


figure(961);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,6), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,6), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,6), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-6');
ylabel('age rank');
ylim([0 4]);
figure(962);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,6), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,6), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,6), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,6), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,6), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,1), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-6');
ylabel('age rank');
ylim([0 7]);
figure(963);
clf;
boxplot(X(:,6), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-6');


figure(971);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,7), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,7), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,7), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-7');
ylabel('age rank');
ylim([0 4]);
figure(972);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,7), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,7), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,7), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,7), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,7), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,7), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-7');
ylabel('age rank');
ylim([0 7]);
figure(973);
clf;
boxplot(X(:,7), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-7');


figure(981);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,8), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,8), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,8), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-8');
ylabel('age rank');
ylim([0 4]);
figure(982);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,8), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,8), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,8), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,8), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,8), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,8), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-8');
ylabel('age rank');
ylim([0 7]);
figure(983);
clf;
boxplot(X(:,8), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-8');


figure(991);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,9), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,9), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,9), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-9');
ylabel('age rank');
ylim([0 4]);
figure(992);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,9), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,9), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,9), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,9), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,9), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,9), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-9');
ylabel('age rank');
ylim([0 7]);
figure(993);
clf;
boxplot(X(:,9), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-9');


figure(9101);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,10), 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,10), 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,10), 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-10');
ylabel('age rank');
ylim([0 4]);
figure(9102);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,10), 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,10), 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,10), 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,10), 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,10), 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,10), 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-10');
ylabel('age rank');
ylim([0 7]);
figure(9103);
clf;
boxplot(X(:,10), age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-10');



v = zeros(size(X,2), 1);
v([3 7]) = 1;
v = v/sqrt(sum(v.^2));
figure(1001);
clf;
hold on;
hit = ismember(age, ua(1:2));
plot(X(hit,:)*v, 1*ones(sum(hit),1), 'or');
hit = ismember(age, ua(3:4));
plot(X(hit,:)*v, 2*ones(sum(hit),1), 'ok');
hit = ismember(age, ua(5:6));
plot(X(hit,:)*v, 3*ones(sum(hit),1), 'ob');
hold off;
xlabel('PC-3,7');
ylabel('age rank');
ylim([0 4]);
figure(1002);
clf;
hold on;
hit = strcmp(age, ua{1});
plot(X(hit,:)*v, 1*ones(sum(hit),1), 'or');
hit = strcmp(age, ua{2});
plot(X(hit,:)*v, 2*ones(sum(hit),1), 'og');
hit = strcmp(age, ua{3});
plot(X(hit,:)*v, 3*ones(sum(hit),1), 'om');
hit = strcmp(age, ua{4});
plot(X(hit,:)*v, 4*ones(sum(hit),1), 'oc');
hit = strcmp(age, ua{5});
plot(X(hit,:)*v, 5*ones(sum(hit),1), 'ob');
hit = strcmp(age, ua{6});
plot(X(hit,:)*v, 6*ones(sum(hit),1), 'ok');
hold off;
xlabel('PC-3,7');
ylabel('age rank');
ylim([0 7]);
figure(1003);
clf;
boxplot(X*v, age, 'orientation', 'horizontal', 'grouporder', ua);
xlabel('PC-3,7');
%}



%{
iscerebellum = strcmp(tissuesubtype, 'brain - cerebellum');
temp = cmcoldiscard(temp, ~iscerebellum);
[~, X] = pca(temp.matrix');

age = cell(temp.numentries, 1);
for i = 1:1:temp.numentries
ui = find(temp.entrydesc{i} == '_');
age{i} = temp.entrydesc{i}(ui(2)+1:ui(3)-1);
end
ua = unique(age);
figure(400);
clf;
hold on;
hit = ismember(age, ua(1:3));
scatter3(X(hit,1), X(hit,2), X(hit,3), 25, 'r');
hit = ismember(age, ua(5:6));
scatter3(X(hit,1), X(hit,2), X(hit,3), 25, 'b');
hold off;

gender = cell(temp.numentries, 1);
for i = 1:1:temp.numentries
ui = find(temp.entrydesc{i} == '_');
gender{i} = temp.entrydesc{i}(ui(3)+1:end);
end
ug = unique(gender);
figure(500);
clf;
hold on;
hit = strcmp(gender, ug{1});
plot(X(hit,1), X(hit,2), 'or');
hit = strcmp(gender, ug{2});
plot(X(hit,1), X(hit,2), 'sb');
hold off;
%}



%{
tissue = cell(gene_atb.cm.numentries, 1);
for i = 1:1:gene_atb.cm.numentries
ui = find(gene_atb.cm.entrydesc{i} == '_');
tissue{i} = gene_atb.cm.entrydesc{i}(ui(1)+1:ui(2)-1);
end

[ut, ui, ri] = unique(tissue);
nt = zeros(numel(ut), 1);
for i = 1:1:numel(ut)
    nt(i) = sum(ri==i);
end

keep = nt > 300 & nt < 400;
kt = ut(keep);
temp = cmcoldiscard(gene_atb.cm, ~ismember(tissue, kt));
[~, X] = pca(temp.matrix');

tissue = cell(temp.numentries, 1);
for i = 1:1:temp.numentries
ui = find(temp.entrydesc{i} == '_');
tissue{i} = temp.entrydesc{i}(ui(1)+1:ui(2)-1);
end
ut = unique(tissue);

figure(300);
clf;
hold on;
hit = strcmp(tissue, ut{1});
plot(X(hit,1), X(hit,2), 'or');
hit = strcmp(tissue, ut{2});
plot(X(hit,1), X(hit,2), 'og');
hit = strcmp(tissue, ut{3});
plot(X(hit,1), X(hit,2), 'ob');
hit = strcmp(tissue, ut{4});
plot(X(hit,1), X(hit,2), 'oc');
hit = strcmp(tissue, ut{5});
plot(X(hit,1), X(hit,2), 'om');
hold off;
%}


