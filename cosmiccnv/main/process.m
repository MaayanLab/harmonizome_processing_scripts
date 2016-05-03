


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




% loci_cl.edges = edgesinit(58316, [], 'ChromosomeLocation', [], 'CNVclass', [], [], [], 'CellLine', [], 'Tissue', [], 'CellLineID', true, []);
loci_cl.edges = edgesinit(58316, [], 'ChromosomeLocation', [], [], [], [], [], 'CellLine', [], 'Tissue', [], 'CellLineID', true, []);

chr = cell(loci_cl.edges.numedges, 1);
start = zeros([loci_cl.edges.numedges 1]);
finish = zeros([loci_cl.edges.numedges 1]);

fid = fopen('input/CosmicCLPCompleteCNA.tsv', 'r');

currline = fgetl(fid);

tic;
for i = 1:1:loci_cl.edges.numedges
    
    currline = upper(fgetl(fid));
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    loci_cl.edges.source{i} = currcells{13};
%     loci_cl.edges.sourcedesc{i} = currcells{11};
    loci_cl.edges.target{i} = currcells{8};
    loci_cl.edges.targetdesc{i} = currcells{4};
    loci_cl.edges.targetid(i) = str2double(currcells{2});
%     loci_cl.edges.weight(i) = str2double(currcells{9}) - 2;
    loci_cl.edges.weight(i) = 2*strcmp(currcells{11}, 'GAIN')-1;
    
    subcells = strsplit(currcells{13}, ':');
    
    chr{i} = ['chr' subcells{1}];
    
    subcells = strsplit(subcells{2}, '..');
    
    start(i) = str2double(subcells{1});
    
    finish(i) = str2double(subcells{2});
    
end

fclose(fid);

loci_cl.edges.targetdesc = lower(loci_cl.edges.targetdesc);
hit = strcmp(loci_cl.edges.targetdesc, 'ns');
loci_cl.edges.targetdesc(hit) = {'-666'};

load('C:/Users/Andrew/Dropbox/GeneInfo/gene_loci_GRCh37_hg19.mat', '-mat', 'GeneSym', 'GeneID', 'ChrNum', 'StartPos', 'EndPos');

[uloci, ui, ri] = unique(loci_cl.edges.source);
uchr = chr(ui);
ustart = start(ui);
ufinish = finish(ui);
unum = numel(uchr);

hit = strcmp(uchr, 'chr23');
uchr(hit) = {'chrX'};
hit = strcmp(uchr, 'chr24');
uchr(hit) = {'chrY'};

loci_gene.lists = listsinit(unum, uloci, 'ChromosomeLocation', [], [], [], [], [], 'GeneSym', [], [], [], 'GeneID', false, []);

for i = 1:1:loci_gene.lists.numterms
    
%     hit = strcmp(ChrNum, uchr{i}) & StartPos >= ustart(i) & EndPos <= ufinish(i);
    hit = strcmp(ChrNum, uchr{i}) & ((StartPos >= ustart(i) & StartPos < ufinish(i)) | (EndPos > ustart(i) & EndPos <= ufinish(i)));
    
    [loci_gene.lists.entries{i}, ui, ri] = unique(GeneSym(hit));
    loci_gene.lists.entryids{i} = GeneID(ui);
    loci_gene.lists.numentries(i) = numel(loci_gene.lists.entries{i});
    
end

discard = loci_gene.lists.numentries == 0;
loci_gene.lists = listsdiscard(loci_gene.lists, discard);

loci_gene.edges = lists2edges(loci_gene.lists);

cl_loci.edges = sourcetargetswap(loci_cl.edges);
cl_loci.lists = edges2lists(cl_loci.edges);

cl_gene.lists = listsinit(cl_loci.lists.numterms, cl_loci.lists.term, cl_loci.lists.termname, cl_loci.lists.termdesc, cl_loci.lists.termdescname, cl_loci.lists.termid, cl_loci.lists.termidname, [], 'GeneSym', [], [], [], [], true, []);

for i = 1:1:cl_gene.lists.numterms
    
    [o1, o2] = ismember(loci_gene.edges.source, cl_loci.lists.entries{i});
    o2(o2 == 0) = [];
    weight = zeros([loci_gene.edges.numedges 1]);
    weight(o1) = cl_loci.lists.weights{i}(o2);
    
    cl_gene.lists.entries{i} = loci_gene.edges.target(o1);
    cl_gene.lists.weights{i} = weight(o1);
    cl_gene.lists.numentries(i) = numel(cl_gene.lists.entries{i});
    
end

discard = cl_gene.lists.numentries == 0;
cl_gene.lists = listsdiscard(cl_gene.lists, discard);

cl_gene.cm = lists2cm(cl_gene.lists);

gene_cl.cm = cmtranspose(cl_gene.cm);

% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
[gene_cl.cm.term, gene_cl.cm.termname, gene_cl.cm.termid, gene_cl.cm.termidname, discard] = genesymlookup(gene_cl.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_cl.cm = cmrowdiscard(gene_cl.cm, discard);

% trim
gene_cl.cm = cmtrim(gene_cl.cm, 1, 100, 1, Inf);  % ***

% cluster and view data
cgo = cm2clustergram(gene_cl.cm, 'none', 'all', 'cosine', 'average');

close force all;

gene_cl.cm = conmatmap(gene_cl.cm, cgo.RowLabels, cgo.ColumnLabels');

HeatMap(gene_cl.cm.matrix, 'Colormap', redbluecmap);

% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_cl');






% get standardized matrix (already thresholded?)

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





%{

load('input/cnv_filenames.mat', '-mat', 'CellLineID');
CellLineID = str2double(CellLineID);

gene_cl = load('output/gene_attribute_matrix_imported.mat', '-mat');

[o1, o2] = ismember(CellLineID, gene_cl.cm.entryid);
o2(o2 == 0) = [];

CellLine = repmat({'-666'}, numel(CellLineID), 1);
CellLine(o1) = gene_cl.cm.entry(o2);

Tissue = repmat({'-666'}, numel(CellLineID), 1);
Tissue(o1) = gene_cl.cm.entrydesc(o2);

CellLine(~o1) = [];
Tissue(~o1) = [];
CellLineID(~o1) = [];

clear gene_cl;

cl_gene.uplists = listsinit(numel(CellLine), CellLine, 'CellLine', Tissue, 'Tissue', CellLineID, 'CellLineID', [], 'GeneSym', [], [], [], [], false, []);
cl_gene.dnlists = cl_gene.uplists;

tic;
for i = 1:1:cl_gene.uplists.numterms
    
    cl_gene.uplists.entries{i} = cell(100000, 1);
    cl_gene.dnlists.entries{i} = cell(100000, 1);
    
    fid = fopen(['input/' num2str(cl_gene.uplists.termid(i)) '.txt'], 'r');
    
    currline = fgetl(fid);
    upj = 0;
    dnj = 0;
    
    while ~feof(fid)
        
        currline = fgetl(fid);
        
        currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
        
        hidx = find(currcells{2} == '_', 1, 'first');
        
        if ~isempty(hidx)
            
            currcells{2} = currcells{2}(1:hidx-1);
            
        end
        
        if strcmpi(currcells{7}, 'gain')
            
            upj = upj + 1;
            cl_gene.uplists.entries{i}{upj} = upper(currcells{2});
            
        elseif strcmpi(currcells{7}, 'loss')
            
            dnj = dnj + 1;
            cl_gene.dnlists.entries{i}{dnj} = upper(currcells{2});
            
        end
        
    end
    
    fclose(fid);
    
    if upj < 100000
        
        cl_gene.uplists.entries{i}(upj+1:end) = [];
        
    end
    
    cl_gene.uplists.numentries(i) = upj;
    
    discard = cellfun(@isempty, cl_gene.uplists.entries{i});
    
    cl_gene.uplists.entries{i}(discard) = [];
    
    cl_gene.uplists.numentries(i) = numel(cl_gene.uplists.entries{i});
    
    if dnj < 100000
        
        cl_gene.dnlists.entries{i}(dnj+1:end) = [];
        
    end
    
    cl_gene.dnlists.numentries(i) = dnj;
    
    discard = cellfun(@isempty, cl_gene.dnlists.entries{i});
    
    cl_gene.dnlists.entries{i}(discard) = [];
    
    cl_gene.dnlists.numentries(i) = numel(cl_gene.dnlists.entries{i});
    
    if mod(i, 100) == 0
        disp(['Finished ' num2str(i) ' of ' num2str(cl_gene.uplists.numterms) ' cell lines.']);
        toc;
        disp(' ');
        tic;
    end
    
end
disp(['Finished ' num2str(i) ' of ' num2str(cl_gene.uplists.numterms) ' cell lines.']);
toc;
disp(' ');

discard = cl_gene.uplists.numentries == 0;
cl_gene.uplists = listsdiscard(cl_gene.uplists, discard);

discard = cl_gene.dnlists.numentries == 0;
cl_gene.dnlists = listsdiscard(cl_gene.dnlists, discard);

cl_gene.upedges = lists2edges(cl_gene.uplists);
cl_gene.dnedges = lists2edges(cl_gene.dnlists);

gene_cl.edges = edgescombineupdn(sourcetargetswap(cl_gene.upedges), sourcetargetswap(cl_gene.dnedges));

gene_cl.cm = edges2cm(gene_cl.edges);

gene_cl = rmfield(gene_cl, 'edges');

save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_cl');
%}


