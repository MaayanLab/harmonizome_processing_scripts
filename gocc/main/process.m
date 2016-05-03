


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
% aspects: P = biological process, F = molecular function, C = cellular component



%{
% initialize edges structure
gene_atb.edges = edgesinit(354285, [], 'UniProtACC', [], 'UniProtACC', [], [], [], 'GOID', [], 'GOID', [], [], false, []);

aspect = cell(gene_atb.edges.numedges, 1);

% read data
fid = fopen('input/gene_association_20150331.goa_ref_human', 'r');

for i = 1:1:34
    currline = fgetl(fid);
end

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb.edges.source{i} = currcells{2};
    
    gene_atb.edges.sourcedesc{i} = currcells{2};
    
    gene_atb.edges.target{i} = currcells{5};
    
    gene_atb.edges.targetdesc{i} = currcells{5};
    
    aspect{i} = currcells{9};
    
end

fclose(fid);

discard = ~strcmp(aspect, 'P');
gene_bp.edges = edgesdiscard(gene_atb.edges, discard);
save('input/gene_bp_go_unmodified_20150331.mat', '-struct', 'gene_bp');

discard = ~strcmp(aspect, 'F');
gene_mf.edges = edgesdiscard(gene_atb.edges, discard);
save('input/gene_mf_go_unmodified_20150331.mat', '-struct', 'gene_mf');

discard = ~strcmp(aspect, 'C');
gene_cc.edges = edgesdiscard(gene_atb.edges, discard);
save('input/gene_cc_go_unmodified_20150331.mat', '-struct', 'gene_cc');
%}


% load edges
gene_atb = load('input/gene_cc_go_unmodified_20150331.mat', '-mat', 'edges');


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% map protein identifiers to NCBI Entrez Gene Symbols and Gene IDs
uniprot_entrez = load([mappingfilespath '\uniprot_entrez_resource_human.mat'], '-mat');
[o1, o2] = ismember(gene_atb.cm.term, uniprot_entrez.edges.source);
o2(o2 == 0) = [];
gene_atb.cm.termname = 'GeneSym';
gene_atb.cm.term(o1) = uniprot_entrez.edges.target(o2);
gene_atb.cm.termid = -666*ones([gene_atb.cm.numterms 1]);
gene_atb.cm.termidname = 'GeneID';
gene_atb.cm.termid(o1) = uniprot_entrez.edges.targetid(o2);
gene_atb.cm = cmrowdiscard(gene_atb.cm, ~o1);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% map attributes to ontology
atb_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');
id_alt = load('output/gene_attribute_matrix_imported.mat', '-mat');

[o1, o2] = ismember(gene_atb.cm.entrydesc, id_alt.edges.targetdesc);
o2(o2 == 0) = [];

[o3, o4] = ismember(gene_atb.cm.entrydesc, atb_atb.cm.entrydesc);
o4(o4 == 0) = [];

gene_atb.cm.entryname = 'Cellular Component';
gene_atb.cm.entry(o1) = id_alt.edges.source(o2);
gene_atb.cm.entrydesc(o1) = id_alt.edges.sourcedesc(o2);
gene_atb.cm.entry(o3) = atb_atb.cm.entry(o4);
gene_atb.cm.entrydesc(o3) = atb_atb.cm.entrydesc(o4);

discard = ~o1 & ~o3;
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);


% merge columns corresponding to the same attribute
if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% ensure propagation of gene associations up ontology to ancestor terms
gene_atb.cm = cm2em_cosine(gene_atb.cm, atb_atb.cm, false, []);
gene_atb.cm.matrix = double(gene_atb.cm.matrix > 0);
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);


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


