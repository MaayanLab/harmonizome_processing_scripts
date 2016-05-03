


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

% initialize edges structure
gene_atb.edges = edgesinit(75830, [], 'UniprotAcc', [], 'Species', [], [], [], 'UniprotAcc', [], 'Species', [], [], false, []);
detectionmethod = cell(gene_atb.edges.numedges, 1);
interactiontype = cell(gene_atb.edges.numedges, 1);
qualitystatus = cell(gene_atb.edges.numedges, 1);


% read human data
fid = fopen('input/dip20160114.txt', 'r');
currline = fgetl(fid);
headers = strsplit(currline, '\t', 'CollapseDelimiters', false)';
for i = 1:1:gene_atb.edges.numedges
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    
    b = find(currcells{1} == ':', 1, 'last');
    if ~isempty(b)
        gene_atb.edges.source{i} = currcells{1}(b+1:end);
    else
        gene_atb.edges.source{i} = currcells{1};
    end
    
    b = find(currcells{2} == ':', 1, 'last');
    if ~isempty(b)
        gene_atb.edges.target{i} = currcells{2}(b+1:end);
    else
        gene_atb.edges.target{i} = currcells{2};
    end
    
    b = find(currcells{15} == ':', 1, 'first');
    e = find(currcells{15} == '|', 1, 'first');
    if ~isempty(e)
        qualitystatus{i} = currcells{15}(b+1:e-1);
    elseif ~isempty(b)
        qualitystatus{i} = currcells{15}(b+1:end);
    else
        qualitystatus{i} = currcells{15};
    end
    
    pl = find(currcells{10} == '(', 1, 'first');
    pr = find(currcells{10} == ')', 1, 'first');
    gene_atb.edges.sourcedesc{i} = currcells{10}(pl+1:pr-1);
    
    pl = find(currcells{11} == '(', 1, 'first');
    pr = find(currcells{11} == ')', 1, 'first');
    gene_atb.edges.targetdesc{i} = currcells{11}(pl+1:pr-1);
    
    pl = find(currcells{7} == '(', 1, 'first');
    pr = find(currcells{7} == ')', 1, 'first');
    detectionmethod{i} = currcells{7}(pl+1:pr-1);
    
    pl = find(currcells{12} == '(', 1, 'first');
    pr = find(currcells{12} == ')', 1, 'first');
    interactiontype{i} = currcells{12}(pl+1:pr-1);
end
fclose(fid);

save('input/gene_gene_dip_unmodified_20160114.mat', '-struct', 'gene_atb');
save('input/supp_dip_unmodified_20160114.mat', '-mat', 'interactiontype', 'detectionmethod');
%}


gene_atb = load('input/gene_gene_dip_unmodified_20160114.mat', '-mat', 'edges');
load('input/supp_dip_unmodified_20160114.mat', '-mat', 'interactiontype', 'detectionmethod');


% discard interactions that are not human or not complete
discard = strcmp(gene_atb.edges.source, '-') | strcmp(gene_atb.edges.target, '-') | strcmp(gene_atb.edges.source, '') | strcmp(gene_atb.edges.target, '') | ~strcmp(gene_atb.edges.sourcedesc, 'Homo sapiens') | ~strcmp(gene_atb.edges.targetdesc, 'Homo sapiens');
% discard = strcmp(gene_atb.edges.source, '-') | strcmp(gene_atb.edges.target, '-') | strcmp(gene_atb.edges.source, '') | strcmp(gene_atb.edges.target, '') | ~ismember(gene_atb.edges.sourcedesc, {'Homo sapiens' 'Mus musculus'}) | ~ismember(gene_atb.edges.targetdesc, {'Homo sapiens' 'Mus musculus'});
detectionmethod(discard) = [];
interactiontype(discard) = [];
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
gene_atb.edges = rmfield(gene_atb.edges, {'sourcedesc' 'sourcedescname' 'targetdesc' 'targetdescname'});


% discard interactions that are not protein-protein
[u_interactiontype, ui, ri] = unique(interactiontype);
c_interactiontype = zeros(size(u_interactiontype));
for i = 1:1:numel(u_interactiontype)
c_interactiontype(i) = sum(ri==i);
end
[c_interactiontype, si] = sort(c_interactiontype, 'descend');
u_interactiontype = u_interactiontype(si);
discard = ismember(interactiontype, {'association' 'genetic interference' 'physical association' 'physical interaction'});
% discard = ismember(interactiontype, {'association' 'genetic interference' 'physical association'});
% discard = ismember(interactiontype, {'association' 'genetic interference'});
detectionmethod(discard) = [];
interactiontype(discard) = [];
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
% physical association = not necessarily in same complex
% physical interaction = in same complex but not necessarily directly interacting


% discard interactions with bad scores
% discard = ~strcmp(qualitystatus, 'core');
% detectionmethod(discard) = [];
% interactiontype(discard) = [];
% gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
uniprot_entrez = load('C:/Users/Andrew/Dropbox/GeneInfo/uniprot_entrez_resource_human_withisoforms.mat', '-mat');
gene_atb.edges.sourcedesc = gene_atb.edges.source;
gene_atb.edges.sourcedescname = gene_atb.edges.sourcename;
for i = 1:1:gene_atb.edges.numedges
    subcells = strsplit(gene_atb.edges.sourcedesc{i}, '-');
    gene_atb.edges.sourcedesc{i} = subcells{1};
end
gene_atb.edges.sourceid = -666*ones([gene_atb.edges.numedges 1]);
gene_atb.edges.sourceidname = 'GeneID';
gene_atb.edges.source = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene_atb.edges.sourcename = 'GeneSym';
[o1, o2] = ismember(gene_atb.edges.sourcedesc, uniprot_entrez.edges.source);
o2(o2 == 0) = [];
gene_atb.edges.source(o1) = uniprot_entrez.edges.target(o2);
gene_atb.edges.sourceid(o1) = uniprot_entrez.edges.targetid(o2);
gene_atb.edges = edgesdiscard(gene_atb.edges, ~o1);
detectionmethod(~o1) = [];
interactiontype(~o1) = [];


% map attribute identifiers to NCBI Entrez Gene Symbols and Gene IDs
uniprot_entrez = load('C:/Users/Andrew/Dropbox/GeneInfo/uniprot_entrez_resource_human_withisoforms.mat', '-mat');
gene_atb.edges.targetdesc = gene_atb.edges.target;
gene_atb.edges.targetdescname = gene_atb.edges.targetname;
for i = 1:1:gene_atb.edges.numedges
    subcells = strsplit(gene_atb.edges.targetdesc{i}, '-');
    gene_atb.edges.targetdesc{i} = subcells{1};
end
gene_atb.edges.targetid = -666*ones([gene_atb.edges.numedges 1]);
gene_atb.edges.targetidname = 'GeneID';
gene_atb.edges.target = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene_atb.edges.targetname = 'GeneSym';
[o1, o2] = ismember(gene_atb.edges.targetdesc, uniprot_entrez.edges.source);
o2(o2 == 0) = [];
gene_atb.edges.target(o1) = uniprot_entrez.edges.target(o2);
gene_atb.edges.targetid(o1) = uniprot_entrez.edges.targetid(o2);
gene_atb.edges = edgesdiscard(gene_atb.edges, ~o1);
detectionmethod(~o1) = [];
interactiontype(~o1) = [];


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% get edges from source to target and from target to source
gene_atb.edges = edges2symedges(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% cluster and view matrix
gene_atb.cm = cmcluster_nogram(gene_atb.cm, false);


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
gene_atb.cm = cmcluster_nogram(gene_atb.cm, false);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


