


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




% initialize lists structure
atb_gene.lists = listsinit(315, [], 'PMID and descriptive information', [], [], [], [], [], 'GeneSym', [], [], [], [], true, []);

% read data
fid = fopen('input/ESCAPE.txt', 'r');

for i = 1:1:atb_gene.lists.numterms
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    atb_gene.lists.term{i} = strrep(strrep(currcells{1}, '_UP', ''), '_DOWN', '');
    atb_gene.lists.entries{i} = currcells(3:end)';
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
    
    if numel(strfind(currcells{1}, 'DOWN')) > 0
        atb_gene.lists.weights{i} = -1*ones(atb_gene.lists.numentries(i), 1);
    else
        atb_gene.lists.weights{i} = ones(atb_gene.lists.numentries(i), 1);
    end
    
end

fclose(fid);


% convert to edges to merge up and down lists
atb_gene.edges = lists2edges(atb_gene.lists);

% convert to matrix format (attribute table)
gene_atb.cm = cmtranspose(edges2cm(atb_gene.edges));
clear atb_gene;


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% reformatting labels
gene_atb.cm.entry = cellfun(@(x) strrep(strrep(strrep(x, '-', '_'), 'RNAiScreen_', 'RNAi_-666_'), 'CHiP', 'ChIP'), gene_atb.cm.entry, 'UniformOutput', false); 

gene_atb.cm.entry{strcmp(gene_atb.cm.entry, 'hESC_H3K27me3_20682450')} = 'ChIP_H3K27me3_20682450_humanESC';
gene_atb.cm.entry{strcmp(gene_atb.cm.entry, 'mESC_H3K27me3_17603471')} = 'ChIP_H3K27me3_17603471_mouseESC';
gene_atb.cm.entry{strcmp(gene_atb.cm.entry, 'mESC_H3K36me3_18692474')} = 'ChIP_H3K36me3_18692474_mouseESC';
gene_atb.cm.entry{strcmp(gene_atb.cm.entry, 'mESC_H3K9me3_19884255')} = 'ChIP_H3K9me3_19884255_mouseESC';
gene_atb.cm.entry{strcmp(gene_atb.cm.entry, 'mMEF_K27me3_17603471')} = 'ChIP_K27me3_17603471_mouseMEF';
gene_atb.cm.entry{strcmp(gene_atb.cm.entry, 'mNPC_K27me3_17603471')} = 'ChIP_K27me3_17603471_mouseNPC';
gene_atb.cm.entry{strcmp(gene_atb.cm.entry, 'mWholeBrain_H3K27me3_18600261')} = 'ChIP_H3K27me3_18600261_mouseWholeBrain';

for i = 1:1:gene_atb.cm.numentries
    si = regexp(gene_atb.cm.entry{i}, '_');
    prefix = gene_atb.cm.entry{i}(1:si(1)-1);
    if ~ismember(prefix, {'RNAi', 'ChIP', 'Protein'});
        gene_atb.cm.entry{i} = ['mRNA_' gene_atb.cm.entry{i}];
    end
end

gene_atb.cm.entryname = 'Evidence Type_Gene Symbol_Pubmed ID';


% view distributions of row and col stats
% [~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


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
% close force all;


% view distributions of row and col stats
% [~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


