


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




% get file list
dircontents = dir('input/hmdb_proteins/');
filename = {dircontents(:).name}';
filename([dircontents(:).isdir]) = [];
clear dircontents;

% initialize lists structure
gene_atb.lists = listsinit(numel(filename), [], 'GeneSym', [], 'HMDB ACC', [], [], [], 'Metabolite', [], 'HMDB ACC', [], [], false, []);
discard = false(gene_atb.lists.numterms, 1);

% read data
for i = 1:1:gene_atb.lists.numterms
    
    fid = fopen(['input/hmdb_proteins/' filename{i}], 'r');
    
    line1 = fgetl(fid);
    line2 = fgetl(fid);
    line3 = fgetl(fid);
    
    target = '  <accession>';
    tlen = numel(target);
    while numel(line3) < tlen || ~strcmp(line3(1:tlen), target)

        line1 = line2;
        line2 = line3;
        line3 = fgetl(fid);
        
    end
    
    gene_atb.lists.termdesc{i} = strrep(strrep(line3, '  <accession>', ''), '</accession>', '');
    
    line1 = line2;
    line2 = line3;
    line3 = fgetl(fid);
    
    target = '  <gene_name>';
    tlen = numel(target);
    while (numel(line3) < tlen || ~strcmp(line3(1:tlen), target)) && ~strcmp(line3, '  <gene_name/>')

        line1 = line2;
        line2 = line3;
        line3 = fgetl(fid);
        
    end
    
    if strcmp(line3, '  <gene_name/>')
        
        gene_atb.lists.term{i} = '-666';
        discard(i) = true;
        
    else
    
        gene_atb.lists.term{i} = strrep(strrep(line3, '  <gene_name>', ''), '</gene_name>', '');
        
    end
    
    line1 = line2;
    line2 = line3;
    line3 = fgetl(fid);
    
    entries = cell(10000, 1);
    entrydescs = cell(10000, 1);
    j = 0;
    
    while ~feof(fid)
        
        while ~strcmp(line1, '    <metabolite>') && ~strcmp(line1, '  </metabolite_associations>')

            line1 = line2;
            line2 = line3;
            line3 = fgetl(fid);

        end

        if strcmp(line1, '  </metabolite_associations>')
            
            break;
            
        else
            
            j = j + 1;

            entries{j} = strrep(strrep(line3, '      <name>', ''), '</name>', '');
            entrydescs{j} = strrep(strrep(line2, '      <accession>', ''), '</accession>', '');

            line1 = line2;
            line2 = line3;
            line3 = fgetl(fid);
            
        end
    
    end
    
    fclose(fid);
    
    gene_atb.lists.entries{i} = entries(1:j);
    gene_atb.lists.entrydescs{i} = entrydescs(1:j);
    
end

gene_atb.lists = listsdiscard(gene_atb.lists, discard);


% convert to matrix format (attribute table)
gene_atb.cm = lists2cm(gene_atb.lists);
gene_atb = rmfield(gene_atb, 'lists');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% view distributions of row and col stats
% [~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
% close force all;


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
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


