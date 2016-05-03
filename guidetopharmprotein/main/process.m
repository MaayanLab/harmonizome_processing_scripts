


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





% initialize edges structure
gene_atb = struct;
gene_atb.edges = edgesinit(15000, [], 'GeneSym', [], [], [], [], [], 'Ligand', [], [], [], [], false, []);


% read data
fid = fopen('input/interactions_20150616.txt', 'r');

currline = fgetl(fid);

i = 0;

while ~feof(fid)
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells) >= 13 && ~isempty(currcells{3}) && ~isempty(currcells{13})
        
        sourcecells = strsplit(currcells{3}, '|');
        targetcells = strsplit(currcells{13}, '|');
        
        for j = 1:1:numel(sourcecells)
            
            for k = 1:1:numel(targetcells)
                
                i = i + 1;
                
                gene_atb.edges.source{i} = sourcecells{j};
                
                gene_atb.edges.target{i} = targetcells{k};
                
            end
            
        end

    end
    
end

fclose(fid);

if i < gene_atb.edges.numedges
    discard = false([gene_atb.edges.numedges 1]);
    discard(i+1:end) = true;
    gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
end


% from website
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL2';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL3';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL4';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL5';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL7';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL8';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL11';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL13';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL14';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL17';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL22';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CXCL5';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CXCL6';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CXCL8';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CXCL11';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL2';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL5';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL7';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL11';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL14';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'ACKR1';
gene_atb.edges.target{gene_atb.edges.numedges} = 'CCL17';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'GHRHR';
gene_atb.edges.target{gene_atb.edges.numedges} = 'GHRH';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'GNRHR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'GNRH1';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'GNRHR2';
gene_atb.edges.target{gene_atb.edges.numedges} = 'GNRH2';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'NPY6R';
gene_atb.edges.target{gene_atb.edges.numedges} = 'NPY';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'NPY6R';
gene_atb.edges.target{gene_atb.edges.numedges} = 'PYY';
gene_atb.edges.numedges = gene_atb.edges.numedges + 1;
gene_atb.edges.source{gene_atb.edges.numedges} = 'NPY6R';
gene_atb.edges.target{gene_atb.edges.numedges} = 'PPY';


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map attribute identifiers that are genes to NCBI Entrez Gene Symbols and
% Gene IDs and discard non-protein ligands
[gene_atb.edges.target, gene_atb.edges.targetname, gene_atb.edges.targetid, gene_atb.edges.targetidname, discard] = genesymlookup(gene_atb.edges.target, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
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
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


