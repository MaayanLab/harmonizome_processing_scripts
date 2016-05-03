


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




gene_cl.edges = edgesinit(1337971, [], 'GeneSym', [], 'GenomeWideScreen', [], [], [], 'CellLine', [], 'Tissue', [], 'CellLineID', false, []);

fid = fopen('input/CosmicCLP_MutantExport.tsv', 'r');

currline = fgetl(fid);

tic;
for i = 1:1:gene_cl.edges.numedges
    
    currline = upper(fgetl(fid));
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    e = find(currcells{1} == '_', 1, 'first') - 1;
    
    if isempty(e)
        gene_cl.edges.source{i} = currcells{1};
    else
        gene_cl.edges.source{i} = currcells{1}(1:e);
    end

    gene_cl.edges.sourcedesc{i} = currcells{12};
    
    gene_cl.edges.target{i} = currcells{5};
    gene_cl.edges.targetdesc{i} = currcells{8};
    gene_cl.edges.targetid(i) = str2double(currcells{6});
    
    if mod(i,130000) == 0
        disp(['Finished ' num2str(i) ' of ' num2str(gene_cl.edges.numedges) '.']);
        toc;
        disp('');
        tic;
    end
    
end
disp(['Finished ' num2str(i) ' of ' num2str(gene_cl.edges.numedges) '.']);
toc;
disp('');
fclose(fid);

gene_cl.edges.targetdesc = lower(gene_cl.edges.targetdesc);
hit = strcmp(gene_cl.edges.targetdesc, 'ns');
gene_cl.edges.targetdesc(hit) = {'-666'};

gene_cl.cm = edges2cm(gene_cl.edges);

gene_cl = rmfield(gene_cl, 'edges');

% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
[gene_cl.cm.term, gene_cl.cm.termname, gene_cl.cm.termid, gene_cl.cm.termidname, discard] = genesymlookup(gene_cl.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_cl.cm = cmrowdiscard(gene_cl.cm, discard);

% cluster and view data
cgo = cm2clustergram(gene_cl.cm, 'none', 'all', 'cosine', 'average');

close force all;

gene_cl.cm = conmatmap(gene_cl.cm, cgo.RowLabels, cgo.ColumnLabels');

HeatMap(gene_cl.cm.matrix, 'Colormap', redbluecmap);

% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_cl');


