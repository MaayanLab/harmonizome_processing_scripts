


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
numcarriers = 24;
numenzymes = 233;
numtargets = 4207;
numtransporters = 118;
gene_atb.lists = listsinit(numcarriers+numenzymes+numtargets+numtransporters, [], 'GeneSym', [], 'UniprotAcc', [], [], [], 'DrugBankID', [], 'DrugBankID', [], [], false, []);
species = repmat({'-666'}, gene_atb.lists.numterms, 1);

% read data
fid = fopen('input/all_carrier_ids_all_20150416.txt', 'r');

currline = fgetl(fid);

j = 0;

for i = 1:1:numcarriers
    
    j = j + 1;
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells{1}) > 0
        gene_atb.lists.termdesc{j} = currcells{1};
    end
    
    if numel(currcells{3}) > 0
        gene_atb.lists.term{j} = currcells{3};
    end
    
    if numel(currcells{3}) > 0
        species{j} = currcells{12};
    end
    
    subcells = strsplit(currcells{13}, '; ');
    
    discard = cellfun(@numel, subcells) == 0;
    
    gene_atb.lists.entries{j} = unique(subcells(~discard))';
    gene_atb.lists.entrydescs{j} = gene_atb.lists.entries{j};
    gene_atb.lists.numentries(j) = numel(gene_atb.lists.entries{j});
    
end

fclose(fid);


fid = fopen('input/all_enzyme_ids_all_20150416.txt', 'r');

currline = fgetl(fid);

for i = 1:1:numenzymes
    
    j = j + 1;
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells{1}) > 0
        gene_atb.lists.termdesc{j} = currcells{1};
    end
    
    if numel(currcells{3}) > 0
        gene_atb.lists.term{j} = currcells{3};
    end
    
    if numel(currcells{3}) > 0
        species{j} = currcells{12};
    end
    
    subcells = strsplit(currcells{13}, '; ');
    
    discard = cellfun(@numel, subcells) == 0;
    
    gene_atb.lists.entries{j} = unique(subcells(~discard))';
    gene_atb.lists.entrydescs{j} = gene_atb.lists.entries{j};
    gene_atb.lists.numentries(j) = numel(gene_atb.lists.entries{j});
    
end

fclose(fid);


fid = fopen('input/all_target_ids_all_20150416.txt', 'r');

currline = fgetl(fid);

for i = 1:1:numtargets
    
    j = j + 1;
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells{1}) > 0
        gene_atb.lists.termdesc{j} = currcells{1};
    end
    
    if numel(currcells{3}) > 0
        gene_atb.lists.term{j} = currcells{3};
    end
    
    if numel(currcells{3}) > 0
        species{j} = currcells{12};
    end
    
    subcells = strsplit(currcells{13}, '; ');
    
    discard = cellfun(@numel, subcells) == 0;
    
    gene_atb.lists.entries{j} = unique(subcells(~discard))';
    gene_atb.lists.entrydescs{j} = gene_atb.lists.entries{j};
    gene_atb.lists.numentries(j) = numel(gene_atb.lists.entries{j});
    
end

fclose(fid);


fid = fopen('input/all_transporter_ids_all_20150416.txt', 'r');

currline = fgetl(fid);

for i = 1:1:numtransporters
    
    j = j + 1;
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells{1}) > 0
        gene_atb.lists.termdesc{j} = currcells{1};
    end
    
    if numel(currcells{3}) > 0
        gene_atb.lists.term{j} = currcells{3};
    end
    
    if numel(currcells{3}) > 0
        species{j} = currcells{12};
    end
    
    subcells = strsplit(currcells{13}, '; ');
    
    discard = cellfun(@numel, subcells) == 0;
    
    gene_atb.lists.entries{j} = unique(subcells(~discard))';
    gene_atb.lists.entrydescs{j} = gene_atb.lists.entries{j};
    gene_atb.lists.numentries(j) = numel(gene_atb.lists.entries{j});
    
end

fclose(fid);

discard = ~strcmpi(species, 'Human') | (strcmp(gene_atb.lists.term, '-666') & strcmp(gene_atb.lists.termdesc, '-666')) | gene_atb.lists.numentries == 0;
gene_atb.lists = listsdiscard(gene_atb.lists, discard);

gene_atb.edges = lists2edges(gene_atb.lists);

gene_atb.lists = edges2lists(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = lists2cm(gene_atb.lists);
gene_atb = rmfield(gene_atb, {'lists' 'edges'});


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% map attribute names to attributes
gene_atb.cm.entryname = 'Drug';
for i = 1:1:gene_atb.cm.numentries
    
    try
        
        u = urlread(['http://www.drugbank.ca/drugs/' gene_atb.cm.entrydesc{i}]);

        b = strfind(u, '<title>DrugBank: ') + 17;
        e = strfind(u, [' (' gene_atb.cm.entrydesc{i} ')</title>']) - 1;

        gene_atb.cm.entry{i} = u(b:e);

        pause(1.5);
        
    catch
        
        try
            
            u = urlread(['http://www.drugbank.ca/drugs/' gene_atb.cm.entrydesc{i}]);

            b = strfind(u, '<title>DrugBank: ') + 17;
            e = strfind(u, [' (' gene_atb.cm.entrydesc{i} ')</title>']) - 1;

            gene_atb.cm.entry{i} = u(b:e);

            pause(1.5);
            
        catch
            
            pause(1.5);
            
        end
        
    end
    
end


% merge columns corresponding to the same attribute
if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


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


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


