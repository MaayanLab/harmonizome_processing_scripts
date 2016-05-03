


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
% RS, rs number
% RSPOS, chromosome position
% GENEINFO, associated genes, formatted as GENESYM:GENEID|GENESYM:GENEID
% SAO, variant origin, 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both
% WGT, weight, 00 - unmapped, 1 - weight 1, 2 - weight 2, 3 - weight 3 or more (what does this mean?)
% CLNORIGIN, allele origin, 0 - unknown; 1 - germline; 2 - somatic; 4 - inherited; 8 - paternal; 16 - maternal; 32 - de-novo; 64 - biparental; 128 - uniparental; 256 - not-tested; 512 - tested-inconclusive; 1073741824 - other
% CLNSIG, variant clinical significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other
% CLNDSDB, variant disease database name
% CLNDSDBID variant disease database ID
% CLNDBN variant disease name

%{

% initialize edges structure
gene_atb.edges = edgesinit(114256, [], 'GeneSym', [], 'rsID_chromosome_position', [], 'GeneID', [], 'Phenotype', [], 'PhenotypeDatabase_PhenotypeID', [], 'ClinicalSignificance', false, []);
discard = false([gene_atb.edges.numedges 1]);
numgenes = zeros([gene_atb.edges.numedges 1]);
gene_atb.edges.targetid = zeros([gene_atb.edges.numedges 1]);


% read data
fid = fopen('input/clinvar.vcf', 'r');

for i = 1:1:67
    currline = fgetl(fid);
end

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb.edges.sourcedesc{i} = [currcells{3} '_' currcells{1} '_' currcells{2}];
    
%     subcells = strsplit(currcells{8}, ';', 'CollapseDelimiters', false);
    
    hidx = strfind(currcells{8}, 'GENEINFO=');
    substr = currcells{8}(hidx:end);
    b = find(substr == '=', 1, 'first');
    e = find(substr == ';', 1, 'first');
    geneinfo = substr(b+1:e-1);
    genes = strsplit(geneinfo, '|');
    numgenes(i) = numel(genes);
    genesym = cell(numgenes(i), 1);
    geneid = zeros([numgenes(i) 1]);
    for j = 1:1:numgenes(i)
        genesym{j} = genes{j}(1:find(genes{j} == ':')-1);
        geneid(j) = str2double(genes{j}(find(genes{j} == ':')+1:end));
    end
    
    hidx = strfind(currcells{8}, 'CLNSIG=');
    substr = currcells{8}(hidx:end);
    b = find(substr == '=', 1, 'first');
    e = find(substr == ';', 1, 'first');
    gene_atb.edges.targetid(i) = str2double(substr(b+1:e-1));
    
    hidx = strfind(currcells{8}, 'CLNDSDB=');
    substr = currcells{8}(hidx:end);
    b = find(substr == '=', 1, 'first');
    e = find(substr == ';', 1, 'first');
    diseasedb = strrep(strrep(strrep(substr(b+1:e-1), '_', ' '), '.', ''), ',', ':');
    
    diseasedbs = strsplit(diseasedb, ':');
    numdbs = numel(diseasedbs);
    hit = strcmp(diseasedbs, 'OMIM');
    if sum(hit) < 1
        hit = strcmp(diseasedbs, 'SNOMED CT');
        if sum(hit) < 1
            hit = strcmp(diseasedbs, 'Orphanet');
            if sum(hit) < 1
                hit = strcmp(diseasedbs, 'MedGen');
                if sum(hit) < 1
                    hit = strcmp(diseasedbs, 'GeneReviews');
                    if sum(hit) < 1
                        hit(1) = true;
                    end
                end
            end
        end
    end
    dbidx = find(hit, 1, 'first');
    diseasedb = diseasedbs{dbidx};
    
    
    hidx = strfind(currcells{8}, 'CLNDSDBID=');
    substr = currcells{8}(hidx:end);
    b = find(substr == '=', 1, 'first');
    e = find(substr == ';', 1, 'first');
    diseaseid = strrep(strrep(substr(b+1:e-1), '.', ''), ',', ':');
    
    diseaseids = strsplit(diseaseid, ':');
    diseaseid = diseaseids{dbidx};
    
    hidx = strfind(currcells{8}, 'CLNDBN=');
    substr = currcells{8}(hidx:end);
    b = find(substr == '=', 1, 'first');
    e = find(substr == ';', 1, 'first');
    disease = strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(substr(b+1:e-1), ',', ':'), '\x2c', ','), '_', ' '), 'not provided:', ''), 'not provided', ''), '.', ''), 'not specified:', ''), 'not specified', '');
    
    if numgenes(i) == 1 && ~isnan(geneid(1)) && numel(diseaseid) > 0 && numel(disease) > 0
        
        gene_atb.edges.source{i} = genesym{1};
        
        gene_atb.edges.sourceid(i) = geneid(1);
        
        gene_atb.edges.target{i} = disease;
        
        gene_atb.edges.targetdesc{i} = [diseasedb '_' diseaseid];
        
    else
        
        discard(i) = true;
        
    end
    
end

fclose(fid);

gene_atb.edges = edgesdiscard(gene_atb.edges, discard);

save('input/gene_phenotype_clinvar_filtered_20150415.mat', '-struct', 'gene_atb');
%}


gene_atb = load('input/gene_phenotype_clinvar_filtered_20150415.mat', '-mat', 'edges');


% % discard variants that are not pathogenic
% keep = gene_atb.edges.targetid == 4 | gene_atb.edges.targetid == 5;
% gene_atb.edges = edgesdiscard(gene_atb.edges, ~keep);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.edges = rmfield(gene_atb.edges, {'sourcedesc' 'sourcedescname' 'targetid' 'targetidname'});
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
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


