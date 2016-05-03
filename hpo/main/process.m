


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
phenotype_phenotype.edges = edgesinit(1000000, [], 'Phenotype', [], 'HPOID', [], [], [], 'Phenotype', [], 'HPOID', [], [], false, []);
hpoid_althpoid.edges = edgesinit(1000000, [], 'Phenotype', [], 'HPOID', [], [], [], 'Phenotype', [], 'HPOID', [], [], false, []);


% read data
fid = fopen('input/hp_20150324.obo', 'r');

i = 0;
h = 0;

while ~feof(fid)
    
    currline = fgetl(fid);
    
    if strcmp(currline, '[Term]')
        
        currline = fgetl(fid);
        
        source = cell(1000, 1);
        sourcedesc = cell(1000, 1);
        j = 0;
        
        alttargetdesc = cell(1000, 1);
        k = 0;
        
        while ~isempty(currline)
            
            currcells = strsplit(currline, ': ');
            
            if numel(currcells) > 2
                currcells{2} = strjoin(currcells(2:end), ': ');
            end
            
            switch currcells{1}
                
                case 'id'
                    
                    targetdesc = currcells{2};
                    
                case 'name'
                    
                    target = lower(currcells{2});
                    
                case 'alt_id'
                    
                    k = k + 1;
                    alttargetdesc{k} = currcells{2};
                    
                case 'is_a'
                    
                    j = j + 1;
                    subcells = strsplit(currcells{2}, ' ! ');
                    source{j} = lower(subcells{2});
                    sourcedesc{j} = subcells{1};
                    
                otherwise
                    
%                     disp('skip');
                    
            end
            
            currline = fgetl(fid);
                    
        end
        
        if j < 1000
            source(j+1:end) = [];
            sourcedesc(j+1:end) = [];
        end
        numsources = j;
        
        if k < 1000
            alttargetdesc(k+1:end) = [];
        end
        numalttargets = k;
        
        for j = 1:1:numsources
            
            i = i + 1;
            
            phenotype_phenotype.edges.source{i} = source{j};
            phenotype_phenotype.edges.sourcedesc{i} = sourcedesc{j};
            phenotype_phenotype.edges.target{i} = target;
            phenotype_phenotype.edges.targetdesc{i} = targetdesc;
            
        end
        
        for k = 1:1:numalttargets
            
            h = h + 1;
            
            hpoid_althpoid.edges.source{h} = target;
            hpoid_althpoid.edges.sourcedesc{h} = targetdesc;
            hpoid_althpoid.edges.target{h} = target;
            hpoid_althpoid.edges.targetdesc{h} = alttargetdesc{k};
            
        end
            
    end
    
end

fclose(fid);

if i < phenotype_phenotype.edges.numedges
    discard = false([phenotype_phenotype.edges.numedges 1]);
    discard(i+1:end) = true;
    phenotype_phenotype.edges = edgesdiscard(phenotype_phenotype.edges, discard);
end

if h < hpoid_althpoid.edges.numedges
    discard = false([hpoid_althpoid.edges.numedges 1]);
    discard(h+1:end) = true;
    hpoid_althpoid.edges = edgesdiscard(hpoid_althpoid.edges, discard);
end

save('output/gene_attribute_matrix_imported.mat', '-struct', 'hpoid_althpoid');



% convert to connectivity matrix of directed "is_a" relationships
phenotype_phenotype.cm = edges2cm(phenotype_phenotype.edges);
phenotype_phenotype = rmfield(phenotype_phenotype, 'edges');



% view distributions of row and col stats
[~] = cmrowcolstats(phenotype_phenotype.cm, true, 1);



% cluster and view matrix
phenotype_phenotype.cm = cmcluster(phenotype_phenotype.cm, true);



% convert to sparse matrix
phenotype_phenotype.cm.matrix = sparse(phenotype_phenotype.cm.matrix);



% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'phenotype_phenotype');



% convert to connectivy matrix of parents and all descendants
phenotype_phenotype.cm = directednetwork2subnodes(phenotype_phenotype.cm);



% view distributions of row and col stats
[~] = cmrowcolstats(phenotype_phenotype.cm, true, 1);



% cluster and view matrix
phenotype_phenotype.cm = cmcluster(phenotype_phenotype.cm, true);



% convert to sparse matrix
phenotype_phenotype.cm.matrix = sparse(phenotype_phenotype.cm.matrix);



% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'phenotype_phenotype');


