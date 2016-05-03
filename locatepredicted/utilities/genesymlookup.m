function [GeneSym, GeneSymName, GeneID, GeneIdName, Discard] = genesymlookup(genesym, geneid, family, species, usesynonyms, convert2human, filepath)

GeneSymName = 'GeneSym';
GeneIdName = 'GeneID';

if strcmp(species, 'both')
    
    gene_synonym_mouse = load([filepath '\' family '_synonym_resource_mouse.mat'], '-mat');
    gene_synonym_human = load([filepath '\' family '_synonym_resource_human.mat'], '-mat');
    
    gene_synonym_mouse.edges.target = upper(gene_synonym_mouse.edges.target);
    gene_synonym_mouse.lists.term = upper(gene_synonym_mouse.lists.term);
    gene_synonym_human.edges.target = upper(gene_synonym_human.edges.target);
    gene_synonym_human.lists.term = upper(gene_synonym_human.lists.term);

    if ~isempty(genesym)
        
        genesym = upper(genesym);

        NumGenes = numel(genesym);

        GeneSym = genesym;
        GeneID = -666*ones([NumGenes 1]);

        if usesynonyms

            [o1, o2] = ismember(genesym, gene_synonym_mouse.edges.target);
            o2(o2 == 0) = [];

            [o3, o4] = ismember(genesym, gene_synonym_mouse.lists.term);
            o4(o4 == 0) = [];
            
            [o5, o6] = ismember(genesym, gene_synonym_human.edges.target);
            o6(o6 == 0) = [];

            [o7, o8] = ismember(genesym, gene_synonym_human.lists.term);
            o8(o8 == 0) = [];

            GeneSym(o1) = gene_synonym_mouse.edges.source(o2);
            GeneID(o1) = gene_synonym_mouse.edges.sourceid(o2);

            GeneSym(o3) = gene_synonym_mouse.lists.term(o4);
            GeneID(o3) = gene_synonym_mouse.lists.termid(o4);
            
            GeneSym(o5) = gene_synonym_human.edges.source(o6);
            GeneID(o5) = gene_synonym_human.edges.sourceid(o6);

            GeneSym(o7) = gene_synonym_human.lists.term(o8);
            GeneID(o7) = gene_synonym_human.lists.termid(o8);

            Discard = ~o1 & ~o3 & ~o5 & ~o7;

        else

            [o1, o2] = ismember(genesym, gene_synonym_mouse.lists.term);
            o2(o2 == 0) = [];
            
            [o3, o4] = ismember(genesym, gene_synonym_human.lists.term);
            o4(o4 == 0) = [];

            GeneSym(o1) = gene_synonym_mouse.lists.term(o2);
            GeneID(o1) = gene_synonym_mouse.lists.termid(o2);
            
            GeneSym(o3) = gene_synonym_human.lists.term(o4);
            GeneID(o3) = gene_synonym_human.lists.termid(o4);

            Discard = ~o1 & ~o3;

        end

    else

        NumGenes = numel(geneid);

        GeneSym = repmat({'-666'}, NumGenes, 1);
        GeneID = geneid;

        [o1, o2] = ismember(geneid, gene_synonym_mouse.lists.termid);
        o2(o2 == 0) = [];
        
        [o3, o4] = ismember(geneid, gene_synonym_human.lists.termid);
        o4(o4 == 0) = [];

        GeneSym(o1) = gene_synonym_mouse.lists.term(o2);
        GeneID(o1) = gene_synonym_mouse.lists.termid(o2);
        
        GeneSym(o3) = gene_synonym_human.lists.term(o4);
        GeneID(o3) = gene_synonym_human.lists.termid(o4);

        Discard = ~o1 & ~o3;

    end
    
else
    
    gene_synonym = load([filepath '\' family '_synonym_resource_' species '.mat'], '-mat');
    
    gene_synonym.edges.target = upper(gene_synonym.edges.target);
    gene_synonym.lists.term = upper(gene_synonym.lists.term);

    if ~isempty(genesym)
        
        genesym = upper(genesym);

        NumGenes = numel(genesym);

        GeneSym = genesym;
        GeneID = -666*ones([NumGenes 1]);

        if usesynonyms

            [o1, o2] = ismember(genesym, gene_synonym.edges.target);
            o2(o2 == 0) = [];

            [o3, o4] = ismember(genesym, gene_synonym.lists.term);
            o4(o4 == 0) = [];

            GeneSym(o1) = gene_synonym.edges.source(o2);
            GeneID(o1) = gene_synonym.edges.sourceid(o2);

            GeneSym(o3) = gene_synonym.lists.term(o4);
            GeneID(o3) = gene_synonym.lists.termid(o4);

            Discard = ~o1 & ~o3;

        else

            [o1, o2] = ismember(genesym, gene_synonym.lists.term);
            o2(o2 == 0) = [];

            GeneSym(o1) = gene_synonym.lists.term(o2);
            GeneID(o1) = gene_synonym.lists.termid(o2);

            Discard = ~o1;

        end

    else

        NumGenes = numel(geneid);

        GeneSym = repmat({'-666'}, NumGenes, 1);
        GeneID = geneid;

        [o1, o2] = ismember(geneid, gene_synonym.lists.termid);
        o2(o2 == 0) = [];

        GeneSym(o1) = gene_synonym.lists.term(o2);
        GeneID(o1) = gene_synonym.lists.termid(o2);

        Discard = ~o1;

    end
    
end

if convert2human  % check this!
    
    if strcmp(species, 'both')
        human_other = load([filepath '\human_mouse_resource.mat'], '-mat');
    else
        human_other = load([filepath '\human_' species '_resource.mat'], '-mat');
    end
    
    [o1, o2] = ismember(GeneID, human_other.edges.targetid);
    o2(o2 == 0) = [];
    
%     uncertain = GeneID == -666;
    
    GeneSym(o1) = human_other.edges.source(o2);
    GeneID(o1) = human_other.edges.sourceid(o2);
    
%     GeneSym(uncertain) = {'-666'};
%     GeneID(uncertain) = -666;
    
%     Discard = GeneID == -666;
    
    Discard = Discard | (GeneID == -666) | strcmp(GeneSym, '-666');
    
end


