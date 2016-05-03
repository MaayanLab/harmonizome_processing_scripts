function [GeneSym, GeneSymName, GeneID, GeneIdName, Discard] = uniprot2entrez(uniprotacc, species, filepath)

N = numel(uniprotacc);
GeneSym = repmat({'-666'}, N, 1);
GeneSymName = 'GeneSym';
GeneID = -666*ones(N, 1);
GeneIdName = 'GeneID';

if strcmp(species, 'both')
    uniprot_entrez = load([filepath '\uniprot_entrez_resource_mouse_withisoforms.mat'], '-mat', 'edges');
    load([filepath '\uniprot_entrez_resource_human_withisoforms.mat'], '-mat', 'edges');
    uniprot_entrez.edges = edgesvertcat(uniprot_entrez.edges, edges);
    clear edges;
else
    uniprot_entrez = load([filepath '\uniprot_entrez_resource_' species '_withisoforms.mat'], '-mat', 'edges');
end

[o1, o2] = ismember(uniprotacc, uniprot_entrez.edges.source);
o2(o2 == 0) = [];
GeneSym(o1) = uniprot_entrez.edges.target(o2);
GeneID(o1) = uniprot_entrez.edges.targetid(o2);

Discard = ~o1;


