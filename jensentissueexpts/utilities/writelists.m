function writelists(lists, filename, formattype)

if isfield(lists, 'termdesc') && isfield(lists, 'termid')
    
    lists.termid = cellfun(@num2str, num2cell(lists.termid), 'UniformOutput', false);
    lists.description = [lists.termdesc lists.termid];
    lists.description = mat2cell(lists.description, ones([lists.numterms 1]));
    lists.description = cellfun(@(x) [x{1} '_' x{2}], lists.description, 'UniformOutput', false);
    
elseif isfield(lists, 'termid')
    
    lists.description = cellfun(@num2str, num2cell(lists.termid), 'UniformOutput', false);
    
elseif isfield(lists, 'termdesc')
    
    lists.description = lists.termdesc;
    
else
    
    lists.description = repmat({'NA'}, lists.numterms, 1);
    
end



wf_cr = fopen([filename '_crisp.gmt'], 'w');

for i = 1:1:lists.numterms

%     fmt = strjoin([{'%s\t'} {'%s\t'} repmat({'%s\t'}, 1, lists.numentries(i)-1) {'%s\r\n'}], '');
%     fprintf(wf_cr, fmt, lists.term{i}, lists.description{i}, lists.entries{i}');
    
    fprintf(wf_cr, '%s\r\n', strjoin([lists.term(i) lists.description(i) lists.entries{i}'], '\t'));

end

fclose(wf_cr);
    
if isfield(lists, 'weights') && ~strcmp(formattype, 'crisp')

    wf_fz = fopen([filename '_fuzzy.gmt'], 'w');
    
    for i = 1:1:lists.numterms
        
%         fmt = strjoin([{'%s\t'} {'%s\t'} repmat({'%s,%f\t'}, 1, lists.numentries(i)-1) {'%s,%f\t'}], '');
% 
%         entryweightarray = reshape([lists.entries{i}'; num2cell(lists.weights{i}')], 1, 2*lists.numentries(i));
% %         entryweightarray = reshape([lists.entries{i}'; cellfun(@num2str, num2cell(lists.weights{i}'), 'UniformOutput', false)], 1, 2*lists.numentries(i));
%         
%         fprintf(wf_fz, fmt, lists.term{i}, lists.description{i}, entryweightarray);
        
        weights = cellfun(@num2str, num2cell(lists.weights{i}), 'UniformOutput', false);
        entryweightpairs = [lists.entries{i} weights];
        entryweightpairs = mat2cell(entryweightpairs, ones([lists.numentries(i) 1]));
        entryweightpairs = cellfun(@(x) [x{1} ',' x{2}], entryweightpairs, 'UniformOutput', false);
        
        fprintf(wf_fz, '%s\r\n', strjoin([lists.term(i) lists.description(i) entryweightpairs'], '\t'));

    end

    fclose(wf_fz);
    
end


