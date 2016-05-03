function writeedges(edges, filename)

if ~isfield(edges, 'sourcedesc')
    
    edges.sourcedescname = 'NA';
    edges.sourcedesc = repmat({'na'}, edges.numedges, 1);
    
end

if ~isfield(edges, 'sourceid')
    
    edges.sourceidname = 'NA';
    edges.sourceid = -666*ones([edges.numedges 1]);
    
end

if ~isfield(edges, 'targetdesc')
    
    edges.targetdescname = 'NA';
    edges.targetdesc = repmat({'na'}, edges.numedges, 1);
    
end

if ~isfield(edges, 'targetid')
    
    edges.targetidname = 'NA';
    edges.targetid = -666*ones([edges.numedges 1]);
    
end

if ~isfield(edges, 'weight')
    
    edges.weight = ones([edges.numedges 1]);
    
end

fid = fopen([filename '_edges.txt'], 'w');

fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n', 'source', 'source_desc', 'source_id', 'target', 'target_desc', 'target_id', 'weight');

fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n', edges.sourcename, edges.sourcedescname, edges.sourceidname, edges.targetname, edges.targetdescname, edges.targetidname, 'weight');

fmt = '%s\t%s\t%i\t%s\t%s\t%i\t%f\r\n';

for i = 1:1:edges.numedges
    
    fprintf(fid, fmt, edges.source{i}, edges.sourcedesc{i}, edges.sourceid(i), edges.target{i}, edges.targetdesc{i}, edges.targetid(i), edges.weight(i));
    
end

fclose(fid);

%%% no target_desc and no weight
% fid = fopen([filename '_edges.txt'], 'w');
% 
% fprintf(fid, '%s\t%s\t%s\t%s\t%s\r\n', 'source', 'source_desc', 'source_id', 'target', 'target_id');
% 
% fprintf(fid, '%s\t%s\t%s\t%s\t%s\r\n', edges.sourcename, edges.sourcedescname, edges.sourceidname, edges.targetname, edges.targetidname);
% 
% fmt = '%s\t%s\t%i\t%s\t%i\r\n';
% 
% for i = 1:1:edges.numedges
%     
%     fprintf(fid, fmt, edges.source{i}, edges.sourcedesc{i}, edges.sourceid(i), edges.target{i}, edges.targetid(i));
%     
% end
% 
% fclose(fid);

%%% no weight
% fid = fopen([filename '_edges.txt'], 'w');
% 
% fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\r\n', 'source', 'source_desc', 'source_id', 'target', 'targetdesc', 'target_id');
% 
% fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\r\n', edges.sourcename, edges.sourcedescname, edges.sourceidname, edges.targetname, edges.targetdescname, edges.targetidname);
% 
% fmt = '%s\t%s\t%i\t%s\t%s\t%i\r\n';
% 
% for i = 1:1:edges.numedges
%     
%     fprintf(fid, fmt, edges.source{i}, edges.sourcedesc{i}, edges.sourceid(i), edges.target{i}, edges.targetdesc{i}, edges.targetid(i));
%     
% end
% 
% fclose(fid);


