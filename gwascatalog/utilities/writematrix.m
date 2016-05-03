function writematrix(filename, matrix, collabelname1, collabel1, collabelname2, collabel2, collabelname3, collabel3, rowlabelname1, rowlabel1, rowlabelname2, rowlabel2, rowlabelname3, rowlabel3)

if issparse(matrix)
    matrix = full(matrix);
end

[Nr, Nc] = size(matrix);

if isempty(collabel1)
    collabelname1 = 'NA';
    collabel1 = repmat({'na'}, Nc, 1);
elseif isa(collabel1, 'numeric')
    collabel1 = cellfun(@num2str, num2cell(collabel1), 'UniformOutput', false);
end

if isempty(collabel2)
    collabelname2 = 'NA';
    collabel2 = repmat({'na'}, Nc, 1);
elseif isa(collabel2, 'numeric')
    collabel2 = cellfun(@num2str, num2cell(collabel2), 'UniformOutput', false);
end

if isempty(collabel3)
    collabelname3 = 'NA';
    collabel3 = repmat({'na'}, Nc, 1);
elseif isa(collabel3, 'numeric')
    collabel3 = cellfun(@num2str, num2cell(collabel3), 'UniformOutput', false);
end

if isempty(rowlabel1)
    rowlabelname1 = 'NA';
    rowlabel1 = repmat({'na'}, Nr, 1);
elseif isa(rowlabel1, 'numeric')
    rowlabel1 = cellfun(@num2str, num2cell(rowlabel1), 'UniformOutput', false);
end

if isempty(rowlabel2)
    rowlabelname2 = 'NA';
    rowlabel2 = repmat({'na'}, Nr, 1);
elseif isa(rowlabel2, 'numeric')
    rowlabel2 = cellfun(@num2str, num2cell(rowlabel2), 'UniformOutput', false);
end

if isempty(rowlabel3)
    rowlabelname3 = 'NA';
    rowlabel3 = repmat({'na'}, Nr, 1);
elseif isa(rowlabel3, 'numeric')
    rowlabel3 = cellfun(@num2str, num2cell(rowlabel3), 'UniformOutput', false);
end

fid = fopen([filename '.txt'], 'w');

fprintf(fid, '%s\r\n', strjoin([{'#'} {'#'} {collabelname1} collabel1'], '\t'));
fprintf(fid, '%s\r\n', strjoin([{'#'} {'#'} {collabelname2} collabel2'], '\t'));
fprintf(fid, '%s\r\n', strjoin([{rowlabelname1} {rowlabelname2} {[rowlabelname3 '/' collabelname3]} collabel3'], '\t'));

fmt = strjoin([{'%s\t'} {'%s\t'} {'%s\t'} repmat({'%f\t'}, 1, Nc-1) {'%f\r\n'}], '');

for i = 1:1:Nr
    
    fprintf(fid, fmt, rowlabel1{i}, rowlabel2{i}, rowlabel3{i}, matrix(i,:));
    
end

fclose(fid);


