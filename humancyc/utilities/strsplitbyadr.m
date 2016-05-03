function outcells = strsplitbyadr(instring, delimiter)

if strcmp(delimiter, '|') || strcmp(delimiter, '.')
    delimiter = ['[' delimiter ']'];
end

si = regexp(instring, delimiter);

N = numel(si);

if N == 0
    
    outcells = {instring};
    
else

    outcells = cell(1, N+1);

    outcells{1} = instring(1:si(1)-1);

    for i = 2:1:N
        outcells{i} = instring(si(i-1)+1:si(i)-1);
    end

    outcells{end} = instring(si(end)+1:end);

end


