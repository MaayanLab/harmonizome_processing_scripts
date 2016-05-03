function edge = sourcetargetcat(source, target)

edge = [source target];
edge = mat2cell(edge, ones([numel(source) 1]));
edge = cellfun(@(x) [x{1} '+' x{2}], edge, 'UniformOutput', false);


