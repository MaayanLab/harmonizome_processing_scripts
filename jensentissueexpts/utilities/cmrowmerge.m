function outcm = cmrowmerge(incm, method)

outcm = struct;

[outcm.term, ~, ri] = unique(incm.term);
outcm.termname = incm.termname;

outcm.numterms = numel(outcm.term);

outcm.entry = incm.entry;
outcm.entryname = incm.entryname;

outcm.numentries = incm.numentries;

if isfield(incm, 'entrydesc')
    outcm.entrydesc = incm.entrydesc;
    outcm.entrydescname = incm.entrydescname;
end

if isfield(incm, 'entryid')
    outcm.entryid = incm.entryid;
    outcm.entryidname = incm.entryidname;
end

outcm.matrix = zeros([outcm.numterms outcm.numentries]);

switch method
    
    case 'mean'
        
        if isfield(incm, 'termdesc') && isfield(incm, 'termid')
            
            outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
            outcm.termdescname = incm.termdescname;
            numtermdescs = zeros([outcm.numterms 1]);
            
            outcm.termid = repmat({'-666'}, outcm.numterms, 1);
            outcm.termidname = incm.termidname;
            numtermids = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termdesc{i} = unique(incm.termdesc(ri==i));
                numtermdescs(i) = numel(outcm.termdesc{i});
                
                outcm.termid{i} = unique(incm.termid(ri==i));
                numtermids(i) = numel(outcm.termid{i});
                
                outcm.matrix(i,:) = mean(incm.matrix(ri==i,:), 1);
            end
            
            if all(numtermdescs == 1)
                outcm.termdesc = vertcat(outcm.termdesc{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), outcm.termdesc, 'UniformOutput', false);
            end
            
            if all(numtermids == 1)
                outcm.termid = vertcat(outcm.termid{:});
            else
                outcm = rmfield(outcm, {'termid' 'termidname'});
            end
            
        elseif isfield(incm, 'termdesc')
            
            outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
            outcm.termdescname = incm.termdescname;
            numtermdescs = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termdesc{i} = unique(incm.termdesc(ri==i));
                numtermdescs(i) = numel(outcm.termdesc{i});
                
                outcm.matrix(i,:) = mean(incm.matrix(ri==i,:), 1);
            end
            
            if all(numtermdescs == 1)
                outcm.termdesc = vertcat(outcm.termdesc{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), outcm.termdesc, 'UniformOutput', false);
            end
            
        elseif isfield(incm, 'termid')
            
            outcm.termid = repmat({'-666'}, outcm.numterms, 1);
            outcm.termidname = incm.termidname;
            numtermids = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termid{i} = unique(incm.termid(ri==i));
                numtermids(i) = numel(outcm.termid{i});
                    
                outcm.matrix(i,:) = mean(incm.matrix(ri==i,:), 1);
            end
            
            if all(numtermids == 1)
                outcm.termid = vertcat(outcm.termid{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), cellfun(@num2str, num2cell(outcm.termid), 'UniformOutput', false), 'UniformOutput', false);
                outcm.termdescname = outcm.termidname;
                outcm = rmfield(outcm, {'termid' 'termidname'});
            end
                
        else
            
            for i = 1:1:outcm.numterms
                outcm.matrix(i,:) = mean(incm.matrix(ri==i,:), 1);
            end
            
        end
        
    case 'nanmean'
        
        if isfield(incm, 'termdesc') && isfield(incm, 'termid')
            
            outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
            outcm.termdescname = incm.termdescname;
            numtermdescs = zeros([outcm.numterms 1]);
            
            outcm.termid = repmat({'-666'}, outcm.numterms, 1);
            outcm.termidname = incm.termidname;
            numtermids = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termdesc{i} = unique(incm.termdesc(ri==i));
                numtermdescs(i) = numel(outcm.termdesc{i});
                
                outcm.termid{i} = unique(incm.termid(ri==i));
                numtermids(i) = numel(outcm.termid{i});
                
                outcm.matrix(i,:) = nanmean(incm.matrix(ri==i,:), 1);
            end
            
            if all(numtermdescs == 1)
                outcm.termdesc = vertcat(outcm.termdesc{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), outcm.termdesc, 'UniformOutput', false);
            end
            
            if all(numtermids == 1)
                outcm.termid = vertcat(outcm.termid{:});
            else
                outcm = rmfield(outcm, {'termid' 'termidname'});
            end
            
        elseif isfield(incm, 'termdesc')
            
            outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
            outcm.termdescname = incm.termdescname;
            numtermdescs = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termdesc{i} = unique(incm.termdesc(ri==i));
                numtermdescs(i) = numel(outcm.termdesc{i});
                
                outcm.matrix(i,:) = nanmean(incm.matrix(ri==i,:), 1);
            end
            
            if all(numtermdescs == 1)
                outcm.termdesc = vertcat(outcm.termdesc{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), outcm.termdesc, 'UniformOutput', false);
            end
            
        elseif isfield(incm, 'termid')
            
            outcm.termid = repmat({'-666'}, outcm.numterms, 1);
            outcm.termidname = incm.termidname;
            numtermids = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termid{i} = unique(incm.termid(ri==i));
                numtermids(i) = numel(outcm.termid{i});
                    
                outcm.matrix(i,:) = nanmean(incm.matrix(ri==i,:), 1);
            end
            
            if all(numtermids == 1)
                outcm.termid = vertcat(outcm.termid{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), cellfun(@num2str, num2cell(outcm.termid), 'UniformOutput', false), 'UniformOutput', false);
                outcm.termdescname = outcm.termidname;
                outcm = rmfield(outcm, {'termid' 'termidname'});
            end
                
        else
            
            for i = 1:1:outcm.numterms
                outcm.matrix(i,:) = nanmean(incm.matrix(ri==i,:), 1);
            end
            
        end
        
    case 'max'
        
        if isfield(incm, 'termdesc') && isfield(incm, 'termid')
            
            outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
            outcm.termdescname = incm.termdescname;
            numtermdescs = zeros([outcm.numterms 1]);
            
            outcm.termid = repmat({'-666'}, outcm.numterms, 1);
            outcm.termidname = incm.termidname;
            numtermids = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termdesc{i} = unique(incm.termdesc(ri==i));
                numtermdescs(i) = numel(outcm.termdesc{i});
                
                outcm.termid{i} = unique(incm.termid(ri==i));
                numtermids(i) = numel(outcm.termid{i});
                
                outcm.matrix(i,:) = max(incm.matrix(ri==i,:), [], 1);
            end
            
            if all(numtermdescs == 1)
                outcm.termdesc = vertcat(outcm.termdesc{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), outcm.termdesc, 'UniformOutput', false);
            end
            
            if all(numtermids == 1)
                outcm.termid = vertcat(outcm.termid{:});
            else
                outcm = rmfield(outcm, {'termid' 'termidname'});
            end
            
        elseif isfield(incm, 'termdesc')
            
            outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
            outcm.termdescname = incm.termdescname;
            numtermdescs = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termdesc{i} = unique(incm.termdesc(ri==i));
                numtermdescs(i) = numel(outcm.termdesc{i});
                
                outcm.matrix(i,:) = max(incm.matrix(ri==i,:), [], 1);
            end
            
            if all(numtermdescs == 1)
                outcm.termdesc = vertcat(outcm.termdesc{:});
            else
                outcm.termdesc = cellfun(@(x) strjoin(x', '+'), outcm.termdesc, 'UniformOutput', false);
            end
            
        elseif isfield(incm, 'termid')
            
            outcm.termid = repmat({'-666'}, outcm.numterms, 1);
            outcm.termdesc = repmat({'-666'}, outcm.numterms, 1);
            outcm.termidname = incm.termidname;
            outcm.termdescname = incm.termidname;
            numtermids = zeros([outcm.numterms 1]);
            
            for i = 1:1:outcm.numterms
                outcm.termid{i} = unique(incm.termid(ri==i));
                outcm.termdesc{i} = strjoin(cellfun(@num2str, num2cell(outcm.termid{i}), 'UniformOutput', false)', '+');
                numtermids(i) = numel(outcm.termid{i});
                    
                outcm.matrix(i,:) = max(incm.matrix(ri==i,:), [], 1);
            end
            
            if all(numtermids == 1)
                outcm.termid = vertcat(outcm.termid{:});
                outcm = rmfield(outcm, {'termdesc' 'termdescname'});
            else
                outcm = rmfield(outcm, {'termid' 'termidname'});
            end
                
        else
            
            for i = 1:1:outcm.numterms
                outcm.matrix(i,:) = max(incm.matrix(ri==i,:), [], 1);
            end
            
        end
        
    otherwise
        
        warning('Invalid method.');
        
end


