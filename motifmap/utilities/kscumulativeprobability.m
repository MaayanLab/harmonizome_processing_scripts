function P = kscumulativeprobability(X, method, support)

if exist('support', 'var') == 0 || isempty(support)
    support = 'unbounded';
end

switch method
    
    case 'cols'
        
        [numobs, numvect] = size(X);
        
        P = zeros([numobs numvect]);
        
        noexistingpool = isempty(gcp('nocreate'));
        
        if numobs <= 1000
        
            parfor i = 1:numvect

                P(:,i) = ksdensity(X(:,i), X(:,i), 'function', 'cdf', 'support', support);

            end
        
        else
            
            parfor i = 1:numvect

                [p, x] = ksdensity(X(:,i), 'npoints', 1000, 'function', 'cdf', 'support', support);
                
                P(:,i) = interp1(x, p, X(:,i));

            end
            
        end
        
        if noexistingpool
            delete(gcp);
        end
        
    case 'rows'
        
        X = X';
        
        [numobs, numvect] = size(X);
        
        P = zeros([numobs numvect]);
        
        noexistingpool = isempty(gcp('nocreate'));
        
        if numobs <= 1000
        
            parfor i = 1:numvect

                P(:,i) = ksdensity(X(:,i), X(:,i), 'function', 'cdf', 'support', support);

            end
        
        else
            
            parfor i = 1:numvect

                [p, x] = ksdensity(X(:,i), 'npoints', 1000, 'function', 'cdf', 'support', support);
                
                P(:,i) = interp1(x, p, X(:,i));

            end
            
        end
        
        if noexistingpool
            delete(gcp);
        end
        
        P = P';
        
    case 'matrix'
        
        [numobs, numvect] = size(X);
        
        if numobs*numvect <= 3000
            
            P = ksdensity(X(:), X(:), 'function', 'cdf', 'support', support);
            
            P = reshape(P, numobs, numvect);
            
        else
        
            [p, x] = ksdensity(X(:), 'npoints', 3000, 'function', 'cdf', 'support', support);
            
            P = interp1(x, p, X(:));
            
            P = reshape(P, numobs, numvect);
            
        end
        
    otherwise
        
        disp('invalid method. input ''cols'' ''rows'' or ''matrix''.');
        
end


