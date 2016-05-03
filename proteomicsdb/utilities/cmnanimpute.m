function cm = cmnanimpute(cm, Method, Dim)

switch Method
    
    case 'mean'
        
        switch Dim

            case 'row'

                rowmean = nanmean(cm.matrix, 2);

                missmat = isnan(cm.matrix);

                cm.matrix(missmat) = 0;

                cm.matrix = cm.matrix + missmat.*repmat(rowmean, 1, cm.numentries);

            case 'column'

                colmean = nanmean(cm.matrix, 1);

                missmat = isnan(cm.matrix);

                cm.matrix(missmat) = 0;

                cm.matrix = cm.matrix + missmat.*repmat(colmean, cm.numterms, 1);

            otherwise

                disp('Dim must be row or column.');

        end
        
    case 'median'
        
        switch Dim

            case 'row'

                rowmedian = nanmedian(cm.matrix, 2);

                missmat = isnan(cm.matrix);

                cm.matrix(missmat) = 0;

                cm.matrix = cm.matrix + missmat.*repmat(rowmedian, 1, cm.numentries);

            case 'column'

                colmedian = nanmedian(cm.matrix, 1);

                missmat = isnan(cm.matrix);

                cm.matrix(missmat) = 0;

                cm.matrix = cm.matrix + missmat.*repmat(colmedian, cm.numterms, 1);

            otherwise

                disp('Dim must be row or column.');

        end
        
    otherwise
        
        disp('Method must be mean or median.');
        
end


