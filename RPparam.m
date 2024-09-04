function [tau, dim, ps, index] = RPparam(signal, mFlag)


                % DELAY (tau): Average Mutual Information
                    mi = mutual(signal); %mutual information index
                    [~, ~, ~, imin] = extrema(mi);
                    tau = imin(end); %pick the LAST value (minima sorted descending)                     
%                     figure(300),clf %plot to double-check
%                     display(['Sorted tau: ' num2str(imin)])
%                     display(['Selected tau: ' num2str(tau)])
%                     plot(mi , 'ks')
%                     xlabel('Tau') , ylabel('Mutual information index')
%                     
                % DIMENSION: False Nearest Neighbors
                    if mFlag == 1
                        maxdim = 10; % set high dimensionality for embedding, to further trim down
                        out = false_nearest(signal, 1, maxdim, tau);
                        fnn = out(:,2);
                        dim = dsearchn(fnn , .1); %number of dimensions with fnn closest to the second argument
                    else
                        dim=[];
                    end
                    display(['M = ' num2str(dim)]);
                    
                % DISTANCE THRESHOLD: 
                    % Get phase space size per trial (set at 10% of phase space)
                        ps=pss(signal, dim, tau, 'maxnorm')*.10;
                        
                    % Old "Mattia"-solution (Fixed amount of nearest neighbours)
                        % fan = length(signal)*.05;  
                    
                index=1;
                    
end
