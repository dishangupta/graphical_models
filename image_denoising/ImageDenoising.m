load('hw1_images.mat');

% Noisy Observed Image
X = noisyImg;

% Original Image 
Zorig = origImg;

% Initialize Z to noisy values
Z = randi([0 1], size(X)); 
Z(find(Z==0)) = -1;

% Energy, gradient and error initialization 
oldE = inf;
E = 0;
deltaEpos = 0;
deltaEneg = 0;
error = 0;

% Parameter initialization
h = -0.75; 
beta = 5;
nu = 7;
maxIter = 30;
epsilon = 0.001;

for k = 1:maxIter   
    
    if (oldE - E < epsilon)
        break;
    end
    
    oldE = E;
    
    % Traverse each Zij
    for i = 1:size(Z, 1)
        for j = 1:size(Z,2)

            Zpos = 1;
            Zneg = -1;

            deltaEpos = h*Zpos - nu*(Zpos*X(i,j));
            deltaEneg = h*Zneg - nu*(Zneg*X(i,j));

            % Corner Case 1
            if (i==1 && j==1) 
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i,j+1) + Z(i+1,j)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i,j+1) + Z(i+1,j)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            %Corner Case 2    
            elseif (i==1 && j==size(Z,2))
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i+1,j) + Z(i,j-1)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i+1,j) + Z(i,j-1)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            %Corner Case 3    
            elseif (i==size(Z,1) && j==1)
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i,j+1) + Z(i-1,j)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i,j+1) + Z(i-1,j)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            %Corner Case 4    
            elseif (i==size(Z,1) && j==size(Z,2))
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i-1,j) + Z(i,j-1)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i-1,j) + Z(i,j-1)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            % Edge Case 1
            elseif (i==1)
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i,j+1) + Z(i,j-1) + Z(i+1,j)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i,j+1) + Z(i,j-1) + Z(i+1,j)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            % Edge Case 2    
            elseif (j==1)
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i,j+1) + Z(i-1,j) + Z(i+1,j)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i,j+1) + Z(i-1,j) + Z(i+1,j)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            % Edge Case 3
            elseif (i==size(Z,1))
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i,j+1) + Z(i,j-1) + Z(i-1,j)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i,j+1) + Z(i,j-1) + Z(i-1,j)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            % Edge Case 4    
            elseif (j==size(Z,2))
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i-1,j) + Z(i+1,j) + Z(i,j-1)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i-1,j) + Z(i+1,j) + Z(i,j-1)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end

            % Normal Case
            else
                deltaEpos = deltaEpos - beta*(Zpos*(Z(i,j+1) + Z(i-1,j) + Z(i,j-1) + Z(i+1,j)));
                deltaEneg = deltaEneg - beta*(Zneg*(Z(i,j+1) + Z(i-1,j) + Z(i,j-1) + Z(i+1,j)));
                if(deltaEpos < deltaEneg)
                    Z(i,j) = 1;
                else
                    Z(i,j) = -1;
                end
            end
        end
    end
    
    % Compute energy
    temp = Z.*X;
    E = h*sum(Z(:)) - nu*(sum(temp(:))); 
        
    % Traverse each Zij again
    for i = 1:size(Z, 1)
        for j = 1:size(Z,2)

            % Corner Case 1
            if (i==1 && j==1) 
               E = E - (beta/2)*(Z(i,j)*(Z(i,j+1) + Z(i+1,j)));

            %Corner Case 2    
            elseif (i==1 && j==size(Z,2))
                E = E - (beta/2)*(Z(i,j)*(Z(i+1,j) + Z(i,j-1)));
                            
            %Corner Case 3    
            elseif (i==size(Z,1) && j==1)
                E = E - (beta/2)*(Z(i,j)*(Z(i,j+1) + Z(i-1,j)));
            
            %Corner Case 4    
            elseif (i==size(Z,1) && j==size(Z,2))
                E = E - (beta/2)*(Z(i,j)*(Z(i-1,j) + Z(i,j-1)));
                
            % Edge Case 1
            elseif (i==1)
                E = E - (beta/2)*(Z(i,j)*(Z(i,j+1) + Z(i,j-1) + Z(i+1,j)));
               
            % Edge Case 2    
            elseif (j==1)
                E = E - (beta/2)*(Z(i,j)*(Z(i,j+1) + Z(i-1,j) + Z(i+1,j)));
                
            % Edge Case 3
            elseif (i==size(Z,1))
                E = E - (beta/2)*(Z(i,j)*(Z(i,j+1) + Z(i,j-1) + Z(i-1,j)));
            
            % Edge Case 4    
            elseif (j==size(Z,2))
                E = E - (beta/2)*(Z(i,j)*(Z(i-1,j) + Z(i+1,j) + Z(i,j-1)));
            
            % Normal Case
            else
                E = E - (beta/2)*(Z(i,j)*(Z(i,j+1) + Z(i-1,j) + Z(i,j-1) + Z(i+1,j)));
            end
        end
    end
end

error = 100*numel(find((Zorig - Z)~=0))/(size(Z,1)*size(Z,2));
            