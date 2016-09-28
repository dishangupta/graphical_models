%%Sample from Z~Bernoulli(p)
function [ z ] = Bernoulli( p )
    
    z = 0;
    
    Z = [p 1-p];
    r = rand(1,1);
    
    cumsum = 0;
    for i = 1:size(Z,2)
        cumsum = cumsum + Z(i);
        if (r < cumsum)
            z = i;
            break;
        end
    end
    
end

