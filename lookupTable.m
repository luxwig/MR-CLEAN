function V = lookupTable(keys,values,K)
    if sum(keys==K,'all') ~=1
        error('More than one key or no keys have been found for %f', K)
    else
        V = values(keys==K);
    end
end

