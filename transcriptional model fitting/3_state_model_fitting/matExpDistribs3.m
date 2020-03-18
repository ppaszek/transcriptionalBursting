function piMat = matExpDistribs3( maxRNA, nLoci, mu, k, delta, times, eps )
%computeDistribs Find the distribution of mRNA number as a function of time.
%   
    % Set the default cutoff for small probabilities
    if( nargin < 7 )
        eps = 0.000001 ;
    end
    
    % Initialise the result
    nTimes = length( times ) ;
    nCellStates = nCellStates3( nLoci ) ;
    nMarkovStates = nCellStates * (maxRNA + 1) ;
    piMat = zeros( nMarkovStates, nTimes ) ;
    
    % Build the rate matrix
    rateMat = buildRateMat3( maxRNA, nLoci, mu, k, delta ) ;
    
    % Set the distribution at t=0. 
    crntPi = zeros( nMarkovStates, 1 ) ;
    crntPi(1) = 1.0 ;
    crntTime = 0.0 ;
     
    % Compute the pi's
    for j = 1:nTimes
        tStep = times(j) - crntTime ;
        [nextPi,s,m,mv,mvd,unA] = expmv( tStep, rateMat, crntPi ) ;
        
        % Set very small entries to zero
        for k = 1:nMarkovStates
            if( nextPi(k) < eps )
                nextPi(k) = 0.0 ;
            end
        end
        
        % Correct the normalisation and
        % record the result
        nextPi = nextPi / sum(nextPi) ;
        piMat(:,j) = nextPi ;
        
        % Get ready for the next iteration.
        crntPi = nextPi ;
        crntTime = times(j) ;
    end
end

