function rateMat = buildRateMat( maxRNA, nLoci, mu, k, delta )
%buildRateMat: transition rate matrix
%   k = [k_0, k_1]
%   mu = [mu_0, mu_1]

    % Prepare the result
    nStates = (nLoci + 1)*(maxRNA + 1) ;
    rateMat = zeros( nStates, nStates ) ;

    % Loop over all states, considering all possible
    % *outgoing* reactions.
    for m = 0:maxRNA % number of mRNA
        for s = 0:nLoci % number of active sites
            myIdx = stateToIdx( s, m, nLoci ) ;
            % Activation
            if( s < 2 ) 
                resultIdx = stateToIdx( s+1, m, nLoci ) ;
                activationRate = (2 - s) * k(2) ; % (2 - s) * k_1
                rateMat(resultIdx, myIdx) = activationRate ;
                rateMat(myIdx, myIdx) = rateMat(myIdx, myIdx) - activationRate ;
            end

            % Deactivation
            if( s > 0 ) 
                resultIdx = stateToIdx( s-1, m, nLoci ) ;
                deactivationRate = s * k(1) ; 
                rateMat(resultIdx, myIdx) = deactivationRate ;
                rateMat(myIdx, myIdx) = rateMat(myIdx, myIdx) - deactivationRate ;
            end

            % Transcription
            transcriptionRate = (2 - s) * mu(1) + s * mu(2) ;
            rateMat(myIdx, myIdx) = rateMat(myIdx, myIdx) - transcriptionRate ;
            if( m < maxRNA )
                resultIdx = stateToIdx( s, m+1, nLoci ) ;
                rateMat(resultIdx, myIdx) = transcriptionRate ;
            end

            % mRNA degradation
            if( m > 0 )
                resultIdx = stateToIdx( s, m-1, nLoci ) ;
                degradationRate = m * delta ; 
                rateMat(resultIdx, myIdx) = degradationRate ;
                rateMat(myIdx, myIdx) = rateMat(myIdx, myIdx) - degradationRate ;
            end
        end
    end
    
    % Sparsify the matrix
    rateMat = sparse( rateMat ) ;
end