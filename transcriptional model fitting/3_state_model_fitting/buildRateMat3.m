function rateMat = buildRateMat3( maxRNA, nLoci, mu, k, delta )
%buildRateMat: transition rate matrix
%   k = [k_0, k_1, k_2, k_3, k_4]
%   mu = [mu_0, mu_1]

    % Prepare the result
    nCellStates = nCellStates3( nLoci ) ;
    nMarkovStates = nCellStates * (maxRNA + 1) ;
    rateMat = zeros( nMarkovStates, nMarkovStates ) ;

    % Loop over all states of the Markov chain, adding entries to 
    % the rate matrix to describe all possible *outgoing* reactions.
    for m = 0:maxRNA % number of mRNA
        for na = 0:nLoci % number of active sites
            for ni = 0:(nLoci-na) % number of inactive, but not closed sites
                fromIdx = stateToIdx3( ni, na, m, nLoci ) ;
                nc = nLoci - (ni + na) ;
                
                rateVec = transitionRates3( nc, ni, na, m, mu, k, delta ) ;

                % Active -> Inactive
                if( na > 0 ) 
                    toIdx = stateToIdx3( ni+1, na-1, m, nLoci ) ;
                    rateMat(toIdx, fromIdx) = rateVec(1) ;
                    rateMat(fromIdx, fromIdx) = rateMat(fromIdx, fromIdx) - rateVec(1) ;
                end
                
                % Inactive -> Active
                if( ni > 0 ) 
                    toIdx = stateToIdx3( ni-1, na+1, m, nLoci ) ;
                    rateMat(toIdx, fromIdx) = rateVec(2) ;
                    rateMat(fromIdx, fromIdx) = rateMat(fromIdx, fromIdx) - rateVec(2) ;
                end
                
                % Inactive -> Closed
                if( ni > 0 ) 
                    toIdx = stateToIdx3( ni-1, na, m, nLoci ) ;
                    rateMat(toIdx, fromIdx) = rateVec(3) ;
                    rateMat(fromIdx, fromIdx) = rateMat(fromIdx, fromIdx) - rateVec(3) ;
                end

                % Closed -> Inactive
                if( nc > 0 ) 
                    toIdx = stateToIdx3( ni+1, na, m, nLoci ) ;
                    rateMat(toIdx, fromIdx) = rateVec(4) ;
                    rateMat(fromIdx, fromIdx) = rateMat(fromIdx, fromIdx) - rateVec(4) ;
                end

                 % Active -> Closed
                if( na > 0 ) 
                    toIdx = stateToIdx3( ni, na-1, m, nLoci ) ;
                    rateMat(toIdx, fromIdx) = rateVec(5) ;
                    rateMat(fromIdx, fromIdx) = rateMat(fromIdx, fromIdx) - rateVec(5) ;
                end

                % Transcription
                rateMat(fromIdx, fromIdx) = rateMat(fromIdx, fromIdx) - rateVec(6) ;
                if( m < maxRNA )
                    toIdx = stateToIdx3( ni, na, m+1, nLoci ) ;
                    rateMat(toIdx, fromIdx) = rateVec(6) ;
                end

                % mRNA degradation
                if( m > 0 )
                    toIdx = stateToIdx3( ni, na, m-1, nLoci ) ;
                    rateMat(toIdx, fromIdx) = rateVec(7) ;
                    rateMat(fromIdx, fromIdx) = rateMat(fromIdx, fromIdx) - rateVec(7) ;
                end
            end
        end
    end
    
    % Sparsify the matrix
    rateMat = sparse( rateMat ) ;
end