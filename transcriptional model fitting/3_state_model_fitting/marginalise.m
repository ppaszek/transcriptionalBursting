function txCountProbMat = marginalise( stateProbMat, nLoci )
%marginalise() Get probs of mRNA counts

nModelStates = size(stateProbMat, 1) ;
nTimes = size(stateProbMat, 2) ;
nCellStates = nCellStates3( nLoci ) ;

maxRNA = (nModelStates/nCellStates) - 1 ;

txCountProbMat = zeros( maxRNA+1, nTimes ) ;

row = 1 ;
for j=1:(maxRNA + 1)
    for s=1:nCellStates
        for k=1:nTimes
            txCountProbMat(j,k) = txCountProbMat(j,k) + stateProbMat( row, k ) ;
        end
        
        % Move on to the next row of the original matrix
        row = row + 1 ;
    end
end

end

