function idx = stateToIdx3( ni, na, m, nLoci )
% stateToIdx get the index for a state
%	nc is the number of loci that are in closed chromatin
%   ni is the number of inactive sites of transcription
%   m is the number of transcripts

nStates = nCellStates3( nLoci ) ;
idx =  ni + na*(nLoci + 1) - (na * (na - 1)/2) + (m * nStates) + 1 ;

end