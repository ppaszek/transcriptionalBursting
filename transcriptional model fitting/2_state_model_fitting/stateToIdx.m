function idx = stateToIdx( s, m, nLoci )
% stateToIdx get the index for a state
%   s is the number of active sites of transcription
%   m is the number of transcripts
    idx = s + (m * (nLoci + 1)) + 1 ;
end