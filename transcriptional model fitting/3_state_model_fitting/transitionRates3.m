function [rateVec, rateSum] = transitionRates3( nc, ni, na, m, mu, k, delta )
%transitionRates3() Conpute transition rates *from* a given state

% Unpack the args
mu_0 = mu(1) ;  % Translation for the inactive state
mu_1 = mu(2) ;  % Translation form the active state

k_0 = k(1) ;    % Active -> Inactive
k_1 = k(2) ;    % Inactive -> Active
k_2 = k(3) ;    % Inactive -> Closed
k_3 = k(4) ;    % Closed -> Inactive
k_4 = k(5) ;    % Active -> Closed 

rateVec = zeros( 7, 1 ) ;

rateVec(1) = na * k_0 ;   % Active -> Inactive
rateVec(2) = ni * k_1 ;   % Inactive -> Active
rateVec(3) = ni * k_2 ;   % Inactive -> Closed
rateVec(4) = nc * k_3 ;   % Closed -> Inactive
rateVec(5) = na * k_4 ;   % Active -> Closed
rateVec(6) = (na*mu_1) + (ni*mu_0) ;  % Transcription
rateVec(7) = m * delta ;              % mRNA Degradation

rateSum = sum( rateVec ) ;

end

