% Test the marginalisation code
function aa =state3_IL1(z)% Set parameters: all rates are in 1/min
    
    nLoci = 2 ;     % total number of sites of transcription
    maxRNA = 2000 ;  % Maximum RNA count

    
 

    k_0 = z(1);       % active -> inactive  /koff
    k_1 = z(2);     % inactive -> active   /kon
    k_2 = z(3);    % inactive -> closed /toff
    k_3 = z(4);      % closed -> inactive  /ton
    k_4 = 0;     % active -> closed  /toff (k2=k4) Additional regulation can be added to the model
    
    mu_0 = z(5);      % mRNA synthesis rate if promoter is inactive
    mu_1 = z(6);    % mRNA synthesis rate if promoter is active

    delta =0.0037;  % mRNA degradation rate

    k = [k_0, k_1, k_2, k_3, k_4] ;
    mu = [mu_0, mu_1] ;

% Specify a list of times at which we make observations and then,
% for each time, do nCells Gillespie simulations and work out the 
% number of mRNA at the desired time.
    observationTimes = [3*180]; 
    nTimes = length(observationTimes) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Do the matrix exponential version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure the function expmv() can be found.
%addpath('HighamAlMohly')

% Do the computation
    stateProbMat = matExpDistribs3( maxRNA, nLoci, mu, k, delta, observationTimes ) ;
    pi=marginalise( stateProbMat, nLoci );
    txCountProbMat = cumsum(pi) ;

% Check whether the columns sum to 1.0
    %colSums = sum( txCountProbMat ) 
    
    xx=[[0:1:10],[20:10:90],[100:50:950],[1000:100:1500]];
    %provide cdf to be fitted to (using checking_dist_tnf code)
    
    %Objective function: minimize the distance between data cdf and
    %simulated cdf
    
    %control dataset
    p1=[0.460882100133027,0.526428024091774,0.591536360748628,0.655155476132467,0.716217338440335,0.773684285462079,0.826500211090624,0.873591802696770,0.913885747651316,0.946308733325063,0.969787447088811,0.998881030652367,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998];


    xxx=xx+1;
    p2=txCountProbMat(xxx); % extract simulated cdf
    
    aa=sum(abs(p1-p2'))/length(p1);
end




