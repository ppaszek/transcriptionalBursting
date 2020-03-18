% Time the function that computes distributions
clear % wipe the slate clean

% Set parameters
maxRNA = 2000;  % Maximum RNA count
nLoci = 2;      % total number of sites of transcription
mu_0 = 0 ;      % mRNA synthesis rate if promoter is inactive
mu_1 = 8.3 ;     % mRNA synthesis rate if promoter is active
k_0 = 0.04 ;    % promoter deactivation rate
k_1 = 0.02 ;    % promoter activation rate
delta = 0.017;   % mRNA degradation rate 

k_vec = [k_0, k_1] ;
mu_vec = [mu_0, mu_1] ;

% Specify a list of times at which we make observations and then,
% for each time, do nCells Gillespie simulations and work out the 
% number of mRNA at the desired time.
timeStep = 5 ;
observationTimes = [180];     % should be multiples of timeStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Time the approach that uses matrix exponentials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure the function expmv() can be found.
addpath('HighamAlMohly')

% The main event
tic
nReps = 1 ;
for j = 1:nReps
    mePiMat = matExpDistribs( maxRNA, nLoci, mu_vec, k_vec, delta, observationTimes ) ;
end

% Report how long it all took
elapsed = toc 
timePerRun = elapsed / nReps 

a=[1:3:3*(maxRNA+1)];
b=[2:3:3*(maxRNA+1)];
c=[3:3:3*(maxRNA+1)];

pi=mePiMat(a,:)+mePiMat(b,:)+mePiMat(c,:);
figure(1);plot(pi);
figure(2);plot(cumsum(pi))


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   Time Gillespie simulations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nCells = 2000 ;
% 
% % The main event
% tic
% nReps = 50 ;
% for j = 1:nReps
%     gPiMat = gillespieDistribs( maxRNA, nLoci, nCells, mu_vec, k_vec, delta, observationTimes ) ;
% end
% 
% % Report how long it all took
% elapsed = toc 
% timePerRun = elapsed / nReps 

