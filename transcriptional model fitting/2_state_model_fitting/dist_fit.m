function aa =dist_fit(z)

    % Set parameters
    maxRNA = 2000;  % Maximum RNA count
    nLoci = 2;      % total number of sites of transcription
    
    %Indivisual paramnetr valueas can be either fitted or fixed
  
    mu_0 = 0 ;      % mRNA synthesis rate if promoter is inactive
    k_1 = z(1) ;    % promoter activation rate
    k_0 = z(2) ;    % promoter deactivation rate
    mu_1 = z(3) ;     % mRNA synthesis rate if promoter is active
    delta =  z(4); % mRNA degradation rate 
    
     

    k_vec = [k_0, k_1] ;
    mu_vec = [mu_0, mu_1] ;

    % Specify a list of times at which we make observations and then,
    % for each time, do nCells Gillespie simulations and work out the 
    % number of mRNA at the desired time.
    timeStep = 5 ;
    observationTimes = [180];     % should be multiples of timeStep



    mePiMat = matExpDistribs( maxRNA, nLoci, mu_vec, k_vec, delta, observationTimes ) ;

    a=[1:3:3*(maxRNA+1)];
    b=[2:3:3*(maxRNA+1)];
    c=[3:3:3*(maxRNA+1)];

    pi=mePiMat(a,:)+mePiMat(b,:)+mePiMat(c,:);
    pii=cumsum(pi);

    
    xx=[[0:1:10],[20:10:90],[100:50:950],[1000:100:1500]];
    
   
    %provide cdf to be fitted to (using checking_dist_tnf code)
    
    %Objective function: minimize the distance between data cdf and
    %simulated cdf
    p1=[0.0307453112911650,0.0398561142589718,0.0502924135186276,0.0621825856481706,0.0755285013653620,0.0902576711360412,0.106230903392499,0.123292355398006,0.141286311084427,0.160069535432610,0.179462460982657,0.375240686651282,0.561936432879814,0.709802415338990,0.805555503562317,0.866071239241690,0.898627054261813,0.929838851103308,0.953614849895031,0.964200902374883,0.988798536807525,0.995833333333331,0.998469720813099,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998,0.999999999999998];
    aa=sum(abs(p1-p2'))/length(p1);

end