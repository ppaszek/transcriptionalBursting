

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Time the approach that uses matrix exponentials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure the function expmv() can be found.

addpath('HighamAlMohly')
load('gaoptions.mat') % this file contains paremetrs of the ga
starttime=clock;

for i=1:100
i


lb=[0;0;0;0.006];   %lower and upper parameter bounds for kon, koff, kt and kd, respectively
ub=[0.2;0.2;20;0.07];

[x, fval, exitflag, output, scores]=ga(@dist_fit,4,[],[],[],[],lb,ub,[],options_gui);

DATA(i).par=x;
DATA(i).fun=fval;
DATA(i).exit=exitflag;
DATA(i).out=output;

end;

etime(clock,starttime)
save TNF_model_fit

