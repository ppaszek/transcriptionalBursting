%load(gaoptions.mat)

%options = gaoptimset(@ga)
% options1=gaoptimset('PopulationSize',200)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Time the approach that uses matrix exponentials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure the function expmv() can be found.
addpath('HighamAlMohly')
load('gaoptions.mat')
starttime=clock;

for i=1:100
i


lb=[0;0;0;0;0;0;0.002];   %lower and upper parameter bounds for kon, koff, kt and kd, respectively
ub=[0.2;0.2;0.2;0.2;20;20;0.006];

[x, fval, exitflag, output, scores]=ga(@state3_IL1,7,[],[],[],[],lb,ub,[],options_gui);

DATA(i).par=x;
DATA(i).fun=fval;
DATA(i).exit=exitflag;
DATA(i).out=output;

end;

etime(clock,starttime)
save IL1_model_fit

