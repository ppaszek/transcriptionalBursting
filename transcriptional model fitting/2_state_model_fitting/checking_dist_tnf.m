clc

[num,txt,raw] = xlsread('tnf_smFISH.xlsx');

%  select specific column in the excel spreadsheet to calculate cdf
x=num(:,6);
r=find(isfinite(x)==1);
data=x(r);


% Plot empirical cdf from the data. 

[f,x,flo,fup]=ecdf(data,'alpha',0.05,'bounds','on')
figure(1)
plot(x,f,'b');
hold on
plot(x, flo,'-b')
plot(x, fup,'-b')


%Use kernel density estimator to approximate cdf at specific points
%p1 is the input for the objective function of the ga
% bandwitch parametr can  be adjected to obtain accurate
% approximation
xx=[[0:1:10],[20:10:90],[100:50:950],[1000:100:1500]];
figure(1)
pd1 = fitdist(data,'Kernel','Kernel','epanechnikov','bandwidth',5);
p1=cdf(pd1,xx);
plot(xx,p1,'o')
hold on
 
