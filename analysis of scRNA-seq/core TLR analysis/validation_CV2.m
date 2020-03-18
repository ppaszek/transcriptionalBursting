[num1,txt1,raw1]=xlsread('summary_mean_var_812_no_zeros_sorted.xls');

[num2,txt2,raw2]=xlsread('induced_genes_outliers_removed_812_no_zeros.xls');

genelist=txt1(2:end,1);
singlenames=txt2(1:end,2);
condlist=txt2(1:end,3);

Pert=['On_Chip_Stimulation_LPS_4h', 'Ifnar1_KO_LPS_4h','IFNB_2h','Stat1_KO_LPS_4h','Tnfr_KO_LPS_4h','gpc_t4','hero_t0','ifnr_t275','gpa_t4','gpb_t4','tr1_t4','tr2_t4','tr3_t2','tr4_t2','wps_t4'];

stim=['Unstimulated','Unstimulated_Replicate','LPS_1h','LPS_2h','LPS_4h','LPS_6h','LPS_4h_Replicate','PAM_1h','PAM_2h','PAM_4h','PAM_6h','PIC_1h','PIC_2h','PIC_4h','PIC_6h'];


F = @(p,xx) p(1)+p(2)*log(xx);
x0=[1000,1];

figure(1)
for i=1:150
   
    
   h=find(strcmp(genelist(i),singlenames)>0);
   x=num2(h,5);
   y=num2(h,11);
   
   subplot(10,15,i)
   scatter(x,y,'k')
   hold on
   title(genelist(i))
   
   [X,hi]=sort(x);
   Y=y(hi);
   X=X';
   Y=Y';
   [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,X,Y);
   plot(X,F(p,X),'b')
   %[p,s] = polyfit(x,y,1);
   %plot(x,p(2)+p(1)*x,'b')
   
   %R2 = corr(y,p(2)+p(1)*x);
   %R2= 1 - (s.normr/norm(y - mean(y)))^2;
    R2=1-sum((Y-F(p,X)).^2)/(length(Y)*var(Y)); 
   %tx=['R^2 ',num2str(round(R2,2))];
   tx=['R^2 ',num2str(round(R2,2))];
   text(10,max(y)-max(y)/5, tx);
   
   xlim([0,max(x)])
   ylim([0,max(y)])
   box on
end;

figure(2)
for i=1:140
   
    
   h=find(strcmp(genelist(i+150),singlenames)>0);
   x=num2(h,5);
   y=num2(h,11);
   
   subplot(10,15,i)
   scatter(x,y,'k')
   title(genelist(i+150))
   hold on
   %[p,s] = polyfit(x,y,1);
   %plot(x,p(2)+p(1)*x,'b')
   
   %R2 = corr(y,p(2)+p(1)*x);
   %R2= 1 - (s.normr/norm(y - mean(y)))^2;
   [X,hi]=sort(x);
   Y=y(hi);
   X=X';
   Y=Y';
   [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,X,Y);
   plot(X,F(p,X),'b')
   %[p,s] = polyfit(x,y,1);
   %plot(x,p(2)+p(1)*x,'b')
   
   %R2 = corr(y,p(2)+p(1)*x);
   %R2= 1 - (s.normr/norm(y - mean(y)))^2;
    R2=1-sum((Y-F(p,X)).^2)/(length(Y)*var(Y)); 
   tx=['R^2 ',num2str(round(R2,2))];
   text(10,max(y)-max(y)/5, tx);
   
   xlim([0,max(x)])
   ylim([0,max(y)])
   box on
   
end;

% figure(3)
% for i=1:84
%    
%     
%    h=find(strcmp(genelist(i+200),singlenames)>0);
%    x=num2(h,5);
%    y=num2(h,6);
%    
%    subplot(10,10,i)
%    scatter(x,y)
%    title(genelist(i+200))
%    hold on
%    [p,s] = polyfit(x,y,1);
%    plot(x,p(2)+p(1)*x,'b')
%    
%    %R2 = corr(y,p(2)+p(1)*x);
%    R2= 1 - (s.normr/norm(y - mean(y)))^2;
%    tx=['R^2 ',num2str(round(R2,2))];
%    text(10,max(y)-max(y)/10, tx);
%    
%    xlim([0,max(x)])
%    ylim([0,max(y)])
%    box on
% end;
% 
% 
% N=204
% genes=genelist(1:204);
% 
% Min=min(num1(1:N,4));
% Max=max(num1(1:N,5));
% MinMax=min(num1(1:N,5));
% int=num1(1:N,2);
% slope=num1(1:N,3);
%     
% xx=[[0:5:100],[120:20:500],[550:50:1000],[1000:250:30000]];
% 
% 
% AA=NaN*zeros(length(xx),N);
% 
% for i=1:N
%     h1=min(find(xx>num1(i,4)));
%     h2=min(find(xx>num1(i,5)));
% 
%     for j=h1:h2
%         AA(j,i)=int(i)+slope(i)*xx(j);
%     end;
% end;


%[C,ia,ic] = unique(txt1);

% nn=2;
% kk='b'
% 
% range=Inf;
% 
% LPS=[23,3,4,5,7];
% PAM=[23,9,10,11,12];
% PIC=[23,13,14,15,16];
% all=[23,3,4,5,7,9,10,11,12,13,14,15,16];
% Pert=[5,1,2,8,17,18];
% 
% All_rel=[1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,23,24];
% 
% 
% 
% P1=zeros(143,5);
% 
% xx=[[1:10:100],[200:50:1000],[1250:250:30000]];
% 
% S1h=[23,3,9,13];
% S2h=[23,4,10,14];
% S4h=[23,5,11,15];
% S6h=[23,7,12,16];


% figure(10)
% for i=1:length(Pert);
%     h=find(ic==Pert(i));
%     data=num1(h);
%     [f,x] = ecdf(data);
%     plot(x,f,'k')
%     hold on
%     
%     pd1 = fitdist(data,'Kernel','Kernel','epanechnikov','bandwidth',0.05);
%     p1=cdf(pd1,xx);
%     
%     plot(xx,p1,'r');
%     
%     P1(:,i)=p1';
%     
% end;
% 
% pd1 = fitdist(data1','Kernel','Kernel','epanechnikov','bandwidth',0.05);
% p1=cdf(pd1,xx);




% hh=[];
% for i=1:length(all)
%     hh=[hh;find(ic==all(i))];
%     end; 
%     
% maxx=max(num1(hh));  
% 
% for i=3:7
%     h=find(ic==i);
%     data=num1(h);
%     [f,x] = ecdf(data);
%     
%     figure(nn)
%     semilogx(x,f,'b')
%     hold on
% end;
% for i=9:12
%     i
%     h=find(ic==i);
%     data=num1(h);
%     [f,x] = ecdf(data);
%     
%     figure(nn)
%     semilogx(x,f,'r')
%     hold on
% end;
% 
% for i=13:17
%     i
%     h=find(ic==i);
%     data=num1(h);
%     [f,x] = ecdf(data);
%     
%     figure(nn)
%     semilogx(x,f,'k')
%     hold on
% end;
%     
% for i=[2,8,11,17,18]
%     i
%     h=find(ic==i);
%     data=num1(h);
%     [f,x] = ecdf(data);
%     
%     figure(nn)
%     semilogx(x,f,'g')
%     hold on
% end;

% means=zeros(1,length(All_rel));
% vars=means;
% 
% for i=1:length(All_rel)
%     i
%     h=find(ic==All_rel(i));
%     data=num1(h);
%     means(i)=mean(data);
%     vars(i)=[std(data)]^2;
%     
% end;
% 
% figure(77)
% scatter(means,vars,kk)
% hold on
% text(means+10, vars+100, C(All_rel));
% 
% names=C(All_rel);
% 
%     'LPS_1h', 'LPS_2h', 'LPS_4h','LPS_4h_Replicate', 'LPS_6h', 'Unstimulated'
%     'Unstimulated_Replicate'
% 
%     
%   
% figure(88)
% 
% 
% for i=1:5;
%     h=find(ic==LPS(i));
%     data=num1(h);
%     [f,x] = ecdf(data);
%     
%     subplot(1,6,i) 
%     semilogx(x,f,'k')
%     hold on
%     xlim([0,range])
%     ylim([0,1])
% end;
%     
% for i=1:5;
%     h=find(ic==PAM(i));
%     data=num1(h);
%     [f,x] = ecdf(data);
%     
%     subplot(1,6,i) 
%     semilogx(x,f,'b')
%     hold on
%     xlim([0,range])
%     ylim([0,1])
% end;
%     
% for i=1:5;
%     h=find(ic==PIC(i));
%     data=num1(h);
%     [f,x] = ecdf(data);
%     subplot(1,6,i) 
%     semilogx(x,f,'r')
%     hold on
%     xlim([0,range])
%     ylim([0,1])
% end;
% 
% 
% for i=1:5;
%     h=find(ic==Pert(i));
%     data=num1(h);
%     [f,x] = ecdf(data);
%     
%     subplot(1,6,6) 
%     if i==1
%     semilogx(x,f,'k')
%     hold on
% 
%     else
%     semilogx(x,f,'g')
%     hold on
%     xlim([0,range])
%     ylim([0,1])
%     end
% end;
%     
% figure(99)
% 
% 
% for i=1:5;
%     h=find(ic==LPS(i));
%     data=num1(h);
%     [f,x] = ksdensity(data);
%     
%     subplot(1,5,i) 
%     plot(x,f,'k')
%     hold on
%     xlim([0,range])
% end;
%     
% for i=1:5;
%     h=find(ic==PAM(i));
%     data=num1(h);
%     [f,x] = ksdensity(data);
%     
%     subplot(1,5,i) 
%     plot(x,f,'b')
%     hold on
%     xlim([0,range])
% end;
%     
% for i=1:5;
%     h=find(ic==PIC(i));
%     data=num1(h);
%     [f,x] = ksdensity(data);
%     
%     subplot(1,5,i) 
%     plot(x,f,'r')
%     hold on
%     xlim([0,range])
% end;
%     
%     

% gene=txt1(2:end,2);
% groups=txt1(2:end,13);
% 
% % index=zeros(1,length(gene));ind=[];
% % 
% % for i= 1:length(groups)
% % %     ind=find(gene==groups(i));
% %     ind=find(strcmp(gene, groups(i)));
% %     ind=[ind,ind]
% %     
% % end;
% %     
% 
% tf=ismember(gene,groups);
% %idx=[1:length(gene)];
% %idx=idx(tf);
% %idx=idx(loc(tf));
% %disp(gene(idx))
% 
% index=find(tf==1);
% 
% final=raw1(index+1,1:12);
