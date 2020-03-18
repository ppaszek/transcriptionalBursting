clear
clc

[num1,txt1,raw1]=xlsread('summary_mean_var_all_cond.xls');

%[num1,txt1,raw1]=xlsread('my_summary_burst_var_all_conds_log_robust.xls');

[num2,txt2,raw2]=xlsread('induced_genes_outliers_removed_all_zero_removed.xls');



genelist=txt1(2:end,1);
% 
singlenames=txt2(1:end,2);
% 
condlist=txt2(1:end,3);

Pert={'On_Chip_Stimulation_LPS_4h', 'Ifnar1_KO_LPS_4h','IFNB_2h','Stat1_KO_LPS_4h','Tnfr_KO_LPS_4h','gpc_t4','hero_t0','ifnr_t275','gpa_t4','gpb_t4','tr1_t4','tr2_t4','tr3_t2','tr4_t2','wps_t4'};

Stim={'Unstimulated','Unstimulated_Replicate','LPS_1h','LPS_2h','LPS_4h','LPS_6h','LPS_4h_Replicate','PAM_1h','PAM_2h','PAM_4h','PAM_6h','PIC_1h','PIC_2h','PIC_4h','PIC_6h'};

%genelist={'TNF','NFKBIA','IL6','CXCL10','CCL3'};

%F = @(p,xx) p(1)*xx.^p(2);
F = @(p,xx) p(1)+p(2)*xx.^p(3);
%x0=[10,0.0001,10];
x0=[1,0.1,0.1];
%R22=num1(:,6);
%h1=find(R22>0.68);
%genes=genelist(h1);

% NN=length(h1);
% Min=min(num1(1:NN,4));
% Max=max(num1(1:NN,5));
% MinMax=min(num1(1:NN,5));
% int=num1(1:NN,2);
% power=num1(1:NN,3);
% 
% xx=[0:0.01:5];
% 
% 
% AA=NaN*zeros(length(xx),NN);
% 
% for i=1:NN
%     h11=min(find(xx>num1(h1(i),4)));
%     h21=min(find(xx>num1(h1(i),5)));
% 
%     for j=h11:h21
%         AA(j,i)=int(i)*exp(xx(j)*power(i));
%     end;
% end;



P=zeros(5,3);

% 
for i=1:length(genelist)
   
   figure(1) 
   h=find(strcmp(genelist(i),singlenames)>0);
   %cond=condlist(h);
   
   x=(num2(h,9));
   y=log10(num2(h,6)+1);
   
   %x0=[max(x),1];
   
   subplot(3,3,i)
   scatter(x,y,'k')
   title(genelist(i))
   hold on
   
   [X,hi]=sort(x);
   Y=y(hi);
   X=X';
   Y=Y';
   
   %x0=[max(X),1];
   [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,X,Y);
   plot(X,F(p,X),'b')
   
   P(i,:)=p;
   
   Residuals=Y-F(p,X);   
   I = abs( Residuals) < 1.5 * std( Residuals );
   XX=X(find(I>0));
   YY=Y(find(I>0));
   scatter(XX,YY,'r')
%    
   [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,XX,YY);
   plot(XX,F(p,XX),'m')
   R2=1-sum((YY-F(p,XX)).^2)/(length(YY)*var(YY)); 
   tx=['R^2 ',num2str(round(R2,2))];
   text(1,max(Y)-max(Y)/5, tx);
%    
   %xlim([1,5])
   ylim([0,max(Y)])
   box on
   
   col={'b','b','b','r','r'}
   figure(2)
   plot(log10(XX),10.^F(p,XX),col{i})
   hold on
   
end;

%    
% %    figure(11)
% %    subplot(12,14,i)
% %    
% %    scatter(x,10.^y,'k')
% %    hold on
% %    plot(X,10.^F(p,X),'b')
% %    scatter(XX,10.^YY,'r')
% %    plot(XX,10.^F(p,XX),'m')
% %       
% %    xlim([0,max(X)])
% %    ylim([0,max(10.^Y)])
% %    
% end;
% 
% 
% figure(2)
% for i=1:155
%    
%     
%    h=find(strcmp(genelist(i+168),singlenames)>0);
%    x=log10(num2(h,9));
%    y=num2(h,6);
%    
%    subplot(12,14,i)
%    scatter(x,y,'k')
%    hold on
%    
%    title(genelist(i+168))
%    [X,hi]=sort(x);
%    Y=y(hi);
%    X=X';
%    Y=Y';
%    [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,X,Y);
%    plot(X,F(p,X),'b')
%    
%    Residuals=Y-F(p,X);   
%    I = abs( Residuals) < 1.5 * std( Residuals );
%    XX=X(find(I>0));
%    YY=Y(find(I>0));
%    scatter(XX,YY,'r')
% %    
%    [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,XX,YY);
%    plot(XX,F(p,XX),'m')
%    %P(i+100,:)=p; 
%  
%    R2=1-sum((YY-F(p,XX)).^2)/(length(YY)*var(YY)); 
%    tx=['R^2 ',num2str(round(R2,2))];
%    text(1,max(Y)-max(Y)/5, tx);
% %    
%    xlim([1,5])
%    ylim([0,max(Y)])
%    box on
%    
% end;
% 
% 
% % 
% % Min=min(num1(1:146,4));
% % Max=max(num1(1:146,5));
% % MinMax=min(num1(1:146,5));
% % int=num1(1:146,2);
% % slope=num1(1:146,3);
% %     
% % xx=[[0:5:100],[120:20:500],[550:50:1000],[1000:250:13000]];
% % 
% % 
% % AA=NaN*zeros(length(xx),146);
% % 
% % for i=1:146
% %     h1=min(find(xx>num1(i,4)));
% %     h2=min(find(xx>num1(i,5)));
% % 
% %     for j=h1:h2
% %         AA(j,i)=int(i)+slope(i)*xx(j);
% %     end;
% % end;
% 
% %[C,ia,ic] = unique(txt1);
% 
% % nn=2;
% % kk='b'
% % 
% % range=Inf;
% % 
% % LPS=[23,3,4,5,7];
% % PAM=[23,9,10,11,12];
% % PIC=[23,13,14,15,16];
% % all=[23,3,4,5,7,9,10,11,12,13,14,15,16];
% % Pert=[5,1,2,8,17,18];
% % 
% % All_rel=[1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,23,24];
% % 
% % 
% % 
% % P1=zeros(143,5);
% % 
% % xx=[[1:10:100],[200:50:1000],[1250:250:30000]];
% % 
% % S1h=[23,3,9,13];
% % S2h=[23,4,10,14];
% % S4h=[23,5,11,15];
% % S6h=[23,7,12,16];
% 
% 
% % figure(10)
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
