clear
clc

[num1,txt1,raw1]=xlsread('summary_mean_noise_all_cond_sort.xls');

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
F = @(p,xx) p(1)*xx.^p(2);
%x0=[10,0.0001,10];
x0=[10,-0.1];


R22=num1(:,6);
h1=find(R22>0.65);
genes=genelist(h1);

NN=length(h1);
Min=min(num1(1:NN,4));
Max=max(num1(1:NN,5));
MinMax=min(num1(1:NN,5));
slope=num1(1:NN,2);
power=num1(1:NN,3);


xx=[[0:5:100],[120:20:500],[550:50:1000],[1000:250:30000]];


AA=NaN*zeros(length(xx),NN);

for i=1:NN
    h11=min(find(xx>num1(h1(i),4)));
    h21=min(find(xx>num1(h1(i),5)));

    for j=h11:h21
        aa=slope(i)*xx(j)^power(i);
        AA(j,i)=aa;
    end;
end;


% NNN=zeros(323,1);
% P=zeros(323,2);
% RR2=zeros(323,1);
% 
% for i=1:168
%    
%    figure(1) 
%    h=find(strcmp(genelist(i),singlenames)>0);
%    %cond=condlist(h);
%    
%    x=(num2(h,5));
%    y=(num2(h,11));
%    
%    %x0=[max(x),1];
%    
%    subplot(12,14,i)
%    scatter(x,y,'k')
%    title(genelist(i))
%    hold on
%    
%    [X,hi]=sort(x);
%    Y=y(hi);
%    X=X';
%    Y=Y';
%    
%    %x0=[max(X),1];
%    [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,X,Y);
%    plot(X,F(p,X),'b')
%    
%    
%    
%    Residuals=Y-F(p,X);   
%    I = abs( Residuals) < 1.5 * std( Residuals );
%    XX=X(find(I>0));
%    YY=Y(find(I>0));
%    scatter(XX,YY,'r')
% %    
%    [p,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,XX,YY);
%    plot(XX,F(p,XX),'m')
%    R2=1-sum((YY-F(p,XX)).^2)/(length(YY)*var(YY)); 
%    tx=['R^2 ',num2str(round(R2,2))];
%    text(max(X)/2,max(Y)-max(Y)/5, tx);
% %   
%    P(i,:)=p;
%    RR2(i)=R2;
%    NNN(i)=length(XX);
%    xlim([1,max(X)])
%    ylim([0,max(Y)])
%    box on
%    
% %    col={'b','b','b','r','r'}
% %    figure(2)
% %    plot(log10(XX),10.^F(p,XX),col{i})
% %    hold on
%    
% end;
% 
% 
% figure(2)
% 
% for i=1:155
%    
%     
%    h=find(strcmp(genelist(i+168),singlenames)>0);
%    x=(num2(h,5));
%    y=(num2(h,11));
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
%    text(max(X)/2,max(Y)-max(Y)/5, tx);
% %  
%    P(i+168,:)=p;
%    RR2(i+168)=R2;
%    NNN(i+168)=length(XX);
%    xlim([1,max(X)])
%    ylim([0,max(Y)])
%    box on
%    
% end;
