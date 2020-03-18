clear
clc

[num1,txt1,raw1]=xlsread('summary_mean_var_all_cond_sorted.xls');

[num2,txt2,raw2]=xlsread('induced_genes_outliers_removed_all_zero_removed.xls');

genelist=txt1(2:324,1);

singlenames=txt2(1:end,2);

condlist=txt2(1:end,3);

Pert={'On_Chip_Stimulation_LPS_4h', 'Ifnar1_KO_LPS_4h','IFNB_2h','Stat1_KO_LPS_4h','Tnfr_KO_LPS_4h','gpc_t4','hero_t0','ifnr_t275','gpa_t4','gpb_t4','tr1_t4','tr2_t4','tr3_t2','tr4_t2','wps_t4'};

Stim={'Unstimulated','Unstimulated_Replicate','LPS_1h','LPS_2h','LPS_4h','LPS_6h','LPS_4h_Replicate','PAM_1h','PAM_2h','PAM_4h','PAM_6h','PIC_1h','PIC_2h','PIC_4h','PIC_6h'};


figure(1)

genes={'TNF','IL1'}


for i=1:168
   
    
   h=find(strcmp(genelist(i),singlenames)>0);
   cond=condlist(h);
   
   x=num2(h,5);
   y=num2(h,6);
   
   subplot(12,14,i)
   scatter(x,y,'k')
   title(genelist(i))
   hold on
   
   h1=h(find(ismember(condlist(h),Pert)>0));
   x1=num2(h1,5);
   y1=num2(h1,6);
   scatter(x1,y1,'r')
   
   [p,s] = polyfit(x,y,1);
   plot(x,p(2)+p(1)*x,'b')
   hold on
   
   R2= 1 - (s.normr/norm(y - mean(y)))^2;
   tx=['R^2 ',num2str(round(R2,2))];
   text(10,max(y)-max(y)/5, tx);
   
   xlim([0,max(x)])
   ylim([0,max(y)])
   box on
   
end;
% % 
% % figure(2)
% % for i=1:155
% %    
% %     
% %    h=find(strcmp(genelist(i+168),singlenames)>0);
% %    cond=condlist(h);
% %    
% %    x=num2(h,5);
% %    y=num2(h,6);
% %    
% %    subplot(12,14,i)
% %    scatter(x,y,'k')
% %    title(genelist(i+168))
% %    hold on
% %    
% %    h1=h(find(ismember(condlist(h),Pert)>0));
% %    x1=num2(h1,5);
% %    y1=num2(h1,6);
% %    scatter(x1,y1,'r')
% %    
% %    [p,s] = polyfit(x,y,1);
% %    plot(x,p(2)+p(1)*x,'b')
% %    hold on
% %    
% %    R2= 1 - (s.normr/norm(y - mean(y)))^2;
% %    tx=['R^2 ',num2str(round(R2,2))];
% %    text(10,max(y)-max(y)/10, tx);
% %    
% %    xlim([0,max(x)])
% %    ylim([0,max(y)])
% %    box on
% %    
% % end;


% N=204;
% genes=genelist(1:N);
% 
% Min=min(num1(1:N,4));
% Max=max(num1(1:N,5));
% MinMax=min(num1(1:N,5));
% int=num1(1:N,2);
% slope=num1(1:N,3);
%     
% xx=[[0:5:100],[120:20:500],[550:50:1000],[1000:250:40000]];
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

