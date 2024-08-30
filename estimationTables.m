function  Res=estimationTables(ta,y,grp)



k=0;
Res=table;
if nargin>2
    %grp(grp<0)=(grp(grp<0)*-1)+10;
    ugrp=categories(grp);
    nugrp=length(ugrp);
else
    ugrp=0;
    nugrp=0;
end


warning('off','MATLAB:table:RowsAddedExistingVars');
for i=1:size(ta,2)
    
    x=ta{:,i};
    
    
    
    if iscategorical(x)
        
        if  isordinal(x)
            x=grp2idx(x);
            
        else
            disp(['ERROR: ' ta.Properties.VariableNames(i) ' is a non ordinal categorical '])
             continue
        end
    end
    
    k=k+1;
 
    [ma ,p,maci(:,1),b ,r , bci(:,1),rci(:,1) , b_sig_p(:,1)]=regressionMeasures(x,y);
   
    rs=ma^2;
    n=sum(isnan(x)==0 & isnan(y)==0);
    
    
    if nugrp>1
        for j=1:nugrp
            xgrp=x(grp==ugrp(j));
            ygrp=y(grp==ugrp(j));
            if ~isempty(xgrp)
                [ma(1,j+1)  p(1,j+1),maci(:,j+1), b(1,j+1),r(1,j+1) , bci(:,j+1),rci(:,j+1), b_sig_p(:,j+1)]=regressionMeasures(xgrp,ygrp);
           
                n(1,j+1)=sum(isnan(xgrp)==0 & isnan(ygrp)==0);
            else
                p(1,j+1)=nan;
                b(1,j+1)=nan;
                ma(1,j+1)=nan;
                n(1,j+1)=sum(isnan(x)==0);
                    n(1,j+1)=nan;


            end
        end
    end
    
    R(i,:)=ma;
    
    
    Res(k,1)={ta.Properties.VariableNames(i)};
    Res{end,2}={''};
    if nugrp>1
        Res{end,2:nugrp+2}={''};
    end
    Res{end+1,1}={' N'};
    
    
    for j=1:nugrp+1
        
        Res{end,1+j}= {num2str(n(j))};
    end
    
    
    Res{end+1,1}={' Bias'};
    
    for j=1:nugrp+1
        
        Res{end,1+j}= {[num2str(b(j),2) ' (' num2str(bci(1,j),2)  ',' num2str(bci(2,j),2) ')' ]};
             if  b_sig_p(:,j)<0.05

 Res{end,1+j}={[ char(Res{end,1+j}) '*']};
     end

    end


    Res{end+1,1}={' MAPE (%)'};
    
    for j=1:nugrp+1
        Res{end,1+j}= {[num2str(ma(j)^1,3)   ' (' num2str(maci(1,j),3)  ',' num2str(maci(2,j),3) '%)']};
    end
    
      Res{end+1,1}={' r'};
    
    for j=1:nugrp+1
        
        Res{end,1+j}= {[num2str(r(j),3)  ' (' num2str(rci(1,j),2)  ',' num2str(rci(2,j),2) ')' ] };
    end

    Res{end+1,1}={' SEE'};
    
    for j=1:nugrp+1
        
        Res{end,1+j}= {num2str(p(j),2)};
    end

    

    
    
    k=size(Res,1);
    
end



grpNames=ugrp;
if nugrp>1
    for i=1:length(ugrp)
        
        
        grpNames{i}=strrep(strrep(grpNames{i},'-','_Minus_'),' ','') ;
        
        grpNames{i}=[grpNames{i}];
    end
    colNames={'Vars','All',grpNames{:}};
    
  
else
    colNames={'Vars','All'}  ;
    
end

Res.Properties.VariableNames=colNames;
%varNames=ta.Properties.VariableNames;

R=array2table(R);
R.Properties.RowNames=ta.Properties.VariableNames;
R.Properties.VariableNames=colNames(2:end);
[~,idx]=sort(R.All);

R=R(idx,:);
if nugrp>1
 
   nan_idx= find(isnan(n));
 
    Res(:, nan_idx+1)=[];

end

end

function [m see m_ci bias r bias_ci rci b_sig_p]=regressionMeasures(predicted,ref)

APE=abs(predicted-ref)./ref;
m=nanmean(APE)*100;
N=sum(~isnan(APE));
                                   
SEM = nanstd(APE)/sqrt(N);                             
ts = tinv([0.025 0.975], N-1);   
m_ci=m+  100*(ts*SEM);
see=sqrt(nansum((predicted-ref).^2)./(sum(~isnan(predicted-ref))-2));
bias=nanmean(predicted-ref);
%bias_ci=[bias-1.96*nanstd(predicted-ref)/sqrt(N) bias+1.96*nanstd(predicted-ref)/sqrt(N)];
[r,~, rciLow, rcihigh]=corrcoef(predicted,ref,'rows','pairwise');
r=r(2,1);
rci=[rciLow(2,1) rcihigh(2,1) ];
[b_sig b_sig_p bias_ci]=ttest(predicted-ref);
R2=1-sum((predicted-ref).^2)/sum((ref-mean(ref)).^2)

end
