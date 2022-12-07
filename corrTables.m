function [ Res,R] = corrTables(ta,y,grp)
% Estimates Pearson correlation between variables in the table ta with the variabe y
%
% [ Res,R ] = corrTables(ta,y,grp)
%  ta:  Table including the variable for correlation analyse with y. Suported datatypes are:
%       numeric and ordinal categoricals
%  y:   Reference Variable for correlation analyse (Vector)
%  grp: Grouping variable
%
%  Res: Correlation results table
%  R:   Raw numeric correlation values. Pretty in Matlab , but not for HTML
%
%
% Examples using Matlab dataset:
%   load patients
%   T = table(Gender,Age,Height,Weight,Smoker,Systolic,Diastolic);
%
%   To test if Age, Heigth or Weigth correlates to Systolic blood pressure:
%
%   [Res]=corrTables(T(:,{'Age','Height','Weight'}),T.Systolic)
%
%
%   To test if the correlations are equal between across Genders:
%
%   T=convertvars(T,{'Gender'},'categorical') % makes 'Gender' categorical
%  [Res]=corrTables(T(:,{'Age','Height','Weight'}),T.Systolic,T.Gender)
%
%   To get the raw correlation values
%  [Res, R]=corrTables(T(:,{'Age','Height','Weight'}),T.Systolic,T.Gender)
%
% Tip:  To replace varibale names with more describing names use
%       VariableDescriptions in the table, like
%
%       T.Properties.VariableDescriptions('Age')={'Subject age'};
%
%       Use VariableUnits to include units in the results
%
%       T.Properties.VariableUnits('Age')={'Years'};
%
%  Author: Samuel E Schmidt sschmidt@hst.aau.dk


k=0;
Res=table;
if nargin>2
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

    [r,p]=corr(x,y,'rows','pairwise');
    rs=r^2;
    n=sum(isnan(x)==0 & isnan(y)==0);


    if nugrp>1
        for j=1:nugrp
            xgrp=x(grp==ugrp(j));
            ygrp=y(grp==ugrp(j));
            if ~isempty(xgrp)
                %[r,p]=corr(x,y,'rows','pairwise');
                [r(1,j+1),p(1,j+1)]=corr(xgrp,ygrp,'rows','pairwise');
                n(1,j+1)=sum(isnan(xgrp)==0 & isnan(ygrp)==0);
            else
                p(1,j+1)=nan;
                r(1,j+1)=nan;
                n(1,j+1)=sum(isnan(x)==0);
                n(1,j+1)=nan;
            end
        end
    end

    R(i,:)=r;



    if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
        Res(k,1)={ta.Properties.VariableNames(i)};
    else
        Res(k,1)={ta.Properties.VariableDescriptions(i)};
    end
    if isempty(ta.Properties.VariableUnits) |  isempty(ta.Properties.VariableUnits{i})

    else

        e=cellstr([ char(Res{k,1}) ' ('  ta.Properties.VariableUnits{i} ')']);
        Res(k,1)={e};
    end

    Res{end,2}={''};
    if nugrp>1
        Res{end,2:nugrp+2}={''};
    end
    Res{end+1,1}={' N'};


    for j=1:nugrp+1

        Res{end,1+j}= {num2str(n(j))};
    end

    Res{end+1,1}={' r'};

    for j=1:nugrp+1
        Res{end,1+j}= {num2str(r(j)^1,3)};
    end

    Res{end+1,1}={' p value'};

    for j=1:nugrp+1

        Res{end,1+j}= {num2str(p(j),2)};
    end




    k=size(Res,1);

end



grpNames=ugrp;
if nugrp>1
    for i=1:length(ugrp)


        grpNames{i}=strrep(grpNames{i},'-','_Minus_') ;

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
