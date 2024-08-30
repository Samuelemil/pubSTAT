function Res=pairedTable(ta,grp)
%
% Makes a paired statistical analyse the variables in the table ta
%
% Res=pairedTable(ta,grp)
%  ta:  Table including the variables for analyse. Suported datatypes are
%  only numerical
%
%  Res: Results table
%
%
% * For continues variables a mean,STD, Median,(IQR), Min & max are calualted
%
%
% Examples using Matlab dataset:
%   load patients
%   T = table(Gender,Age,Height,Weight,Smoker,Systolic,Diastolic);
%   T=convertvars(T,{'Gender'},'categorical') % makes 'Gender' categorical
%
%   To get summery statistics of the whole dataset
%
%   Res=summeryTable(T)
%
%   To get the summery statistics of the whole dataset acording to
%   gender
%
%   Res=summeryTable(T,T.Gender)
%
%
%   To get the summery statistics of only Age, Diastolic & Systolic acording to
%   gender
%
%   Res=summeryTable(T(:,{'Age','Systolic','Diastolic'}),T.Gender)
%
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
arguments
    ta
    grp =[]


end

k=0;









i=1;

k=k+1;

warning('off')
Res=addanalyse(ta);



if ~isempty(grp)


    if ~iscategorical(grp)
        grp=categorical(grp);

    else
        grp=removecats(grp);

    end
    if sum(isundefined(grp) )>0
        idx_ex=isundefined(grp) ;
        disp(['Warning! ' num2str(sum(idx_ex) )  ' records excluded from analyse due to missing grp data'])
        ta(idx_ex,:)=[];

        grp(idx_ex)=[];
    end
    ugrp=categories(grp);
    nugrp=length(ugrp);
else
    ugrp=0;
    nugrp=0;
end

if nugrp>0

for i=1:nugrp
Res(end+1,1)={ugrp(i)};
Res(end,2:end)={'- - - -'};
Res=[Res ;addanalyse(ta(grp==ugrp(i),:))];
end

end

warning('on')

%end





colNames={'Vars','All'}  ;

%end
Res.Properties.VariableNames(1)={' '};
Res.Properties.VariableNames(2:end)=ta.Properties.VariableNames ;
end

function Res=addanalyse(ta)

warning('off');
Res=table;

n_all=sum(isnan(ta{:,:})==0,1);

touse=~isnan(sum(ta{:,:},2));
% estimate statistics for all data
u=nanmean(ta{  touse,:},1);
sd=nanstd(ta{touse,:},[],1);
n=sum(isnan(ta{touse,:})==0,1);



if ~isempty(ta(  touse,:))
pValue=rep1wayanova(ta(  touse,:));
else

pValue='';
end
if pValue<0.001
    pValueStr= '<0.001';
else

    pValueStr= num2str(pValue(1),3);
end


Res(1,1)={'n'};
Res(2,1)={'n (In analyse)'};
Res(3,1)={['Mean: (p='    pValueStr ')']};
for i=1:size(ta,2)

     Res(1,i+1)={[num2str(n_all(i))    ]};
    Res(2,i+1)={[num2str(n(i))    ]};
    Res(3,i+1)={[ num2str(u(i),3)  setstr(177)  num2str(sd(i),3)  ]};
end

end
