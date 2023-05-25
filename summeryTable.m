function Res=summeryTable(ta,grp,options)

arguments

    ta
    grp =[]
    options.Compact = false
end

%
% Makes a summery table of the variables in the table ta. If grp is assigned the
% summery statistics are done in each grp
%
% Res=summeryTable(ta,grp)
%  ta:  Table including the variable for analyse. Suported datatypes are: categorical, Logical and numeric
%  grp: Grouping variable (optional)
%
%  Res: Results table
%
% * For categorical variables numbers and propotions are counted
% * Logical variable numbers and propotions are counted
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



Res=table;



% in case defined groups prepared the grp variale ot categorical
if nargin>1


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


warning('off','MATLAB:table:RowsAddedExistingVars');
if options.Compact
    k=1;
    Res(k,1)={'n'};
    n=size(ta,1);
    Res(k,2)= {num2str(n)};

    for j=1:nugrp

        n(j)=sum(grp==ugrp(j));


        Res{k,j+2}={num2str(n(j))};
    end
else
    k=0;
end

% bulid the table
for i=1:size(ta,2)


    % check if  pseudoLogical
    if isnumeric(ta{:,i}) & sum( ~ismember(ta{:,i},[ 0 1]) &  ~isnan(ta{:,i}))<=0
        if sum(isnan(ta{:,i}))==0

            ta.(ta.Properties.VariableNames{i})=logical( ta{:,i});

        else
            ta.(ta.Properties.VariableNames{i})=categorical( ta{:,i},[nan 0 1],{'Undefined','N' ,'Y'});

        end
    end


    % Make table headings
    if   strcmp(ta{1,i},'Heading')
        k=k+1;
        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
            Res(k,1)={ta.Properties.VariableNames(i)};
        else
            Res(k,1)={ta.Properties.VariableDescriptions(i)};
        end


        Res{end,2}={''};

        if nugrp>0
            for j=1:nugrp
                Res{end,j+2}={''};
            end

        end

    elseif isnumeric(ta{:,i}) % in case of numeric data calculate statistics

        k=k+1;

        % estimate statistics for all data
        u=nanmean(ta{:,i}) ;
        m=nanmedian(ta{:,i});
        mi=nanmin(ta{:,i});
        ma=nanmax(ta{:,i});
        sd=nanstd(ta{:,i});
        n=sum(isnan(ta{:,i})==0);
        q(1:2,1)=quantile(ta{:,i},[0.25 0.75])';

        % estimate statistics for groups data
        if nugrp>0
            for j=1:nugrp
                x=ta{grp==ugrp(j),i};
                if ~isempty(x)
                    u(1,j+1)=nanmean(x) ;
                    m(1,j+1)=nanmedian(x);
                    mi(1,j+1)=nanmin(x);
                    ma(1,j+1)=nanmax(x);
                    sd(1,j+1)=nanstd(x);
                    n(1,j+1)=sum(isnan(x)==0);
                    q(1:2,j+1)=quantile(x,[0.25 0.75])';
                else
                    u(1,j+1)=nan;
                    m(1,j+1)=nan;
                    mi(1,j+1)=nan;
                    ma(1,j+1)=nan;
                    sd(1,j+1)=nan;
                    q(1:2,j+1)=nan;
                    n(1,j+1)=sum(isnan(x)==0);

                end
            end
        end

        % add statics for numeric data to table
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
        if nugrp>0
            Res{end,2:nugrp+2}={''};
        end
        if ~options.Compact
            Res{end+1,1}={' N'};

            for j=1:nugrp+1
                Res{end,1+j}= {num2str(n(j))};
            end


        end


        % estiate p-value
        if nugrp>1
            [p]=anova1(ta{:,i},grp,'off');

            if ~options.Compact

                Res{end+1,1}={[' Mean' setstr(177)  'SD (p=' num2str(p,3) ')']};
            else
                Res{end,1}={[char(Res{end,1}) ' (p=' num2str(p,3) ')']};
            end

        else
            if ~options.Compact


                Res{end+1,1}={[ ' Mean' setstr(177)  'SD'  ]};
            else

            end
        end


        for j=1:nugrp+1
            Res{end,1+j}= {[ num2str(u(j),'%.1f') setstr(177) num2str(sd(j),3)]};
        end
        if ~options.Compact
            Res{end+1,1}={' Median(IQR)'};
            for j=1:nugrp+1
                Res{end,1+j}= {[num2str(m(j),3) '(' num2str(q(1,j),3) '-' num2str(q(2,j),3) ')']};
            end
            Res{end+1,1}={' Min-Max'};
            for j=1:nugrp+1
                cif=floor( log10(ma(j)));

                if cif>=2
                    deci=cif+2;
                else
                    deci=3;
                end
                Res{end,1+j}= {['(' num2str(mi(j),deci) '-' num2str(ma(j),deci) ')']};
            end
        end

    elseif iscategorical(ta{:,i}) % In case of categorical datatype
        k=k+1;

        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
            Res(k,1)={cellstr([ta.Properties.VariableNames{i} '(N,(%))'])};
        else
            Res(k,1)={cellstr([ta.Properties.VariableDescriptions{i}  '(N,(%))'])};
        end

        if nugrp>1 % estimate p-value for categorical group data
            x=removecats(ta{:,i});
            [c ,ch,p]=crosstab(x,grp);
            temp_str=Res{k,1}{:};
            Res(k,1)={cellstr([temp_str ' (p=' num2str(p) ')'] )};

        end

        Res{end,2}={''};

        C = categories(ta{:,i});
        N=countcats(ta{:,i});
        if   ~options.Compact
            Res{end,2}={[num2str(sum(N))]};
        else

            Res{end,2}={''};
        end

        if nugrp>0 % estimate n for categorical group data
            for j=1:nugrp
                n=countcats(ta{ grp==ugrp(j),i});
                N(:,j+1)=n;
                if   ~options.Compact
                    Res{end,j+2}={[num2str(sum(N(:,j+1)))]};
                else
                    Res{end,j+2}={''};
                end
            end
        end

        % Set n in Res
        for k=1:length(C)

            Res{end+1,1}={[' ' C{k}]};

            for j=1:nugrp+1
                if j>1
                    Res{end,1+j}= {[num2str(N(k,j)) ' (' num2str(100*N(k,j)/sum(grp==ugrp(j-1) & ~isundefined( ta{:,i})),3) '%)']};
                    Res{end,1+j}=strrep( Res{end,1+j},'(NaN%)','');
                else
                    Res{end,1+j}= {[num2str(N(k,j)) ' (' num2str(100*N(k,j)/sum(N(:,j)),3) '%)']};

                end
            end

        end


    elseif islogical(ta{:,i})  % In case of logical datatype
        k=k+1;
        % write variable name in Res
        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})

            Res(k,1)={cellstr([ta.Properties.VariableNames{i} '(N,(%))'])};
        else
            str=cellstr([ta.Properties.VariableDescriptions{i} ' (N,(%))']);
            Res(k,1)={str};
        end

        n=sum(ta{:,i}==1);
        Res{end,2}={[num2str(n) ' (' num2str(100*n/length(ta{:,i}),3) '%)']};


        if nugrp>1

            for j=1:nugrp
                n=sum(ta{ grp==ugrp(j),i}==1);

                Res{end,j+2}={[num2str(n) ' (' num2str(100*n/sum(grp==ugrp(j)),3) '%)']};
            end


        end

    else
        disp(['Warning datatype not recognized for variable: -' ta.Properties.VariableNames{i} '- Consider to convert to categorical ']  )
    end
    k=size(Res,1);
end
warning('on','MATLAB:table:RowsAddedExistingVars')



% set table col names
grpNames=ugrp;
if nugrp>0
    for i=1:length(ugrp)
        %  grpNames{i}=strrep(grpNames{i},'-','_');
    end

    colNames={'Vars','All',grpNames{:}};
else
    colNames={'Vars','All'}  ;

end

Res.Properties.VariableNames=colNames;



