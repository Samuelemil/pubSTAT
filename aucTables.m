function [ Res, AUC,pairSTAT] = aucTables(ta,t,Threshold,nboot,grp, plotROC,scoresOnSamePlot,pairedAnalysis, options)
%[ Res , AUC] = aucTables(ta,t,Threshold,nboot,grp)
%
% Estimates the AUC NPV  of the table ta.
% ta:           Tables with the scores for analysis (NxM)
% t:            Patient type (Dichotomous vector, binary or categorical
%               (Nx1)), I cases of categorical the 2. ordinal value is used
%               as a positve class
% Threshold:    Threshold for classification calcualtions (1x1 if one score or the
%               same for all scores or 1xM vector if different thresholds for each score )
% nboot:        Number of bootstraping runs for AUC CI calcualtion.
%               In case of 0 the parametic method of Hanley, J.A. et al from 1982 is used (default 0)
% grp:          Groupping variable, categorical vector (Nx1))
% plotROC:      Plot roc curve (default=1)
% scoresOnSamePlot:  Plot AUC from different scores in the same plot (default=0)
% pairedAnalysis:    compare significans levels across scores (note no grp's can be defined)
% Name-Value Arguments:  
%
% positiveClass {'posclass'} defienction of the postive class in t 
% 
% Res:          Results table
% AUC:          Table with numeric AUC values
%
% Examples using MATLAB dataset:
%   load patients
%   T = table(Gender,Age,Height,Weight,Smoker,Systolic,Diastolic);
%   T=convertvars(T,{'Gender'},'categorical') % makes 'Gender' categorical
%
%   To test of if systolic or diastolic blood pressure can predict if you are
%   smoking
%
%   Res = aucTables(T(:,{'Systolic','Diastolic'}),T.Smoker)
%
%   To test of if systolice blood pressure above 120 mmHg og Diastolice blood pressure above 90 mmHg can predict if you are
%   smoking
%
%   Res = aucTables(T(:,{'Systolic','Diastolic'}),T.Smoker,[120 90])
%
%   To test of if systolice blood pressure and Diastolice blood pressure predict
%   smoking in both males & females
%
%   Res = aucTables(T(:,{'Systolic','Diastolic'}),T.Smoker,[120 90],[],T.Gender)
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
    ta {mustBeTableOrnumeric(ta)}
    t {mustBeLogicalOrCategorical(t),mustHaveSameNrOfRows(ta,t)}
    Threshold =[]
    nboot=0
    grp {mustBeLogicalOrCategoricalEmpty(grp),mustHaveSameNrOfRows(ta,grp)} =[]
    plotROC=1
    scoresOnSamePlot=0
    pairedAnalysis=0
    options.positiveClass=[]
end

pairSTAT=[];
% if input ta is not a table, but a vector make it a table
if ~istable(ta)
    varName = inputname(1);
    ta_temp=table;
    ta_temp.(varName)=ta;
    ta=ta_temp;
end


if iscategorical(t)==0
    t=categorical(t);
end


ut=categories(t);

if length(ut)>2
    disp('Error : patient type has to be Dichotomous!')
    Res=table;
    return
end


% postive outcome is not number 2 revers the order of categories
if ~isempty(options.positiveClass)

     if ut(2)~=categorical({options.positiveClass})
        ut=ut([2 1]);
     end
end

has_t=~isundefined(t);
% prepare variables
hasRocPlot=1;
k=0;
Res=table;

if plotROC
    clf
end

if  ~isempty(grp)

    if islogical(grp)
        grp=categorical(grp) ;
    end
    grp=removecats(grp);
    ugrp=categories(grp);

    nugrp=length(ugrp);
else
    ugrp=0;
    nugrp=0;
end


if ~isempty(    Threshold )
    if length(Threshold)==1 & size(ta,2)>1
        % replicate thresholds in onlye one threshold is given
        Threshold=Threshold.*ones(1,size(ta,2));
    end
end


warning('off','MATLAB:table:RowsAddedExistingVars');
for i=1:size(ta,2)

    x=ta{:,i};
    k=k+1;

    if ~isempty(    Threshold )
        % statistic
        [npv(1,1), npv(2:3,1)] = binofit(sum(x(t==ut(1) )<=Threshold(i)),sum(x(t==ut(1) | t==ut(2) )<=Threshold(i)));
        [ppv(1,1), ppv(2:3,1)] = binofit(sum(x(t==ut(2) )>Threshold(i)),sum(x(t==ut(1) | t==ut(2) )>Threshold(i)));
        [sens(1,1), sens(2:3,1)] = binofit(sum(x(t==ut(2)) >Threshold(i)),sum(t==ut(2) & ~isnan(x)));
        [specs(1,1), specs(2:3,1)] = binofit(sum(x(t==ut(1)) <=Threshold(i)),sum(t==ut(1) & ~isnan(x)));
        [ruleOut(1,1), ruleOut(2:3,1)] = binofit(sum(x(has_t) <=Threshold(i)),sum(has_t & ~isnan(x)));
        [ruleIn(1,1), ruleIn(2:3,1)] = binofit(sum(x(has_t) >Threshold(i)),sum(has_t & ~isnan(x)));
        t2x2(1,1)= sum(x(t==ut(1) )<=Threshold(i));
        t2x2(2,1)= sum(x(t==ut(2) )<=Threshold(i));
        t2x2(1,2)= sum(x(t==ut(1) )>Threshold(i));
        t2x2(2,2)= sum(x(t==ut(2) )>Threshold(i));

        % Likelihood ratio
        LRP(1,1)=sens(1,1)/(1-specs(1,1));
        % Likelihood ratio positive
        LRN(1,1)=(1-sens(1,1))/specs(1,1);
        if pairedAnalysis
            topair(:,i)=x>Threshold(i);
        end
    end

    % estimate AUC
    [x_a y_a T au]=perfcurve(t,x,ut(2));
    auc=au(1);



    % plot ROC
    if plotROC

        if scoresOnSamePlot

            if nugrp>2
                subplot(ceil((nugrp+1)/3),3,1)
            elseif   nugrp>1
                subplot(1,  nugrp+1,1)
            else
                subplot(1,1,1)
            end

            if i>1
                hold on

            end

            plot(x_a(:,1),y_a(:,1))
            pbaspect([1 1 1])
        else
            if size(ta,2)>3
                subplot(ceil(size(ta,2)/3),3,i)
            else
                subplot(1,size(ta,2),i)
            end
            plot(x_a(:,1),y_a(:,1),'k')
            pbaspect([1 1 1])
        end

        hold off
        xlabel('1-Specificity')
        ylabel('Sensitivity')

        if scoresOnSamePlot

            if nugrp>1
                title( 'ROC curves: All' )
            else

                title( 'ROC curves' )

            end
        else

            if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
                title( ['ROC curve: ' ta.Properties.VariableNames(i)]);
            else

                title(['ROC curve: ' ta.Properties.VariableDescriptions(i)]);
            end

        end



    end


    n0(i,1)=sum(isnan(x)==0 & t==ut(1));

    n1(i,1)=sum(isnan(x)==0 & t==ut(2));

    if nboot>0


        auc_l=au(2);
        auc_h=au(3);
    else
        ac=aucConf(auc(1,1),n1(i,1),n0(i,1),2);
        auc_l=ac(1);
        auc_h=ac(2);
    end
    prv(i,1)=n1(i,1)/(n0(i,1)+n1(i,1));

    %% if more that one groups estimate all stat. for each group
    if nugrp>1

        for j=1:nugrp
            hasRocPlot(j)=0;

            xgrp=x(grp==ugrp(j));
            tgrp=t(grp==ugrp(j));
            has_t_grp=has_t(grp==ugrp(j));
            if ~isempty(xgrp) & sum(tgrp==ut(1))>=1 & sum(tgrp==ut(2))>=1


                n0(i,j+1)=sum(isnan(xgrp)==0 & tgrp==ut(1));
                n1(i,j+1)=sum(isnan(xgrp)==0 & tgrp==ut(2));

                [x_a y_a T au]=perfcurve(tgrp,xgrp,ut(2));

                auc(1,j+1)=au(1);

                if plotROC  & ~scoresOnSamePlot
                    hasRocPlot(j)=1;
                    hold on
                    plot(x_a(:,1),y_a(:,1))

                    hold off

                elseif plotROC  & scoresOnSamePlot
                    hasRocPlot(j)=1;
                    if nugrp>2
                        subplot(ceil( (nugrp+1)/3),3,j+1)

                    elseif   nugrp>1
                        subplot(1,  nugrp+1,j+1)
                    end
                    if i>1
                        hold on
                    end


                    plot(x_a(:,1),y_a(:,1))
                    if size(ta,2)>1
                        title(['ROC curves: '  ugrp{j}])
                    else
                        title(['ROC curve: '  ugrp{j}])
                    end
                    if i>1
                        hold on
                    end


                    if i==size(ta,2)
                        hold on
                        plot([ 0 1],[0 1],'k:')
                        hold off
                    end

                    hold off


                end

                prv(i,j+1)=binofit(n1(i,j+1),(n0(i,j+1)+n1(i,j+1)));

                ac=aucConf(auc(1,j+1), n1(i,j+1), n0(i,j+1),2);
                auc_l(1,j+1)=ac(1);
                auc_h(1,j+1)=ac(2);


                if ~isempty(    Threshold )
                    n=sum(xgrp(tgrp==ut(1) | tgrp==ut(2)) <=Threshold(i));

                    [npv(1,j+1), npv(2:3,j+1)] = binofit(sum(xgrp(tgrp==ut(1)) <=Threshold(i)),n)  ;
                    [ppv(1,j+1), ppv(2:3,j+1)] = binofit(sum(xgrp(tgrp==ut(2) )>Threshold(i)),sum(xgrp(tgrp==ut(1) | tgrp==ut(2) )>Threshold(i)));
                    [sens(1,j+1), sens(2:3,j+1)] = binofit(sum(xgrp(tgrp==ut(2)) >Threshold(i)),sum(tgrp==ut(2) & ~isnan(xgrp)));
                    [specs(1,j+1), specs(2:3,j+1)] = binofit(sum(xgrp(tgrp==ut(1)) <=Threshold(i)),sum(tgrp==ut(1) & ~isnan(xgrp)));

                    [ruleOut(1,j+1), ruleOut(2:3,j+1)] = binofit(sum(xgrp(has_t_grp) <=Threshold(i)),sum(has_t_grp & ~isnan(xgrp)));
                    [ruleIn(1,j+1), ruleIn(2:3,j+1)] = binofit(sum(xgrp(has_t_grp) >Threshold(i)),sum(has_t_grp & ~isnan(xgrp)));



                    t2x2(1,1,j+1)= sum(xgrp(tgrp==ut(1) )<=Threshold(i));
                    t2x2(2,1,j+1)= sum(xgrp(tgrp==ut(2) )<=Threshold(i));
                    t2x2(1,2,j+1)= sum(xgrp(tgrp==ut(1) )>Threshold(i));
                    t2x2(2,2,j+1)= sum(xgrp(tgrp==ut(2) )>Threshold(i));


                    % Likelihood ratio
                    LRP(1,j+1)=sens(1,j+1)/(1-specs(1,j+1));
                    % Likelihood ratio positive
                    LRN(1,j+1)=(1-sens(1,j+1))/specs(1,j+1);

                end


            else
                auc(1,j+1)=nan;
                n0(i,j+1)=sum(isnan(xgrp)==0 & tgrp==ut(1));
                n1(i,j+1)=sum(isnan(xgrp)==0 & tgrp==ut(2));
                auc_l(1,j+1)=nan;
                auc_h(1,j+1)=nan;
                npv(1,j+1)=nan;
                ppv(i,j+1)=nan;
                prv(i,j+1)=nan;
                sens(i,j+1)=nan;
                specs(i,j+1)=nan;
                ruleOut(i,j+1)=nan;
                ruleIn(i,j+1)=nan;
                t2x2(1,1,j+1)=nan;
                LRP(1,1,j+1)=nan;
                LRN(1,1,j+1)=nan;
            end
        end
    end

    %% Here all data are written into the result table
    % If VariableDescriptions is not empty use VariableDescriptions as
    % variable name
    if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
        Res(k,1)={ta.Properties.VariableNames(i)};
    else
        Res(k,1)={ta.Properties.VariableDescriptions(i)};
    end

    Res{end,2}={''};
    if nugrp>1
        Res{end,2:nugrp+2}={''};
    end

    Res{end+1,1}={['  N: ' ut{1} ]};


    for j=1:nugrp+1

        Res{end,1+j}= {num2str(n0(i,j))};
    end


    Res{end+1,1}={['  N: ' ut{2} ]};
    for j=1:nugrp+1

        Res{end,1+j}= {num2str(n1(i,j))};
    end

    if nugrp>1

        [ct,chi2,p_prev,labels]=crosstab(t,grp);
        Res{end+1,1}={['  Prevalence of ' ut{2}  ' (p=' num2str(p_prev,4)  ')']};
    else
        Res{end+1,1}={['  Prevalence of ' ut{2}  '' ]};
    end



    for j=1:nugrp+1

        Res{end,1+j}= {[num2str(prv(i,j)*100,3) '%']};
    end

    if nugrp>1


        tmpIdx=0;
        for j1=2:nugrp+1
            for j2=j1+1:nugrp+1

                tmpIdx=tmpIdx+1;
                %  [auc(j1) n0(i,j1) n1(i,j1) auc(j2) n0(i,j2) n1(i,j2)]
                p(tmpIdx)=aucCompareIndependent(auc(j1) ,n0(i,j1),n1(i,j1),auc(j2),n0(i,j2),n1(i,j2));
            end
        end

        p_auc=min(p);
        Res{end+1,1}={['  AUC (p=' num2str(p_auc,3)  ')']};
    else

        Res{end+1,1}={['  AUC' ]};
    end



    for j=1:nugrp+1
        Res{end,1+j}= {[num2str(auc(j)*100,3)  '% (' num2str(auc_l(j)*100,3) '-' num2str(auc_h(j)*100,3) ')']};
    end

    if ~isempty(    Threshold )

        if nugrp>1
            pt= x>Threshold(i);
            I=(t==ut(2));
            [ct,chi2,p_sens,labels]=crosstab(grp(I),pt(I));
            I=(t==ut(1));
            [ct,chi2,p_spec,labels]=crosstab(grp(I),pt(I));


            I=pt==0;
            [ct,chi2,p_npv,labels]=crosstab(grp(I),t(I));
            I=pt==1;
            [ct,chi2,p_ppv,labels]=crosstab(grp(I),t(I));

            [~ ,~,p_rule]=crosstab(grp(has_t),pt(has_t));

        end


        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
            if nugrp>1
                Res{end+1,1}={['  Negative predictive value (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')(p=' num2str(p_npv,4)  ')'  ]};
            else

                Res{end+1,1}={['  Negative predictive value (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')'  ]};

            end
        else

            if nugrp>1
                Res{end+1,1}={['  Negative predictive value (' ta.Properties.VariableDescriptions{i} '<=' num2str(Threshold(i)) ') (p=' num2str(p_npv,4)  ')'  ]};
            else
                Res{end+1,1}={['  Negative predictive value (' ta.Properties.VariableDescriptions{i} '<=' num2str(Threshold(i)) ')'  ]};
            end
        end

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(npv(1,j)*100,3)  '% (' num2str(npv(2,j)*100,3) '-' num2str(npv(3,j)*100,3) '%)']};
        end

        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})

            if nugrp>1
                Res{end+1,1}={['  Positive predictive value (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ') (p=' num2str(p_ppv,4)  ')'  ]};


            else
                Res{end+1,1}={['  Positive predictive value (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ')'  ]};


            end
        else
            if nugrp>1
                Res{end+1,1}={['  Positive predictive value (' ta.Properties.VariableDescriptions{i} '>' num2str(Threshold(i)) ') (p=' num2str(p_ppv,4)  ')' ]};
            else
                Res{end+1,1}={['  Positive predictive value (' ta.Properties.VariableDescriptions{i} '>' num2str(Threshold(i)) ')' ]};
            end
        end

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(ppv(1,j)*100,3)  '% (' num2str(ppv(2,j)*100,3) '-' num2str(ppv(3,j)*100,3) '%)']};
        end



        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
            if nugrp>1
                Res{end+1,1}={['  Sensitivity (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ') (p=' num2str(p_sens,4)  ')'   ]};
            else
                Res{end+1,1}={['  Sensitivity (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ')'   ]};
            end
        else
            if nugrp>1
                Res{end+1,1}={['  Sensitivity (' ta.Properties.VariableDescriptions{i} '>' num2str(Threshold(i)) ') (p=' num2str(p_sens,4)  ')'   ]};
            else

                Res{end+1,1}={['  Sensitivity (' ta.Properties.VariableDescriptions{i} '>' num2str(Threshold(i)) ')'   ]};

            end
        end

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(sens(1,j)*100,3)  '% (' num2str(sens(2,j)*100,3)  '-' num2str(sens(3,j)*100,3) '%)']};
        end

        % specificity
        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
            if nugrp>1
                Res{end+1,1}={['  Specificity (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')  (p=' num2str(p_spec,4)  ')'   ]};
            else
                Res{end+1,1}={['  Specificity (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')'   ]};
            end
        else
            if nugrp>1
                Res{end+1,1}={['  Specificity (' ta.Properties.VariableDescriptions{i} '<=' num2str(Threshold(i)) ')  (p=' num2str(p_spec,4)  ')' ]};
            else

                Res{end+1,1}={['  Specificity (' ta.Properties.VariableDescriptions{i} '<=' num2str(Threshold(i)) ')' ]};
            end
        end

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(specs(1,j)*100,3)  '% (' num2str(specs(2,j)*100,3)  '-' num2str(specs(3,j)*100,3) '%)']};
        end

        % rule-out
        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
            if nugrp>1
                Res{end+1,1}={['  Rule out (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')  (p=' num2str(p_rule,4)  ')'   ]};
            else
                Res{end+1,1}={['  Rule out (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')'   ]};
            end
        else
            if nugrp>1
                Res{end+1,1}={['  Rule out (' ta.Properties.VariableDescriptions{i} '<=' num2str(Threshold(i)) ')  (p=' num2str(p_rule,4)  ')' ]};
            else

                Res{end+1,1}={['  Rule out (' ta.Properties.VariableDescriptions{i} '<=' num2str(Threshold(i)) ')' ]};
            end
        end

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(ruleOut(1,j)*100,3)  '% (' num2str(ruleOut(2,j)*100,3)  '-' num2str(ruleOut(3,j)*100,3) '%)']};
        end

        % rule in
        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})
            if nugrp>1
                Res{end+1,1}={['  Rule in (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ')  (p=' num2str(p_rule,4)  ')'   ]};
            else
                Res{end+1,1}={['  Rule in (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ')'   ]};
            end
        else
            if nugrp>1
                Res{end+1,1}={['  Rule in (' ta.Properties.VariableDescriptions{i} '>' num2str(Threshold(i)) ')  (p=' num2str(p_rule,4)  ')' ]};
            else

                Res{end+1,1}={['  Rule in (' ta.Properties.VariableDescriptions{i} '>' num2str(Threshold(i)) ')' ]};
            end
        end

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(ruleIn(1,j)*100,3)  '% (' num2str(ruleIn(2,j)*100,3)  '-' num2str(ruleIn(3,j)*100,3) '%)']};
        end



        % TN
        Res{end+1,1}={['  TN (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')'     ]};

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(t2x2(1,1,j)) ]};
        end

        % FN
        Res{end+1,1}={['  FN (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')'     ]};

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(t2x2(2,1,j)) ]};
        end

        % FP
        Res{end+1,1}={['  FP (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ')'     ]};

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(t2x2(1,2,j)) ]};
        end

        % TP
        Res{end+1,1}={['  TP (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ')'     ]};

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(t2x2(2,2,j)) ]};
        end



        % Likelihood ratios
        Res{end+1,1}={['  Likelihood ratio positive (' ta.Properties.VariableNames{i} '>' num2str(Threshold(i)) ')'     ]};

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(LRP(1,j),4) ]};
        end


        Res{end+1,1}={['  Likelihood ratio negative (' ta.Properties.VariableNames{i} '<=' num2str(Threshold(i)) ')'     ]};

        for j=1:nugrp+1
            Res{end,1+j}= {[num2str(LRN(1,j),4) ]};
        end



 
    end

       if i<size(ta,2)
            Res{end+1,:}={' '};
       end



    if plotROC & scoresOnSamePlot==0
        hold on
        plot([ 0 1],[0 1],'k:')

        hold off

    elseif  plotROC & scoresOnSamePlot & i==size(ta,2)
        hold on
        plot([ 0 1],[0 1],'k:')

        hold off

    end

    k=size(Res,1);
    AUC(i,:)=auc;

end



grpNames=ugrp;
if nugrp>1
    for i=1:length(ugrp)
        grpNames{i}=strrep(grpNames{i},'-','_');
        grpNames{i}=[grpNames{i}];
    end

    colNames={'Vars','All',grpNames{:}};
else
    colNames={'Vars','All'}  ;

end



%% do paired statics analyse
if  nugrp<2 & size(ta,2)>1 &  pairedAnalysis
    t_bin=t==ut(2);
    notnan=~isnan(sum(ta{:,:},2));

    comb=nchoosek(1:size(ta,2),2);
    pairSTAT=table;
    for idx=1:size(comb,1)
        tmpcell= AUC_compare_correlated(ta{notnan & has_t,[comb(idx,:) ]}, t_bin(notnan & has_t));
        auc_p(idx)= tmpcell{end};
        if ~isempty(    Threshold )
            % McNemars Sensitivty
            ctab = crosstab(topair(notnan & t==ut(2),comb(idx,1)),topair(notnan & t==ut(2),comb(idx,2)));
            p_sen_pair(idx)=1-chi2cdf(((ctab(2,1)-ctab(1,2))^2)./(ctab(2,1)+ctab(1,2)),1);

            % McNemars Specificity
            ctab = crosstab(topair(notnan & t==ut(1),comb(idx,1)),topair(notnan & t==ut(1),comb(idx,2)));
            p_spec_pair(idx)=1-chi2cdf(((ctab(2,1)-ctab(1,2))^2)./(ctab(2,1)+ctab(1,2)),1);
            
             % chisq npv
           ns1=t_bin(topair(:,comb(idx,1))==0 & notnan  & has_t);
           ns2=t_bin(topair(:,comb(idx,2))==0 & notnan  & has_t);
           [c ,chi, p_npv_pair(idx)]=crosstab([ns1; ns2],[zeros(length(ns1),1) ; ones(length(ns2),1)]);


           % chisq ppv
           ns1=t_bin(topair(:,comb(idx,1))==1 & notnan  & has_t);
           ns2=t_bin(topair(:,comb(idx,2))==1 & notnan  & has_t);
           [c ,chi, p_ppv_pair(idx)]=crosstab([ns1; ns2],[zeros(length(ns1),1) ; ones(length(ns2),1)]);


             % chisq rule-out
           ns1=topair( notnan  & has_t,comb(idx,1));
           ns2=topair( notnan  & has_t,comb(idx,2));

            
            % McNemars Specificity
            ctab = crosstab(topair(notnan & has_t,comb(idx,1)),topair(notnan & has_t,comb(idx,2)));
             p_rule_pair(idx)=1-chi2cdf(((ctab(2,1)-ctab(1,2))^2)./(ctab(2,1)+ctab(1,2)),1);


        end
        if isempty(ta.Properties.VariableDescriptions)
        score1=ta.Properties.VariableNames(comb(idx,1));
        score2=ta.Properties.VariableNames(comb(idx,2));

        else
        score1=ta.Properties.VariableDescriptions(comb(idx,1));
        score2=ta.Properties.VariableDescriptions(comb(idx,2));
        end
        if isempty( score1{:})
            score1=ta.Properties.VariableNames(comb(idx,1));
        end


        if isempty( score2{:})
            score2=ta.Properties.VariableNames(comb(idx,2));
        end
        pairSTAT.('Paired analyse')(idx,1)={[score1{:} ' - ' score2{:}]};



    end


    pairSTAT.('auc p-value')(:,1)=auc_p;
    if ~isempty(    Threshold )
        pairSTAT.('Sensitivity  p-value')(:,1)=p_sen_pair;
        pairSTAT.('Specificity  p-value')(:,1)=p_spec_pair;
        pairSTAT.('NPV  p-value')(:,1)=p_npv_pair;
        pairSTAT.('PPV  p-value')(:,1)=p_ppv_pair;
        pairSTAT.('Rule IN/OUT p-value')(:,1)=p_rule_pair;
    end
end

%% plot ROC curve

if plotROC
    if scoresOnSamePlot
        for i=1:length(ta.Properties.VariableNames)

            if ~isempty(ta.Properties.VariableDescriptions) & ~isempty(ta.Properties.VariableDescriptions{i})
                leg(i)=ta.Properties.VariableDescriptions(i);
            else
                leg(i)=ta.Properties.VariableNames(i);
            end

        end

        legend(leg, 'Location','southeast')
    else
        if  nugrp>1
            legend(colNames([2 find(hasRocPlot)+2]), 'Location','southeast')
        end
    end
end
Res.Properties.VariableNames=colNames;


%% Rearrange the table to paired presentation
if  nugrp<2 & size(ta,2)>1 &  pairedAnalysis
    nrscore=size(ta,2);
    nrlines=(size(Res,1)-nrscore+1)/nrscore;
    r=1:nrlines;
    r=Res(1:nrlines,1);
    for i=1:nrscore


        tmp=Res([1:nrlines]+(i-1)*nrlines+i-1,2);
        if isempty(ta.Properties.VariableDescriptions) |  isempty(ta.Properties.VariableDescriptions{i})

            tmp.Properties.VariableNames(1)=ta.Properties.VariableNames(i);

        else


            tmp.Properties.VariableNames(1)=ta.Properties.VariableDescriptions(i);
        end
        r=[r tmp];
    end


    for i=1:size(r,1)

        idx_trm=  find(r{i,1}{:}== '(', 1, 'last');
        if ~isempty(idx_trm)

            r{i,1}{:}(idx_trm-1:end)=[];
        end

    end
       if ~isempty(    Threshold )
    r(end+1,1)={' Binary classification threshold'};
    for i=1:nrscore
        r(end,i+1)={num2str(Threshold(i))};
    end

       end
       r=[r(1:5,:); r(end,:); r(6:end-1,:)];
 Res=r;
end


%%  AUC store AUC data in  array 
Res.Properties.VariableNames(1)={'-'};
AUC=array2table(AUC);
AUC.Properties.RowNames= ta.Properties.VariableNames;
AUC.Properties.VariableNames= colNames(2:end);



% if plotROC
%     % Resize figure to double size
%     hFig = gcf;
%     pos=hFig.Position;
%     set(hFig, 'Position', [pos(1:2) pos(3:4)*2 ])
% end

end

% Custom validation function
function mustBeTableOrnumeric(x)
if ~istable(x) & ~isnumeric(x)
    ME = MException('aucTables:inputError','ta type must be table or numeric vector.');
    throwAsCaller(ME)
end
end

function mustBeLogicalOrCategorical(x)
if ~iscategorical(x) & ~islogical(x)
    ME = MException('aucTables:inputError','Type must be categorical or logical');
    throwAsCaller(ME)
end
end

function mustBeLogicalOrCategoricalEmpty(x)
if ~iscategorical(x) & ~islogical(x) & ~isempty(x)
    ME = MException('aucTables:inputError','Type must be categorical or logical');
    throwAsCaller(ME)
end
end

function mustHaveSameNrOfRows(ta,t)

if ~isempty(t)
    if (size(ta,1)~=size(t,1))
        ME = MException('aucTables:inputError','The number of rows must match the number of rows in position 1');
        throwAsCaller(ME)
    end
end

end