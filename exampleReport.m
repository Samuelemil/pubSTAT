%% Demonstration of pubSTAT 
% The current script document the use of pubSTAT using the MATLAB
% patients dataset. 
%  Author: Samuel E Schmidt sschmidt@hst.aau.dk

% publish the current document using:
%
%   publish('exampleReport.m','showCode',true)
% 
% output can be found in \html
%




load patients
T = table(Gender,Age,Height,Weight,Smoker,Systolic,Diastolic,Location,SelfAssessedHealthStatus);
% convert relevant variable to categorical
T=convertvars(T,{'Gender','SelfAssessedHealthStatus','Location'},'categorical');

% Define units 
T.Properties.VariableUnits('Height')={'Inch'};
T.Properties.VariableUnits('Weight')={'Pounds'};
T.Properties.VariableUnits('Age')={'Years'};
T.Properties.VariableUnits('Systolic')={'mmHg'};
T.Properties.VariableUnits('Diastolic')={'mmHg'};

%% Table 1: Data summary
t1=summaryTable(T);


plotTable(t1)


%% Table 2: Self Assessed Health Status vs blood pressure 
t2=summaryTable(T(:,{'Systolic','Diastolic'}),T.SelfAssessedHealthStatus);
plotTable(t2)


%% Table 3: Smoking vs blood pressure 
t3=summaryTable(T(:,{'Systolic','Diastolic'}),T.Smoker);

plotTable(t3)



%% Table 4: ROC curves & AUC. Can Systolic or Diastolic blood pressure predict if a patient smokes


t4 = aucTables(T(:,{'Systolic','Diastolic'}),T.Smoker);
plotTable(t4)

%% Table 5: ROC curves,AUC & classification. Can Systolic blood pressure >120 mmHg & Diatolic blood pressure >90 mmHg predict if a patient smokes

t5 = aucTables(T(:,{'Systolic','Diastolic'}),T.Smoker,[120 90]);
plotTable(t5)


%% Table 6: ROC curves & AUC & classification. Gender differences in AUC 

t6 = aucTables(T(:,{'Systolic','Diastolic'}),T.Smoker,[120 90],[],T.Gender);
plotTable(t6)



%% Table 7: Correlation between Age, Heigth or Weigth and Systolic blood pressure:
t7=corrTables(T(:,{'Age','Height','Weight'}),T.Systolic);
plotTable(t7)


%% Table 8: Correlation between Age, Heigth or Weigth and Systolic blood pressure in Genders:
t7=corrTables(T(:,{'Age','Height','Weight'}),T.Systolic,T.Gender);
plotTable(t7)
