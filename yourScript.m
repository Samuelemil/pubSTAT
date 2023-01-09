% make a table with age and Sex from 5 subjects
tbl=table([20 33 45 23 57]',categorical({'Female','Male','Female','Male','Male'})','VariableNames',{'Age','Sex'});

% do summary analysis   
result_table=summeryTable(tbl);

% print the table 
plotTable(result_table) 