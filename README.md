# pubSTAT: publication of Statitics
pubSTAT is aimed publication of statistical results in MATLAB. When doing statistical analyze results are often copy pasted manually into Word or other editors for publication. However this might easily introduce errors and is time consuming. Therefore, I have developed pubSTAT, that is a tool for simple statical analysis where the results are stored in tables. These tables can ealsy be converted to HTML tables using the function plotTable.m and the publish function in MATLAB.

The current version includes 3 types of statistical tables:
 * summeryTable: For basic statistical summery 
 * aucTables: For classification performance 
 * corrTables: For correlation  analyze 
 
When a table is constructed, it can be included in a script which can be published using MATLAB's publish function. Central to convert the tables into HTML tables are the function plottable.m, that when using MATLAB’s publishing function converts MATLAB tables into HTML tables.

## Get started 
The basic idea is to make scripts that does the statistical analysis and then published these using MATLAB publish.
### Example 1
So, if you store your data in a table called "tbl" and you want to get summery statics, make a script called 'yourScript.m' including the code below:
```
% make a table with age and Sex from 5 subjects
tbl=table([20 33 45 23 57 ]',categorical({'Female','Male','Female','Male','Male'})','VariableNames',{'Age','Sex'});

% do summary analysis  
result_table=summeryTable(tbl);
% print the table
 plotTable(result_table)
```

Then publish the script MATLAB publishing function  
```
publish('yourScript.m')
```
The resulting HTML file is HTML\yourScript.html

### Example 2
Publish the file 'exampleReport.m', using:
```
publish('exampleReport.m','showCode',true)
```
The resulting HTML file is HTML\exampleReport.html

'exampleReport.m' includes use of all the differente statistical tables using the MATLAB patient dataset.


## Citation

Pleas cite:  Samuel E Schmidt. pubSTAT: publication of Statitics, github.com/Samuelemil/pubSTAT ,2022 


[![View pubSTAT on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/121832-pubstat)

![image](https://user-images.githubusercontent.com/14206853/206318233-2e121f3c-29f8-4735-a2b3-751fbab92dcb.png)







