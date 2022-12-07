# pubSTAT: publication of Statitics
The pubSTAT are aimed publication of statistical results in MATLAB. When doing statistical analyze results are often copy pasted manually into word or other editors for publication. However this might easily introduce errors and is time consuming. Therefore, I have developed pubSTAT that is a tool for simple statical analysis where the results are stored in tables. These tables can ealsy be converted to HTML tables using the publish function in MATLAB.

The current version includes 3 types of tables:
 * summeryTable: For basic statistical summery 
 * aucTables: For classification performance 
 * corrTables: For correlation  analyze 
 
When a table is contracted, it can be included in script which can be published using MATLAB's publish function. Central to convert the tables into HTML tables are the function plottable.m, that using MATLAB’s publishing function are converting MATLAB tables into HTML tables.

## Get started 
The basic idea is to make scripts that does the statistical analysis and then published these using MATLAB publish.
So, if you have stored you data in a table called "tbl" and you want to get summery statics, make a script like:
```
result_table=summeryTable(tbl);
 plotTable(result_table)
```

Then publish the script MATLAB publishing function  

Try to start to publish the file exampleReport.m, using:
```
the publish('exampleReport.m','showCode',true)
```
The resulting HTML file is HTML\exampleReport.html


##Citation:
Pleas cite:  Samuel E Schmidt. pubSTAT: publication of Statitics, github.com/Samuelemil/pubSTAT ,2022 


![image](https://user-images.githubusercontent.com/14206853/206318233-2e121f3c-29f8-4735-a2b3-751fbab92dcb.png)







