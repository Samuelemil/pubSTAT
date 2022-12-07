# pubSTAT: publication of Statitics
The pubSTAT are aimed publication of statistical results in MATLAB. When doing statistical analayse results are offent copy pasted manually into word or other editors for publication. However this migth easly introduce errors and is time consuming. Therefore I have devloped pubSTAT that is a tool for simple statical analysis where the results are stored in tables. These tables can ealsy be converet to HTML tables using the publish function in MATLAB.

The current version include 3 types of tables:
 * summeryTable: For basic statistical summery 
 * aucTables: For classification performance 
 * corrTables: For  correlation  analyse 
 
When a table is contructed, it can be included in script whihc can be published using MATLAB's publish function.

## Get started 
To get started publish the file exampleReport.m using the publish('exampleReport.m','showCode',true). The resulting HTML file is HTML\exampleReport.html

##Citation:
Pleas cite:  Samuel E Schmidt. pubSTAT: publication of Statitics, github.com/Samuelemil/pubSTAT ,2022 









