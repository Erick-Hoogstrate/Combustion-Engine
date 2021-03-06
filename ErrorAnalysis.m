%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\20193197\OneDrive - TU Eindhoven\ErrorNL.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 22-Mar-2021 12:36:47

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A1:B10";

% Specify column names and types
opts.VariableNames = ["cycle", "maxpressure"];
opts.VariableTypes = ["string", "double"];

% Specify variable properties
opts = setvaropts(opts, "cycle", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "cycle", "EmptyFieldRule", "auto");

% Import the data
ErrorFLE15 = readtable("C:\Users\20193197\OneDrive - TU Eindhoven\ErrorFLE15.xlsx", opts, "UseExcel", false);

% Maxpressure = data.maxpressure;
% Cycle = data.cycle;
maxpressure=ErrorFLE15.maxpressure;
cycle=ErrorFLE15.cycle;

%% Clear temporary variables
clear opts

%% Error Analysis
rows = 1:10; %number of rows

% figure()
% plot(maxpressure,'o') title('Index plot of dissolution fraction'); xlabel('Index'); ylabel('Result'); grid on; 
% 
% figure()
% histogram(maxpressure); title('Histogram of dissolution fraction'), xlabel('Result'), ylabel('Frequency'); grid on; 
% 
% figure()
% ksdensity(maxpressure); title('Density trace of dissolution fraction'), xlabel('Result'), ylabel('Density'); grid on; 
% 
% figure()
% boxplot(maxpressure); title('Boxplot of dissolution fraction'), xlabel('Result'); grid on; 


x = maxpressure;
EDA_x = [length(x),min(x),max(x),range(x)]
EDA_y=[mean(x),median(x),length(x)] 
EDA_z=[std(x),var(x),iqr(x)] 


% b)
mean_value = mean(maxpressure)
standard_diviation = std(maxpressure) %standard_diviation
length_value = length(maxpressure) 


x = mean(maxpressure) %mean value
u = (1+standard_diviation)/(sqrt(10)) %So a 68% confidence interval leads to x +/- u, (sqrt(10) because of 10 total values and 1 because that is the z-value for a 68% confidence interval)

se = std(maxpressure)/sqrt(length(maxpressure)) % 68% confidence interval 
