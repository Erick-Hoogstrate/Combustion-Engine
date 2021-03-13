clc; clear all; close all

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 41);

% Specify sheet and range
opts.Sheet = "Sheet2";
opts.DataRange = "A2:AO89";

% Specify column names and types
opts.VariableNames = ["Lot", "BlendSpeed", "Compressor", "Force", "SprayRate", "AtomizerPressure", "Dissolution", "obs1", "obs2", "obs3", "obs4", "obs5", "obs6", "obs7", "obs9", "obs10", "obs11", "obs12", "obs13", "obs14", "obs15", "obs16", "obs17", "obs18", "obs19", "obs20", "obs21", "obs22", "obs23", "obs24", "obs25", "obs26", "obs27", "obs28", "obs29", "obs30", "obs31", "obs32", "obs33", "obs34", "obs35"];
opts.SelectedVariableNames = ["Lot", "BlendSpeed", "Compressor", "Force", "SprayRate", "AtomizerPressure", "Dissolution", "obs1", "obs2", "obs3", "obs4", "obs5", "obs6", "obs7", "obs9", "obs10", "obs11", "obs12", "obs13", "obs14", "obs15", "obs16", "obs17", "obs18", "obs19", "obs20", "obs21", "obs22", "obs23", "obs24", "obs25", "obs26", "obs27", "obs28", "obs29", "obs30", "obs31", "obs32", "obs33", "obs34", "obs35"];
opts.VariableTypes = ["string", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 3], "EmptyFieldRule", "auto");

% Import the data
data = readtable("C:\Users\20193197\Documents\Group08.xlsx", opts, "UseExcel", false); %file location on computer

lot               = data.Lot;
blend_speed       = data.BlendSpeed;
compressor        = data.Compressor;
force             = data.Force;
spray_rate        = data.SprayRate;
atomizer_pressure = data.AtomizerPressure;
dissolution       = data.Dissolution;

% Clear temporary variables
clear opts

%% 
rows = 1:88 %number of rows

%% Task 2.d

% scatter(blend_speed,dissolution)                    %make scatter plots to see correlation
% scatter(force,dissolution)
%scatter(spray_rate,dissolution);
% scatter(atomizer_pressure,dissolution)


a=corrcoef(dissolution,blend_speed); %0.2116          %calculating correlation coefficient
b=corrcoef(dissolution,force); %0.1597
c=corrcoef(dissolution,spray_rate); %-0.3009
d=corrcoef(dissolution,atomizer_pressure); %0.1009

Simpel_reg = fitlm(spray_rate,dissolution)          %fitting a simple regression model
figure()
plot(Simpel_reg)

cof=coefCI(Simpel_reg)                              %confirmation that B0 and B1 are significant in 95% confidence interval

%% Task 2.e residual plot
%plotResiduals(Simpel_reg, 'fitted');
%plotResiduals(Simpel_reg, 'caseorder');

%plotResiduals(Simpel_reg)
%plotResiduals(Simpel_reg,'probability')

%% Task 2.f
disspray=table(spray_rate, dissolution);    
quad_reg=fitlm(disspray,'dissolution ~ spray_rate + spray_rate^2')  %make quadratic regression model
figure()
plot(quad_reg)

%Rsquared higher so better



