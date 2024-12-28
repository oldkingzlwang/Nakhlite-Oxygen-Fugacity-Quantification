% Title: Nakhlite Parental Magma Composition Generator
%
% Description:
% This script generates parental magma compositions (PMC) of nakhlites 
% based on the compositional range compiled from literature. It ensures 
% that the generated values maintain the same distribution as the real values 
% using a Gaussian Mixture Model (GMM).
%
% Requirements:
% - MATLAB R2023b or later (previous versions may also be applicable)
% - Statistics and Machine Learning Toolbox
% - Input file: 'PMC_input.xlsx' (must be in the same directory as the script)
%   The input file contains literature compilations of real PMC of nakhlites.
%
% Outputs:
% - 'PMC_output.pdf': A graphical comparison between generated and original PMCs.
% - 'PMC_output.xlsx': A dataset containing the generated PMCs.
%
% Usage:
% - Place the script and 'PMC_input.xlsx' in the same directory.
% - Run the script to generate PMC data and output files.
%
% Author: Zilong Wang
% Date: December 28, 2024
% Tested with: MATLAB R2023b
%
% Notes:
% - Ensure that MATLAB's Statistics and Machine Learning Toolbox is installed 
%   for the Gaussian Mixture Model implementation.
% - Modify the script to adjust the number of generated compositions if needed.

clear;clc;

data=readmatrix('PMC_input.xlsx','Sheet','PMC');
data=data(:,3:13);

[m,n] = size(data);

% 1) Fit a Gaussian Mixture Model (GMM) to approximate the data distribution.
%    Choose the number of mixture components 'k' (a positive integer).
k = 3;  % For example, try 3 components. Adjust as needed.

% Fit the GMM:
% 'RegularizationValue' helps avoid singular covariance matrices.
gmModel = fitgmdist(data, k, 'RegularizationValue', 1e-5);

% Suppose we want 100 new random points (also 11-dimensional).
numNewSamples = 100;

% 2) Generate random points from this fitted GMM.

xValidGenerated = [];
numAttempts = 0; % Initialize variables

while size(xValidGenerated, 1) < numNewSamples
    num_samples = numNewSamples * (numAttempts+1);  % Adjust factor as needed
    samples = random(gmModel, num_samples);

    % Delete the data that do not look like PMC
    % Condition 1: Non-negative values in all dimensions
    nonNegativeSamples = all(samples >= 0, 2);

    % Condition 2: Sum within [99.5, 100.5]
    rowSums = sum(samples, 2, 'omitnan');
    sumWithinRange = (rowSums >= 99.5) & (rowSums <= 100.5);

    % Combine both conditions
    validSamples = nonNegativeSamples & sumWithinRange;
    xValidGenerated = samples(validSamples, :);
    numAttempts = numAttempts + 1;

    % Optionally, prevent infinite loop
    if numAttempts > 1000  % Adjust as needed
        warning('Could not generate enough valid samples after multiple attempts.');
        break;
    end
end

% Select the required number of samples
if size(xValidGenerated, 1) >= numNewSamples
    newData = xValidGenerated(1:numNewSamples, :);
else
    error('Not enough valid samples generated. Try increasing numNewSamples or numAttempts.');
end

% newData is now a 100-by-11 matrix drawn from the fitted GMM.

% 3) (Optional) Statistical check that newData ~ original data
%    For instance, you could do a quick dimension-by-dimension check 
%    using a 2-sample Kolmogorov-Smirnov test or a 2-sample t-test. 
%    Here is a simple example using kstest2 for each dimension:
alpha = 0.01;
for dim = 1:n
    [h,p] = kstest2(data(:,dim), newData(:,dim), 'Alpha', alpha);
    fprintf('Dimension %d: p = %.3f -> %s\n', ...
        dim, p, ...
        ternary(h==0, 'no significant difference','significant difference'));
end

% Qualitatively comparison: visual comparison using diagrams
combined_data = [newData; data];
group_labels = [ones(size(newData, 1), 1); 2 * ones(size(data, 1), 1)];
colors = 'rb';  markers = '.'; % Red dots are generated new data, and
                               % blue dots are original PMC data
figure;
[h, ax, bigax] = gplotmatrix(combined_data, [], group_labels, colors, markers, [], 'on', '', '');
axgd=["SiO_2","TiO_2","Al_2O_3","Cr_2O_3","FeO","MgO","MnO","CaO","Na_2O","K_2O","P_2O_5"];
numVars = size(combined_data, 2);
for i = 1:numVars
    xlabel(ax(numVars, i), sprintf(axgd(i)));
    ylabel(ax(i, 1), sprintf(axgd(i)));
end

% Save the figures
figureUnits = 'centimeters';
figureHandle = get(groot,'CurrentFigure');
figW = 1300; figH = 1300;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'Position',[50 50 figW figH]);
set(gcf, 'PaperPositionMode', 'auto');
figureHandle.Renderer='Painters';
exportgraphics(gcf,'PMC_Output.pdf', 'ContentType', 'vector');
% print(figureHandle,'PMC_Output.png','-dpng','-r900');
close all

% Save the results
results=array2table(round([newData,sum(newData,2)],3),'VariableNames',...
        {'SiO2','TiO2','Al2O3','Cr2O3','FeO','MgO','MnO','CaO',...
        'Na2O','K2O','P2O5','Total'});
% writetable(results,'PMC_Output.xlsx');

function out = ternary(condition, valTrue, valFalse)
    if condition
        out = valTrue;
    else
        out = valFalse;
    end
end