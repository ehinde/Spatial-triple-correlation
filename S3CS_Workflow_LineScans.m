%% Julissa Sanchez Velasquez, Elizabeth Hinde - March,2026
%
% S3CS_Workflow_LineScans.m - Code for performing cross-triple correlation analysis
% Last update: 03/2026
%______________________________________________________________________________________________________________________________________________
%
% OUTPUTS
%
% corrBins:           Correlation time in seconds.
% avg_BspecCellArray: A cell array where each cell contains the average S3CS profile for all tau1 and tau2 values from the 16 pixels analyzed.
%                     Each cell contains the average profile calculated for all tau1 and tau2 values.
% S3CS_avg:           A matrix containing the average S3CS values (along the tau1=tau2 line) from the 16 pixels analyzed.
% Fraction_S3CS:      A matrix containing the fraction of triple correlation 
% fitted_S3CS:        A matrix containing the fitted cross-S3CS curves using a Gaussian function 
% nPeakTable:         Table containing the translocation time (first column) and amplitude (second column)
% Adit_fvalues_pcf:   A matrix containing the additional fitted cross-S3CS curves using optimized starting values from the curveFitter MATLAB toolbox 
% nPeakTable2:        Table containing the translocation time (first column) and amplitude (second column) from the additional S3CS fitting

%% Step 1. Define the folder with all MATLAB functions

selected_folder = uigetdir('Select the folder containing the functions');
selected_path = genpath(selected_folder);
addpath(selected_path);


%% Step 2. Open TIFF files in a batch mode
%  Each TIFF file from the LSM880 microscope contains data for the three channels, the
%  code will assign each channel to Data1, Data2, and Data3 cell arrays, 
%  and then it will remove the first 10,000 lines from each line scan
     
[Data1, Data2, Data3, FileNames] = readTIFFFilesZeiss_remove10K_rows();

FileNames=FileNames';


%% Step 3. Fragment the data from each channel 
%  Change the value in 'fragmentSize' to increase or decrease the size of the fragments to be analyzed. 
%  A fragment size of 20,000 lines corresponding to ~7 s acquisition time is recommended. 

fragmentSize = 20000; 

[fragmentAverages_Data1, fragmentAverages_Data2, fragmentAverages_Data3] =...
    fragment_time_series2_cellArrays(Data1, Data2, Data3, fragmentSize);


%% Step 4. Calculate the ACF for each fragment to identify outliers 
%  Change the value in 'SampFreq' with the correct sampling frequency 

SampFreq = 616;
[corrBins, A_avg_data1, A_avg_data2, A_avg_data3] = Calculate_ACF(fragmentAverages_Data1, fragmentAverages_Data2,...
    fragmentAverages_Data3,SampFreq);


%% Step 5. Establish the threshold for the standard deviation to systematically discard fragments with artifacts 
%  Change the value in 'mumberSD' to indicate how many standard deviations (SD) from the mean SD we want to set as threshold. 
%  'mumberSD' = 3 is recommended as starting point.
%  Evaluate the 'Before' and 'After' ACF plots to find the correct threshold

% Threshold for Channel 1
mumberSD     = 3; 
delInd_Data1 = Identify_bad_ACF_curves(A_avg_data1, corrBins, mumberSD); 

%%
% Threshold for Channel 2
mumberSD     = 3; 
delInd_Data2 = Identify_bad_ACF_curves(A_avg_data2, corrBins, mumberSD);

%% 
% Threshold for Channel 3
mumberSD     = 3; 
delInd_Data3 = Identify_bad_ACF_curves(A_avg_data3, corrBins, mumberSD);


%% Step 6. Delete those fragments that have strong photobleaching or cell movement (thosed identified in Step 5)  
%  This function will identify the indices for those fragments with
%  artifacts and create a table with all the indices from the three channels.
%  This means that the fragments deleted from Ch_1 will be also deleted from Ch_2 and Ch_3 
%  to keep the temporal information.

uniqueDelInd = combine_deleted_indices(delInd_Data1, delInd_Data2, delInd_Data3);


%% Step 7. Order the channels data

[Filt_fragmentAverages_Data1, Filt_fragmentAverages_Data2, Filt_fragmentAverages_Data3] = filter_fragments_by_indices...
    (fragmentAverages_Data1, fragmentAverages_Data2, fragmentAverages_Data3, uniqueDelInd);


%% Step 8. Get the average matrix to be used for the triple correlation analysis
%  Answer 'yes' to round the average values or answer 'no' to consider
%  decimal points in the average data (answer 'no' or press enter for the latter).

[Fragmented_Data1_avg, Fragmented_Data2_avg, Fragmented_Data3_avg] = average_columns_in_cells_v2...
    (Filt_fragmentAverages_Data1, Filt_fragmentAverages_Data2, Filt_fragmentAverages_Data3);


%% OPTIONAL: Reverse the order of the average matrix for bidirectional analysis 
%  First column becomes last column, and last column becomes first column

Rev_Fragmented_Data1 = reverse_columns_in_cells(Fragmented_Data1_avg);
Rev_Fragmented_Data2 = reverse_columns_in_cells(Fragmented_Data2_avg);
Rev_Fragmented_Data3 = reverse_columns_in_cells(Fragmented_Data3_avg);


%% Step 9. CALCULATE TRIPLE CORRELATION 

Dist     = 0;    % 'Dist' = 0 for triple correlation analysis. 'Dist' > 0 to calculate cross-triple correlation
SampFreq = 3415; % 3415 for 300,000 lines using 7.81 us pixel dwell time. Change sampling frequency as needed.
m        = 220;  % BinSize for multiple tau appraoch. Recommended value. 
FirstCol = 1;    % First column to be analyzed
LastCol  = 16;   % Last column to be analyzed

[BspecCellArray, diagonalElementsAbs,corrBins, corrResult_array, S3CS_avg, avg_BspecCellArray,tau1,tau2] = plot3FCS_MulTau3_v3...
    (Fragmented_Data1_avg, Fragmented_Data2_avg, Fragmented_Data3_avg, FirstCol, LastCol, Dist, m, SampFreq);

% Activate the following lines if you want to perform bidirectional analysis (reverse orientation)

%[BspecCellArray, diagonalElementsAbs,corrBins, corrResult_array, S3CS_avg, avg_BspecCellArray,tau1,tau2] = plot3FCS_MulTau3_v2...
%    (Rev_Fragmented_Data1, Rev_Fragmented_Data2, Rev_Fragmented_Data3, FirstCol, LastCol, Dist, m, SampFreq);


%% Step 10. Plot carpets - Carpet representation 
%  Plot the S3CS output for each pixel along the line scan
%  Change 'smoothingVal' to change level of smoothing

smoothingVal = 3;
plotCorrResultCells(corrResult_array,smoothingVal)


%% Step 11. Plot the average 3D surfaces (for all tau1 and tau2 delay times) 
%  Plot the average S3CS profile from all pixels in the line scan

plotAvgBspecCellArray(avg_BspecCellArray, tau1, tau2)


%% OPTIONAL: Calculate and plot PCAs 
%  Calculate principal component analysis to distinguish data from independent proteins and 
%  heterotrimeric protein complexes (for visualization purposes only).
%  This function uses as input the output from S3CS analysis 

PCA_F3CS(BspecCellArray)

%% Step 12. Calculate the ACF for each channel to get the Fraction of triple correlation 
%  These functions will calculate the ACF for Ch1, Ch2, and Ch3 using the multiple tau approach. 
%  They will then identify the limiting channel and use the limiting channel's values to calculate 
%  the fraction of triple correlation. 

firstCol          = 1;         % First column to be analyzed
lastCol           = 16;        % Last column to be analyzed
sampleFreq        = 3415;      % Change the sampling frequency if necessary 
correction_factor = 0.0000215; % Adjust the correction factor according to calibration analysis in your microscope

Fraction_S3CS = ...
    FractionTriple2(firstCol,lastCol,sampleFreq,Fragmented_Data1_avg,Fragmented_Data2_avg,Fragmented_Data3_avg,S3CS_avg,correction_factor);

disp(Fraction_S3CS);



%% Step 13. Fit cross-S3CS data to a Gaussian 
%  READ: This subsection will fit the cross-S3CS to a Gaussian model using the fit function in MATLAB. 
%  Indicate if the fit will be performed using either one or two-term Gaussian models 
%  answering 'gauss1' or 'gauss2' in the command window. If no answer is provided, the code will fit the pCF to 
%  a one-term Gaussian model. We recommend using one Gaussian model, visualizing the fitting, 
%  and deciding if the data requires a two-term Gaussian model.

fitted_S3CS = GaussianPCF2(S3CS_avg,corrBins);


%% Step 14. Plot raw pCF and fitted data, get peak values
%  Get the translocation time (in seconds) from the fitted cross-S3CS data

peakTable       = plotGaussian1(fitted_S3CS,S3CS_avg,corrBins);
nPeakTable      = peakTable;
nPeakTable.x_pk = 10.^nPeakTable.x_pk;

disp(nPeakTable);

%% Step 15. Fit cross-S3CS data using optimized starting values 

% READ:

% The fit function in MATLAB uses default starting values and constraint bounds. 
% The default options usually produce an excellent fit; however, if the fit is not adequate, 
% you can specify suitable starting values, especially for the centroid (location), 
% using the MATLAB Curve Fitting Toolbox (a complete users guide can be found at 
% https://au.mathworks.com/help/pdf_doc/curvefit/index.html).
% The optimized starting values generated by this toolbox can then be inserted in Step 16. 
% To obtain those values, select the option "Export" from the Toolbox.
% By running Steps 16-18, the code will generate tables for the additional fitted values and parameters.

% We're showing an example about how to get individual S3CS profiles for the new fitting
% The indices in the example correspond to the S3CS profiles (from the S3CS_avg matrix) we want to fit again.

pCF_Cell1 = S3CS_avg(:,1);
pCF_Cell2 = S3CS_avg(:,2);
pCF_Cell3 = S3CS_avg(:,3);

Adit_pCF  = table(pCF_Cell1,pCF_Cell2,pCF_Cell3);
Adit_pCF  = table2array(Adit_pCF);


%% Step 16. Indicate StartPoints obtained from the curveFitter toolbox

% Values for pCF_Cell1, 2, 3, and so on

ft0 = fittype( 'gauss1' );
opts0 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts0.Display = 'Off';
opts0.Lower = [-Inf -Inf 0];
opts0.StartPoint = [1.76503316403817e-05 -3.06980451469502 0.40];
[f1] = fit(log10(corrBins),pCF_Cell1,ft0,opts0);

ft2 = fittype( 'gauss1' );
opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Display = 'Off';
opts2.Lower = [-Inf -Inf 0];
opts2.StartPoint = [1.76503316403817e-05 -3.06980451469502 0.40];
[f2] = fit(log10(corrBins),pCF_Cell2,ft2,opts2);

ft3 = fittype( 'gauss1' );
opts3 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts3.Display = 'Off';
opts3.Lower = [-Inf -Inf 0];
opts3.StartPoint = [1.76503316403817e-05 -3.06980451469502 0.40];
[f3] = fit(log10(corrBins),pCF_Cell3,ft3,opts3);


%% Step 17. Get additional f-values

numDatasets = 3; %Indicate the number of new data sets being evluated. Change accordingly

Adit_fvalues_pcf = zeros(length(corrBins), numDatasets);

for i = 1:numDatasets
    variableName = ['f', num2str(i)];  
    f = eval(variableName);  
    f = f(log10(corrBins));
   
    Adit_fvalues_pcf(:, i) = f;
end


%% Step 18. Plot the fitting and get new peak values

numDatasets = 3; %Indicate the number of new data sets being evluated. Change accordingly
peakTable2  = table();

for i = 1:numDatasets
    f = Adit_fvalues_pcf(:, i); 
    figure;  
    [y_pk, x_pk] = plotFittedCurves(corrBins, Adit_pCF(:, i), f);
    
    % Check if peaks were found
    if isempty(x_pk) || isempty(y_pk)
        x_pk = NaN;
        y_pk = NaN;
    end
    
    peakTable2 = [peakTable2; table(x_pk, y_pk)];
end

nPeakTable2      = peakTable2;
nPeakTable2.x_pk = 10.^nPeakTable2.x_pk;
disp(nPeakTable2);
