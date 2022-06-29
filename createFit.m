function [fitresult, gof] = createFit(z)
%CREATEFIT(Z)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      Y Output: z
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 13-Aug-2019 15:22:11


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( [], z );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.974942583142458 -0.00442780083297767];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


