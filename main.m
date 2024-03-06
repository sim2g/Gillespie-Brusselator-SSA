%--------------------------------------------------------------------------
% This is the main script which calls functions from the class
% ChemicalReaction, for all questions in the assignment. This algorithm was
% written in an object-oriented fashion for robustness.
%--------------------------------------------------------------------------
% clearing values and figures for each run to reduce runtime and memory allocation as
% well as clearing the command window
clear all; clc; close all;
% beginning the run timer to find total run time
tic
% the default font size in figures is set to 11 to match the latex document
set(0,'defaultAxesFontSize',11);
% background of figures is set to white for aesthetic purposes.
set(0,'defaultfigurecolor','w');
% 6 objects are created. R1(1) and R1(2) are for question 1(a) and 1(b),
% R2(1), R2(2) and R2(3) are for question 2 with sizes S=0.1,1,10
% respectively. R_prob_dist is for the longer simulation in question 2.
R1 = [ChemicalReaction(1,1) ChemicalReaction(1,2)];
R2 = [ChemicalReaction(2,1) ChemicalReaction(2,2) ChemicalReaction(2,3)];
R_prob_dist = ChemicalReaction(2,1);
% Calculating and plotting the 10 realisations for question 1(a) as well as
% the deterministic solutions and outputting the NRMSD values for X and Y.
for i = 1:10
    Ra = R1(1).calculate(1);
    [nrmsd_ay(i),nrmsd_ax(i)] = Ra.realisation1(i);
    Ra.deterministic(1);
end
hold off
% new figure is created for question 1(b)
figure
% Calculating and plotting the 10 realisations for question 1(b) as well as
% the deterministic solutions and outputting the NRMSD values for X and Y.
for j = 1:10
    Rb = R1(2).calculate(2);
    [nrmsd_by(j), nrmsd_bx(j)] = Rb.realisation1(j);
    Rb.deterministic(2);
end
% displaying the mean NRMSD values for X and Y in question 1(a) and 1(b)
ChemicalReaction.mean_nrmsd(nrmsd_ax,nrmsd_ay,nrmsd_bx,nrmsd_by);
% specifying a figure for the sub-plots of the PSD function to be plotted
% on as they are plotted in a loop
fig = figure();
% k = 1,2,3 corresponds to S=0.1,1,10 in question 2a
for k = 1:3
    % the realisations of X and Y are calculated and plotted in separate
    % figures for different S values
    figure
    Rc = R2(k).calculate();
    Rc.realisation2();
    % The PSD is calculated for the iteration's S value, the PSD figure is selected and the subplot for the S value specified
    % in the iteration number is plotted.
    figure(fig);
    Rc.psd(k);
end
% The t_max and delta t_sample values are changed for the longer
% simulation
R_prob_dist.t_max = 100;
R_prob_dist.t_delta= 0.001;
% The X and Y realisations are calculated
R_pd = R_prob_dist.calculate();
% A new figure is created for the probability distributions to be plotted on
figure
% The probability distributions are plotted and the chi-square test results
% are outputted.
R_pd.prob_dist();
% The estimated mean, estimated standard deviation and skewness are plotted
% and final values displayed.
figure
R_pd.quantities();
% The proportion of the data above and below the mean, within 1, 2 and 3
% standard deviations of the mean and the median are calculated.
R_pd.additional_quantities();
% the run timer is ended and total time elapsed is displayed.
toc