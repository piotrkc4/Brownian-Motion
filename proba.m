clc
clear all
close all
N = 1E4; %number of experiments, eg Monte Carlo
T = 3E3; %number of timestamps
wHat= randn(T,N); %at each point in time you get a distribution of probabilities. 
ksData = NaN(T,N);
binVals = linspace(0,1, N); %assign the probabilities to one of these buckets
    for t = 1: T
        ksData(t,:) = ksdensity(wHat(t,:), binVals);
    end
    surf(1:T, binVals, ksData');
    colormap jet;
    shading interp;
    colorbar;
    view(-37.50, 30);
    set(gca, 'XLIM', [0 T]);
    set(gcf, 'color', 'white');
    grid on;
    

    