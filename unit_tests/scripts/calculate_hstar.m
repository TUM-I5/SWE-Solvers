% Author: Alexander Breuer, breuera AT in.tum.de
%
% This function computes the middle state of a given homogenous
% Riemann-problem

function[hStar] = calculate_hstar(hLow, hHigh, huLow, huHigh)
%%% constants
zeroTol = 0.00000001;
dryTol = 0.001;
maxWaterHeight = 10000^3;

g = 9.81;

%%% Function phi

phi = @(h) (h <= hLow) .* 2 .* ( sqrt(g.*h) - sqrt(g.*hLow) ) ...
    +      (h  > hLow) .* (h - hLow) .* sqrt( g/2 .* ( 1 ./ h + 1 ./ hLow ) ) ...
    ...
    +      (h <= hHigh) .* 2 .* ( sqrt(g.*h) - sqrt(g.*hHigh) ) ...
    +      (h  > hHigh) .* (h - hHigh) .* sqrt( g/2 .* ( 1 ./ h + 1 ./ hHigh ) ) ...
    ...
    + huHigh ./ hHigh - huLow ./ hLow;

%%%
hMin = min(hLow, hHigh);
%hMax = max(hLow, hHigh)
if  phi(hMin) >= 0 && ...
    huLow / hLow - huHigh / hHigh + 2*(sqrt(g.*hLow)+sqrt(g.*hHigh)) < zeroTol
    hStar = 0.; %large rarefaction
    phi(hMin);
    return;
end
%tic
hStar = fzero( phi, [zeroTol, maxWaterHeight]);
%toc

%h = zeroTol: 0.01:max(hMax, hStar)+10;

%plot(h, phi(h), hMin, phi(hMin), '-o', hMax, phi(hMax), '-o', hStar, 0, '-o')

%text(hMin, phi(hMin), '  h_{min}')
%text(hMax, phi(hMax), '  h_{max}')
%text(hStar, 0, '  h_{star}')

%xlabel('h')
%ylabel('\phi(h)')

%grid on