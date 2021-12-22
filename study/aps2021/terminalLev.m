clear;
thisFile = "terminalLev";
tic();
numFigs = 0;
%
ax = [ -0.3, 1.2, -0.3, 1.2 ];
sizeX1 = 100;
sizeX2 = 100;
%
%vecG = [ 1; 0 ];
%matH = [ 1, 0; 0, -1 ];
funchOmega = @(x1,x2) 0.5*(x1.^2 - x2.^2);
multiArgLevel_omega = 2;
[ gridX1, gridX2, gridOmega ] = gridfunch( funchOmega, multiArgLevel_omega, ax, sizeX1, sizeX2 );
%
muCrit = 1.0;
muMin = muCrit + 0.0001;
muMax = 100.0;
numMuVals = 100;
muVals = muMin+(muMax-muMin)*linspace( 0.0, 1.0, numMuVals ).^2;
curve1X10 = 1.0;
curve1X20 = 0.0;
curve2X10 = 1.0;
curve2X20 = 0.01;
curve1X1Vals = curve1X10*muVals./(muVals+1.0);
curve1X2Vals = curve1X20*muVals./(muVals-1.0);
curve2X1Vals = curve2X10*muVals./(muVals+1.0);
curve2X2Vals = curve2X20*muVals./(muVals-1.0);
%
numColors = 51;
cMap = 0.6 + (0.4*jet(numColors));
%
%
numFigs++; figure(numFigs);
hold off;
plot( ...
  curve1X1Vals, curve1X2Vals, 'k-', 'linewidth', 3, ...
  curve2X1Vals, curve2X2Vals, 'r-', 'linewidth', 3 );
axis(ax);
axis equal;
grid on;
title( "Levenberg curves for omega = 0.5 * ( x_1^2 - x_2^2 )" );
xlabel( "x1" );
ylabel( "x2" );
legend( ...
  "LevC from x_1=1.0, x_2=0.0", ...
  "LevC from x_1=1.0, x_2=0.01" );
hold on;
contourf( gridX1, gridX2, gridOmega, 31 );
colormap(cMap);
plot( ...
  curve1X1Vals, curve1X2Vals, 'k-', 'linewidth', 3, ...
  curve2X1Vals, curve2X2Vals, 'r-', 'linewidth', 3 );
plot( ...
  curve1X10, curve1X20, 'ko', 'linewidth', 5, 'markersize',  5, ...
  curve1X10, curve1X20, 'ko', 'linewidth', 5, 'markersize', 10, ...
  curve1X10, curve1X20, 'ko', 'linewidth', 5, 'markersize', 15, ...
  curve1X10/2.0, 0.0, 'ko', 'linewidth', 3, 'markersize', 20 );
plot( 
  curve2X10, curve2X20, 'ro', 'linewidth', 5, 'markersize',  5, ...
  curve2X10, curve2X20, 'ro', 'linewidth', 5, 'markersize', 10, ...
  curve2X10, curve2X20, 'ro', 'linewidth', 5, 'markersize', 15 );
hold off




%
sMin = 0.0001;
numSVals = 100;
sVals = sMin+(1.0-sMin)*linspace( 0.0, 1.0, numSVals );
curve1X10 = 1.0;
curve1X20 = 0.0;
curve2X10 = 1.0;
curve2X20 = 0.01;
curve1X1Vals = curve1X10.*sVals;
curve1X2Vals = curve1X20./sVals;
curve2X1Vals = curve2X10.*sVals;
curve2X2Vals = curve2X20./sVals;
%
%
numFigs++; figure(numFigs);
hold off;
plot( ...
  curve1X1Vals, curve1X2Vals, 'k-', 'linewidth', 3, ...
  curve2X1Vals, curve2X2Vals, 'r-', 'linewidth', 3 );
axis(ax);
axis equal;
grid on;
title( "Gradient curves for omega = 0.5 * ( x_1^2 - x_2^2 )" );
xlabel( "x1" );
ylabel( "x2" );
legend( ...
  "GradC from x_1=1.0, x_2=0.0", ...
  "GradC from x_1=1.0, x_2=0.01" );
hold on;
contourf( gridX1, gridX2, gridOmega, 31 );
colormap(cMap);
plot( ...
  curve1X1Vals, curve1X2Vals, 'k-', 'linewidth', 3, ...
  curve2X1Vals, curve2X2Vals, 'r-', 'linewidth', 3 );
plot( ...
  curve1X10, curve1X20, 'ko', 'linewidth', 5, 'markersize',  5, ...
  curve1X10, curve1X20, 'ko', 'linewidth', 5, 'markersize', 10, ...
  curve1X10, curve1X20, 'ko', 'linewidth', 5, 'markersize', 15, ...
  0.0, 0.0, 'ko', 'linewidth', 3, 'markersize', 5, ...
  0.0, 0.0, 'ko', 'linewidth', 3, 'markersize', 10, ...
  0.0, 0.0, 'ko', 'linewidth', 3, 'markersize', 15 );
plot( 
  curve2X10, curve2X20, 'ro', 'linewidth', 5, 'markersize',  5, ...
  curve2X10, curve2X20, 'ro', 'linewidth', 5, 'markersize', 10, ...
  curve2X10, curve2X20, 'ro', 'linewidth', 5, 'markersize', 15 );
hold off
