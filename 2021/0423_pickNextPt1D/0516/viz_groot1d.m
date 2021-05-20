thisFile = "viz_groot1d";
numFigs = 0;

xVals = datOut.xVals_raw;
fVals = datOut.fVals_raw;
numPts = datOut.numPts;
evalIndex = (1:numPts);

xVizBumper = (abs(x2-x1) + max(xVals) - min(xVals))/5.0;
xVizLo = min(xVals) - xVizBumper;
xVizHi = max(xVals) + xVizBumper;
xVals_viz = linspace( xVizLo, xVizHi, 1001 );
fVals_viz = funchF( xVals_viz );
%
%
numFigs++; figure(numFigs);
plot( ...
  xVals(3:end), fVals(3:end), 'o-', 'markersize', 10, 'linewidth', 2, ...
  xVals(1:2), fVals(1:2), '+-', 'markersize', 20, 'linewidth', 3, ...
  xVals(end), fVals(end), 'x', 'markersize', 25, 'linewidth', 3, ...
  xVals(2:3), fVals(2:3), '-', 'linewidth', 2, ...
  xVals_viz, fVals_viz, 'ko-', 'linewidth', 1, 'markersize', 2 );
title("F vs X");
xlabel("X");
ylabel("F");
fVizBumper = (max(fVals)-min(fVals))/2.0;
fVizLo = 0.0;
fVizHi = 0.0;
if ( min(fVals) < 0.0 )
	fVizLo = min(fVals)-fVizBumper;
end
if ( max(fVals) > 0.0 )
	fVizHi = max(fVals)+fVizBumper;
end
axis_fvx = axis();
axis([ axis_fvx(1), axis_fvx(2), fVizLo, fVizHi ]);
axis_fvx = axis();
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  xVals(3:end), -evalIndex(3:end), 'o-', 'markersize', 10, 'linewidth', 2, ...
  xVals(1:2), -evalIndex(1:2), '+-', 'markersize', 20, 'linewidth', 3, ...
  xVals(end), -evalIndex(end), 'x', 'markersize', 25, 'linewidth', 3, ...
  xVals(2:3), -evalIndex(2:3), '-', 'linewidth', 2 );
axis([ axis_fvx(1), axis_fvx(2), -numPts-1, 0 ]);
title("X Drip");
xlabel("X");
ylabel("Neg Eval Index");
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  evalIndex(3:end), fVals(3:end), 'o-', 'markersize', 10, 'linewidth', 2, ...
  evalIndex(1:2), fVals(1:2), '+-', 'markersize', 20, 'linewidth', 3, ...
  evalIndex(end), fVals(end), 'x', 'markersize', 25, 'linewidth', 3, ...
  evalIndex(2:3), fVals(2:3), '-', 'linewidth', 2 );
title("F vs Feval");
xlabel("Feval");
ylabel("F");
grid on;
%
numFigs++; figure(numFigs);
semilogy( ...
  evalIndex(1:end-1), abs(xVals(1:end-1)-xVals(end)), 's-', 'linewidth', 2,...
  evalIndex, abs(fVals), 'x-', 'linewidth', 2 );
title("Cnvg vs Feval");
xlabel("Feval");
legend( ...
  "|X-X_{final}|", ...
  "|F|", ...
  "location", "northeast" );
grid on;

numFigs++; figure(numFigs);
semilogy( ...
  0.5+evalIndex(1:end-1), abs(diff(xVals)), 's-', 'linewidth', 2, ...
  0.5+evalIndex(1:end-1), abs(diff(fVals)), '+-', 'linewidth', 3, ...
  0.5+evalIndex(1:end-1), abs(diff(abs(fVals))), 'x-', 'linewidth', 2 );
title("Delta vs Feval");
xlabel("Feval");
legend( ...
  "|X_{n+1} - X_{n}|", ...
  "|F_{n+1} - F_{n}|", ...
  "| |F_{n+1}| - |F_{n}| |", ...
  "location", "northeast" );
grid on;
