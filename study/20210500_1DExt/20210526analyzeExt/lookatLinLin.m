clear;
commondefs;
thisFile = "lookatLinLin";
setprngstates(3);
numFigs = 0;

epsFD = sqrt(eps);

bigX = 0.0;
bigP = 20.0;
bigA = 1.0;
bigB = abs(randn());
funchF = @(x)( bigA + bigB*abs(x-bigX).^bigP );

%xVals = sort(randn(1,3));
xVals = [-1, 0.3, 0.4];
xVals = [-1, 0.7, 0.8];
fVals = funchF(xVals);

vecX = xVals';
vecF = fVals';
matX = [ ones(3,1), vecX, vecX.^2 ];
vecC = matX\vecF;
funchG = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );

funchDF = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDG = @(x)( (funchG(x+epsFD)-funchG(x-epsFD))/(2.0*epsFD) );

xQuad = -vecC(2)/(2.0*vecC(3))
dfVals = diff(fVals)./diff(xVals);
cxVals = cent(xVals);

viz_numPts = 1000;
viz_xVals = linspace(min(xVals),max(xVals),viz_numPts);
viz_fVals = funchF(viz_xVals);
viz_gVals = funchG(viz_xVals);
viz_dfVals = funchDF(viz_xVals);
viz_dgVals = funchDG(viz_xVals);

numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_fVals, 'o-', ...
  viz_xVals, viz_gVals, 'x-', ...
  xVals, fVals, 'k*', 'linewidth', 4, 'markersize', 25, ...
  bigX, funchF(bigX), 'ks', 'linewidth', 4, 'markersize', 25, ...
  xQuad, funchG(xQuad), 'r^', 'linewidth', 4, 'markersize', 25 );
grid on;

numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_dfVals, 'o-', ...
  viz_xVals, viz_dgVals, 'x-', ...
  cxVals, dfVals, 'cv-', 'linewidth', 10, 'markersize', 40, ...
  bigX, 0.0, 'ks', 'linewidth', 4, 'markersize', 25, ...
  xQuad, 0.0, 'r^', 'linewidth', 4, 'markersize', 25 );
grid on;
