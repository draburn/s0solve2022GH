clear;
%setprngstates();
setprngstates(40793712);
numFigs = 0;
%
secret_bigX = randn*exp(2*randn)
secret_bigP = abs(2+randn)
%secret_bigP = 2.0
secret_bigF0 = randn
secret_bigF1 = randn
funch_f = @(x)( secret_bigF0 + secret_bigF1 * abs( x - secret_bigX ).^secret_bigP );
%
numPts = 7;
xVals = sort( randn(1,numPts) );
fVals = funch_f(xVals);
%
findExt
echo__bigX = bigX
echo__bigP = bigP
echo__bigF0 = bigF0
echo__bigF1 = bigF1
funch_fModel = @(x)( bigF0 + bigF1*abs(x-bigX).^bigP );
%
viz_numPts = 100;
viz_x = linspace( ...
  min([min(xVals), secret_bigX])-1.0, ...
  max([max(xVals), secret_bigX])+1.0, viz_numPts );
viz_f = funch_f(viz_x);
viz_fModel = funch_fModel(viz_x);
%
cviz_x = cent(viz_x);
ccviz_x = cent(cviz_x);
%
cviz_f = cent(viz_f);
dviz_f = diff(viz_f)./diff(viz_x);
cdviz_f = cent(dviz_f);
ddviz_f = diff(dviz_f)./diff(cviz_x);
cviz_h1 = cviz_f./dviz_f;
ccviz_h2 = cdviz_f./ddviz_f;
%
cviz_fModel = cent(viz_fModel);
dviz_fModel = diff(viz_fModel)./diff(viz_x);
cdviz_fModel = cent(dviz_fModel);
ddviz_fModel = diff(dviz_fModel)./diff(cviz_x);
cviz_h1Model = cviz_fModel./dviz_fModel;
ccviz_h2Model = cdviz_fModel./ddviz_fModel;
%
numFigs++; figure(numFigs);
plot( ...
  viz_x, viz_f, 'x-', 'linewidth', 2, 'markersize', 5, ...
  viz_x, viz_fModel, '+-', 'linewidth', 1, 'markersize', 15, ...
  secret_bigX, secret_bigF0, 's', 'linewidth', 4, 'markersize', 20, ...
  xVals, fVals, 'o', 'linewidth', 3, 'markersize', 30 );
grid on;
ax0 = axis();
%
numFigs++; figure(numFigs);
plot( ... 
  cviz_x, cviz_h1, 'o-', 'linewidth', 2, 'markersize', 5, ...
  cviz_x, cviz_h1Model, 's-', 'linewidth', 1, 'markersize', 3, ...
  secret_bigX, 0.0, 'x', 'linewidth', 4, 'markersize', 20 );
grid on;
axis([ax0(1),ax0(2),-5.0,5.0]);
%
numFigs++; figure(numFigs);
plot( ... 
  ccviz_x, ccviz_h2, 'o-', 'linewidth', 2, 'markersize', 5, ...
  ccviz_x, ccviz_h2Model, 's-', 'linewidth', 1, 'markersize', 3, ...
  secret_bigX, 0.0, 'x', 'linewidth', 4, 'markersize', 20 );
grid on;
axis([ax0(1),ax0(2),-5.0,5.0]);
