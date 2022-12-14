clear;
thisFile = "findExt_drive";
%
%%
%%%
for trialIndex=1:1
%for trialIndex=1:1000
%%%
%%
%

%setprngstates();
% These cases were before fixing "-f_n" in "rho_n"...
%setprngstates(40793712); % A very bad case.
%setprngstates(72305920); % A slightly bad case.
%setprngstates(67766576); % Another hard case
%%%setprngstates(43447584); % A nicely converging case.
%setprngstates(97827072); % Needs BT.
% After fix, before 05-25-1640
%setprngstates(64596560); % Very slowly converging case.
%setprngstates(98430512); % Bad case - all pts very far from ext.
%setprngstates(70837456); % Woah! Need multi-interval analysis here!
%setprngstates(74903232); % Very bad initial points.
%setprngstates(3158448); % Bad initial points.
%setprngstates(61627072); % P is 0.0627.
% After 05-25-1640 param changes
%setprngstates(43458496); % Triggers a new path.
%setprngstates(73688096); % P is 0.458 and one of the pts in nearly on top of secret_bigX.
%setprngstates(74718336); % P is 0.0009639.
setprngstates(7254960); % Initial guess is totally wrong; P=0.1245.

numFigs = 0;
%
secret_bigX = randn*exp(2*randn)
secret_bigP = abs(2+randn)
secret_bigF0 = randn
secret_bigF1 = randn
funch_f = @(x)( secret_bigF0 + secret_bigF1 * abs( x - secret_bigX ).^secret_bigP );
%%%funch_f = @(x)( secret_bigF0 + secret_bigF1 * exp( -(x-secret_bigX).^(-2) ) );
%
numPts = 7;
xVals = sort( randn(1,numPts) );
if (1)
	% 05-25-1640 modifications.
	xVals += secret_bigX + randn;
	%xVals(1:2) -= 0.01 % Use with 73688096 to cnvg.
end
%%%xVals = sort( secret_bigX + randn(1,numPts) );
fVals = funch_f(xVals);
%
tic();
findExt; thisFile = "findExt_drive";

if ( res > 1e-20 )
	msg( thisFile, __LINE__, "" );
	msg( thisFile, __LINE__, "!!!!!!!!!!!!!!!!!!!" );
	msg( thisFile, __LINE__, "! FOUND A BAD CASE!" );
	msg( thisFile, __LINE__, "!!!!!!!!!!!!!!!!!!!" );
	msg( thisFile, __LINE__, "" );
	break;
end
%
%%
%%%
end
%%%
%%
%

toc();
msg( thisFile, __LINE__, sprintf( ...
  "BigX:    Initial was %10.3e;   Converged to %10.3e;   Secret is %10.3e;   Res is %10.3e.", ...
  bigX_initial, bigX, secret_bigX, bigX-secret_bigX ) );
msg( thisFile, __LINE__, sprintf( ...
  "BigP:    Initial was %10.3e;   Converged to %10.3e;   Secret is %10.3e;   Res is %10.3e.", ...
  bigP_initial, bigP, secret_bigP, bigP-secret_bigP ) );
msg( thisFile, __LINE__, sprintf( ...
  "BigF0:   Initial was %10.3e;   Converged to %10.3e;   Secret is %10.3e;   Res is %10.3e.", ...
  bigF0_initial, bigF0, secret_bigF0, bigF0-secret_bigF0 ) );
msg( thisFile, __LINE__, sprintf( ...
  "BigF1:   Initial was %10.3e;   Converged to %10.3e;   Secret is %10.3e;   Res is %10.3e.", ...
  bigF1_initial, bigF1, secret_bigF1, bigF1-secret_bigF1 ) );
%
funch_fModel_initial = @(x)( bigF0_initial + bigF1_initial*abs(x-bigX_initial).^bigP_initial );
funch_fModel = @(x)( bigF0 + bigF1*abs(x-bigX).^bigP );
%
viz_numPts = 100;
viz_x = linspace( ...
  min([min(xVals), secret_bigX])-1.0, ...
  max([max(xVals), secret_bigX])+1.0, viz_numPts );
viz_f = funch_f(viz_x);
viz_fModel = funch_fModel(viz_x);
viz_fModel_initial = funch_fModel_initial(viz_x);
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
  bigX, bigF0, '*', 'linewidth', 4, 'markersize', 25, ...
  viz_x, viz_fModel_initial, 'v-', 'linewidth', 1, 'markersize', 10, ...
  bigX_initial, bigF0_initial, 'v', 'linewidth', 4, 'markersize', 25, ...
  xVals, fVals, 'o', 'linewidth', 3, 'markersize', 30 );
grid on;
ax0 = axis();
%
numFigs++; figure(numFigs);
plot( ... 
  cviz_x, cviz_h1, 'o-', 'linewidth', 2, 'markersize', 5, ...
  cviz_x, cviz_h1Model, 's-', 'linewidth', 1, 'markersize', 3, ...
  secret_bigX, 0.0, 'x', 'linewidth', 4, 'markersize', 20, ...
  bigX, 0.0, '*', 'linewidth', 4, 'markersize', 25 );
grid on;
axis([ax0(1),ax0(2),-5.0,5.0]);
%
numFigs++; figure(numFigs);
plot( ... 
  ccviz_x, ccviz_h2, 'o-', 'linewidth', 2, 'markersize', 5, ...
  ccviz_x, ccviz_h2Model, 's-', 'linewidth', 1, 'markersize', 3, ...
  secret_bigX, 0.0, 'x', 'linewidth', 4, 'markersize', 20, ...
  bigX, 0.0, '*', 'linewidth', 4, 'markersize', 2, ...
  super_xVals, super_hVals, 'o', 'markersize', 25, ...
  super_xVals_probablyBad, super_hVals, 'v', 'markersize', 20 );
grid on;
axis([ax0(1),ax0(2),-5.0,5.0]);
