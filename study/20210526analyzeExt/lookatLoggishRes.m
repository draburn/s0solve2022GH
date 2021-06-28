clear;
commondefs;
thisFile = "lookatLoggishRes";
%setprngstates(3);
setprngstates(33754416);
numFigs = 0;
%
modelOrder = 2;
epsFD = sqrt(sqrt(eps));
vecP = randn(modelOrder+2,1);
vecP(end) *= 0.1;
switch (modelOrder)
case 0
	funchF = @(x)( vecP(1) + vecP(2)*x );
case 1
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 );
case 2
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 + vecP(4)*x.^3 );
case 3
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 + vecP(4)*x.^3 + vecP(5)*x.^4 );
case 4
	funchF = @(x)( vecP(1) + vecP(2)*x + vecP(3)*x.^2 + vecP(4)*x.^3 + vecP(5)*x.^4 + vecP(6)*x.^5 );
otherwise
	error(["Invalid value of modelOrder (", num2str(modelOrder), ")."]);
end

if (1)
	msg( thisFile, __LINE__, "OVER-RIDING FUNCTION." );
	%funchF = @(x)( exp(-1.0./(x.^2)) );
	%funchF = @(x)( 1E-4 + abs(x-0.0).^20.5 );
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
end

funchDF  = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );
%
numPts = modelOrder+1;
xVals = sort(randn(1,numPts));
if (1)
	msg( thisFile, __LINE__, "OVER-RIDING X-VALUES." );
	%xVals = [-2.1288  -1.1109   1.4579];
	% Zero is at 0.173.
	% [-2.1288  -1.1109   0.173    1.4579] = [ 5.3305e+06   8.6366e+00   1.0000e-04   2.2721e+03 ];
	% ptwise min is 3. drop 1.
	%xVals = [-1.1109    0.173   1.4579];
	% Zero at -0.464.
	%[-1.1109  -0.464   0.173   1.4579] = [8.6366e+00   1.0015e-04   1.0000e-04   2.2721e+03];
	% ptwise min is 3. drop 1.
	%xVals = [-0.464   0.173   1.4579];
	% zero at -0.145
	%[-0.464   -0.145  0.173   1.4579] = [ 1.0015e-04   1.0000e-04   1.0000e-04   2.2721e+03];
	% = 1.00145753368791e-04   1.00000000000006e-04   1.00000000000240e-04   2.27205377617069e+03
	% ptwise min is 2. drop 4.
	%xVals = [-0.464   -0.145  0.173 ];
	% zero at 0.014;
	% [-0.464   -0.145  0.014  0.173 ] =
	% 1.00145753368791e-04   1.00000000000006e-04   1.00000000000000e-04   1.00000000000240e-04
	% 2 is ptwise min.
	% And, it's pretty clear the ext val is 1E-4.
	% And, this happens simply with quad model?!?!
	%
	%
	%xVals = [-1.5,-1.0,-0.5]; %Next is -0.85 as pt 3. min is 3. Drop 1.
	%xVals = [-1.0,-0.85,-0.5]; %Next is -0.917 as pt 2. min is 3. Drop 1.
	%xVals = [-0.917,-0.85,-0.5]; %Next is -0.8835 as pt 2. min is 2. Drop 4.
	%xVals = [-0.917, -0.8835, -0.85]; %Next is -0.8823 as pt 3. min is 3. Drop 1.
	%xVals = [-0.8835, -0.8823, -0.85]; %Next is -0.8828...
	%%%
	% Originally looked like bad pts,
	% with pt on right being continually used.
	% Maybe this only happens some of the time?
	%xVals = [ -0.9, -0.85, -0.5]; %-0.875 as 2. Below tol.
	% Try again...
	%xVals = [ -1.4, -0.82, -0.8 ]; %-0.81! On wrong side!
	%xVals = [ -1.4, -0.82, -0.815 ]; %-0.821! Barely on other side!
	% Look at this case with correct points!
	%
	%
	%
	%xVals = [ -1.0, -0.971949673093843, 0.0 ]
	%xVals = [ -0.971949673093843, -0.956782, 0.0 ]
	%xVals = [ -0.956782, -0.945902, 0.0 ]
	%
	% With 10x8 bracketing...
	%xVals = [ -1.0, -0.971949673093843, 0.0 ]
	%xBracketR = -0.747547057844587
	%fBracketR =  0.00940338904975590
	%xVals = [ -1.0, -0.971949673093843, -0.747547057844587 ]
	%xQuadterp = -0.879145655368713
	%fQuadterp =    5.49162809593981e-08
	%xVals = [ -0.971949673093843, -0.879145655368713, -0.747547057844587 ]fVals =
	%   2.31234543096845e-03   5.49162809591900e-08   9.40338904975573e-03
	%xQuadterp = -0.896538890816371
	%fQuadterp =    4.80236286887677e-07
	%xVals = [ -0.896538890816371, -0.879145655368713, -0.747547057844587 ]
	%   4.80236286887677e-07   5.49162809591900e-08   9.40338904975573e-03
	%xQuadterp = -0.887816787854369
	%fQuadterp =    4.67253217977081e-08
	%xVals = [ -0.896538890816371, -0.887816787854369, -0.879145655368713 ]
	%   4.80236286887677e-07   4.67253217975161e-08   5.49162809591900e-08
	%xQuadterp = -0.883643422750339
	%fQuadterp =    1.11644986491803e-07
	% Cnvg seems to have slowed...
end

fVals = funchF(xVals)


%2021.06.27: Use lin-interp when unbalanced?
xLinterp = ( xVals(2)*fVals(1) - xVals(1)*fVals(2) ) / ( fVals(1) - fVals(2) )
fLinterp = funchF(xLinterp)

%
vecX = xVals';
matX = ones(numPts,1);
for m=1:modelOrder
	matX = [ matX, vecX.^m ];
end
vecF = fVals';
vecC = matX\vecF;
switch (modelOrder)
case 0
	funchFModel = @(x)( vecC(1) );
case 1
	funchFModel = @(x)( vecC(1) + vecC(2)*x );
case 2
	funchFModel = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 );
case 3
	funchFModel = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 + vecC(4)*x.^3 );
case 4
	funchFModel = @(x)( vecC(1) + vecC(2)*x + vecC(3)*x.^2 + vecC(4)*x.^3 + vecC(5)*x.^4 );
otherwise
	error(["Invalid value of modelOrder (", num2str(modelOrder), ")."]);
end
assert( 2==modelOrder )
xQuadterp = -vecC(2)/(2.0*vecC(3))
fQuadterp = funchF(xQuadterp)
funchDFModel  = @(x)( (funchFModel(x+epsFD)-funchFModel(x-epsFD))/(2.0*epsFD) );
funchDDFModel = @(x)( (funchFModel(x+epsFD)+funchFModel(x-epsFD)-2*funchFModel(x))/(epsFD^2) );
%
viz_numPts = 5000;
viz_xVals = linspace(-1.6,1.6,viz_numPts);
%viz_xVals = linspace(min(xVals),max(xVals),viz_numPts);
%viz_xVals = linspace(-0.90,-0.87,viz_numPts);
viz_fVals   = funchF(   viz_xVals );
viz_dfVals  = funchDF(  viz_xVals );
viz_ddfVals = funchDDF( viz_xVals );
viz_fModelVals   = funchFModel(   viz_xVals );
viz_dfModelVals  = funchDFModel(  viz_xVals );
viz_ddfModelVals = funchDDFModel( viz_xVals );
%
xAvg = sum(xVals)/numPts;
xMatchDDFforMO2 = xAvg;
xSqAvg = sum(xVals.^2)/numPts;
xVarSq = xSqAvg - (xAvg^2);
assert( xVarSq > 0.0 );
xVar = sqrt(xVarSq);
xMatchDFforMO2m = xAvg - xVar/sqrt(2.0);
xMatchDFforMO2p = xAvg + xVar/sqrt(2.0);
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_fVals, '^-', ...
  viz_xVals, viz_fModelVals, 'v-', ...
  xVals, fVals, 'ko', 'linewidth', 4, 'markersize', 25, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_dfVals, '^-', ...
  viz_xVals, viz_dfModelVals, 'v-', ...
  xMatchDFforMO2m, funchDFModel(xMatchDFforMO2m), 'cs', 'linewidth', 4, 'markersize', 25,
  xMatchDFforMO2p, funchDFModel(xMatchDFforMO2p), 'cs', 'linewidth', 4, 'markersize', 25, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f'" );
title( "f' vs x" );
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, viz_ddfVals, '^-', ...
  viz_xVals, viz_ddfModelVals, 'v-', ...
  xMatchDDFforMO2, funchDDFModel(xMatchDDFforMO2), 'cs', 'linewidth', 4, 'markersize', 25, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f''" );
title( "f'' vs x" );
%
numFigs++; figure(numFigs);
%plot( ...
%  viz_xVals, viz_dfVals./(0.001*max(abs(viz_dfVals)) + abs(viz_ddfVals)), '^-', ...
%  viz_xVals, viz_dfModelVals./(0.001*max(abs(viz_dfModelVals)) + abs(viz_ddfModelVals)), 'v-' );
plot( ...
  viz_xVals, viz_dfVals./(0*max(abs(viz_dfVals)) + abs(viz_ddfVals)), '^-', ...
  viz_xVals, viz_dfModelVals./(0*max(abs(viz_dfModelVals)) + abs(viz_ddfModelVals)), 'v-', ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f'/|f''|" );
title( "f'/|f''| vs x" );
ax0 = axis();
ax0(3) = max([ax0(3),-2.0]);
ax0(4) = min([ax0(4),+2.0]);
axis(ax0);

% Look at simple bracket-balancing...
% 10x3 happens to hit.
% Use 10x8 instead.
if ( xVals(1) < xVals(2) - 10.0*(xVals(3)-xVals(2)) )
	xBracketL = xVals(2) - 8.0*(xVals(3)-xVals(2))
	fBracketL = funchF(xBracketL)
end
if ( xVals(3) > xVals(2) + 10.0*(xVals(2)-xVals(1)) )
	xBracketR = xVals(2) + 8.0*(xVals(2)-xVals(1))
	fBracketR = funchF(xBracketR)
end


% Try f'/|f''|...
%right_x = sum(xVals)/3.0;
%right_x = (xVals(2)+xVals(3))/2.0;
right_x = xVals(3);
right_f = vecC(1) + vecC(2)*right_x + vecC(3)*(right_x^2);
right_df = vecC(2) + 2.0*vecC(3)*right_x;
right_absddf = abs(2.0*vecC(3));
right_dfoabsddf = right_df / right_absddf;
%
%
%
%zzz
%xVals = [ -1.00057993392590, -1.0, -0.971949673093843 ]
%xVals = [ -0.971949673093843, -0.956782, 0.0 ]
%xVals = [ -0.956782, -0.945902, 0.0 ]
%
% With 10x8 bracketing...
%xVals = [ -1.00057993392590, -1.0, -0.971949673093843 ]
xVals = [ -1.0, -0.971949673093843, -0.747547057844587 ]
fVals = funchF(xVals)
vecX = xVals';
vecF = fVals';
matX = [ ones(numPts,1), vecX, vecX.^2 ];
vecC = matX\vecF;
%
%left_x = sum(xVals)/3.0;
%left_x = (xVals(1)+xVals(2))/2.0;
left_x = xVals(1);
left_f = vecC(1) + vecC(2)*left_x + vecC(3)*(left_x^2);
left_df = vecC(2) + 2.0*vecC(3)*left_x;
left_absddf = abs(2.0*vecC(3));
left_dfoabsddf = left_df / left_absddf;
plot( ...
  [left_x, right_x], [left_dfoabsddf, right_dfoabsddf], 'o-', ...
  [left_x, right_x], 0*[left_dfoabsddf, right_dfoabsddf], 'k-' );
grid on;
%
return;

-0.956782
-0.945902

echo__xNew = 0
echo__xNew =  1
echo__xNew = -1
echo__xNew = -1.00057993392590
echo__xNew = -0.971949673093843
echo__xNew = -0.958662835519853
echo__xNew = -0.944751664283285
echo__xNew = -0.934247587337200
echo__xNew = -0.938570872594207
echo__xNew = -0.923499963559985
echo__xNew = -0.928399100635113
echo__xNew = -0.915577850979833
echo__xNew = -0.919312681802154
echo__xNew = -0.917260427655525
echo__xNew = -0.908624951405276
echo__xNew = -0.911993611185633
echo__xNew = -0.910222926388977
echo__xNew = -0.903456164678508
echo__xNew = -0.905991121284682
echo__xNew = -0.904683604616406
echo__xNew = -0.904033937667877
echo__xNew = -0.899355165070689
echo__xNew = -0.901382906799019
echo__xNew = -0.900350758059296
echo__xNew = -0.899836638519024
echo__xNew = -0.896409944246665
echo__xNew = -0.897872192993241
echo__xNew = -0.897132678927650
echo__xNew = -0.896763794021399
echo__xNew = -0.894347000734750
echo__xNew = -0.895373638126355
echo__xNew = -0.894856375194142
echo__xNew = -0.894598143690240
echo__xNew = -0.894469217012304
echo__xNew = -0.892924361057104
echo__xNew = -0.893633389849797
echo__xNew = -0.893277004116489
echo__xNew = -0.893099003045993
echo__xNew = -0.893010094007102
echo__xNew = -0.892034616778702
echo__xNew = -0.892478390791841
echo__xNew = -0.892255599034039
echo__xNew = -0.892144293615096
echo__xNew = -0.892088684600178
echo__xNew = -0.891508263223970
echo__xNew = -0.891770900486891
echo__xNew = -0.891639134303237
echo__xNew = -0.891573295183846
echo__xNew = -0.891540397109713
echo__xNew = -0.891523958687692
echo__xNew = -0.891207843116097
echo__xNew = -0.891357784545455
echo__xNew = -0.891282590733574
echo__xNew = -0.891245015966321
echo__xNew = -0.891226239494477
echo__xNew = -0.891216856672398
echo__xNew = -0.891046979095507
echo__xNew = -0.891127276918810
echo__xNew = -0.891087016207936
echo__xNew = -0.891066896870080
echo__xNew = -0.891056842665360
echo__xNew = -0.891051818282005
echo__xNew = -0.890962873008389
echo__xNew = -0.891004858534046
echo__xNew = -0.890983809418775
echo__xNew = -0.890973290395576
echo__xNew = -0.890968033639927
echo__xNew = -0.890965406630694
echo__xNew = -0.890919519691522
echo__xNew = -0.890941162246764
echo__xNew = -0.890930312477630
echo__xNew = -0.890924890386046
echo__xNew = -0.890922180735251
echo__xNew = -0.890920826608169
echo__xNew = -0.890897344959103
echo__xNew = -0.890908415064481
echo__xNew = -0.890902865582667
echo__xNew = -0.890900092256819
echo__xNew = -0.890898706296490
echo__xNew = -0.890898013660549
echo__xNew = -0.890886048936731
echo__xNew = -0.890891688203020
echo__xNew = -0.890888861250873
echo__xNew = -0.890887448488647
echo__xNew = -0.890886742476340
echo__xNew = -0.890886389649031
echo__xNew = -0.890880306749172
