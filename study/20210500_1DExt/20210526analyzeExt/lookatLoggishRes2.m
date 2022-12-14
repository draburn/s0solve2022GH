clear;
commondefs;
thisFile = "lookatLoggishRes2";
numFigs = 0;
%
%
caseNum = 0;
switch (caseNum)
case 0
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
	x_secret = -0.890874387338734; % This is the left abs min.
otherwise
	error(["Invalid value of caseNum (", num2str(caseNum), ")."]);
end
%
%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.6, 1.2 ])
%xVals = sort([ -1.00057993392590, -1.0, 0.0 ])
%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ])
% Add balance.
%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, -0.831698038563058, 0.0 ])
% Iterate
%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, -0.897832, -0.831698038563058, 0.0 ])
% Iterate
%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, -0.897832, -0.868436, -0.831698038563058, 0.0 ])
% And... Here, this new tech does something different! -0.883699
%xVals = sort([ -0.883699, -1.00057993392590, -1.0, -0.971949673093843, -0.897832, -0.883145836017872, -0.868436, -0.831698038563058, 0.0 ])
% -0.883555201402426
%xVals = sort([ -0.883555201402426, xVals ]);
%xVals = sort([ -0.884206, xVals ]);
%xVals = sort([ -0.884889, xVals ]);
%
%
%xVals = sort([ -0.883555201402426, -0.883699, -1.00057993392590, -1.0, -0.971949673093843, -0.897832, -0.883145836017872, -0.868436, -0.831698038563058])
%xVals = sort([ -0.886, -0.880, xVals ]); % Force an absMax.
%
% Lets re-do...
xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ]);
xVals = sort([ -0.831698038563058, xVals ]); %xbalanceR
xVals = sort([ -0.897832, xVals ]);
xVals = sort([ -0.868437, xVals ]);
xVals = sort([ -0.883142, xVals ]);
xVals = sort([ -0.883701, xVals ]);
xVals = sort([ -0.886496000000000, xVals ]); %xBalanceL
xVals = sort([ -0.886247, xVals ]);
xVals = sort([ -0.887741000000000, xVals ]);
xVals = sort([ -0.888869, xVals ]);
xVals = sort([ -0.889642, xVals ]);
xVals = sort([ -0.893507000000000, xVals ]);
xVals = sort([ -0.8890615, xVals ]);
xVals = sort([ -0.8906340, xVals ]);
%xVals = sort([ -0.8906723, xVals ]);
%xVals = sort([ -0.89086380000000, xVals ]); % ptweAbsF =    8.81227741557147e-13
%
% Try with avg...
xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ]);
xVals = sort([ -0.831698038563058, xVals ]);
xVals = sort([ -0.897620231966365, xVals ]);
xVals = sort([ -0.867992251154230, xVals ]);
xVals = sort([ -0.88284017185230, xVals ]);
xVals = sort([ -0.887137697860992, xVals ]); % xB goes past data.
xVals = sort([ -0.886078607798801, xVals ]);
xVals = sort([ -0.890138678468405, xVals ]); % xB goes way past data!
xVals = sort([ -0.889487776317036, xVals ]); % must balance...
xVals = sort([ -0.893393189225250, xVals ]);
%xVals = sort([ -0.890922389166615, xVals ]); % xB goes past...
%xVals = sort([ -0.890775831691018, xVals ]); % must balance...
%xVals = sort([ -0.891655176544600, xVals ]);
%xVals = sort([ -0.890874186962029, xVals ]); %1E-16.
% Allowing going past by a little bit, as long as avg is well-behaved.
%xVals = sort([ , xVals ]);
%
msg( thisFile, __LINE__, "THIS IS CONVERGING TO THAT LOCAL MAX!!!" );
msg( thisFile, __LINE__, "CHECK FOR PTWISE LOCAL MAX AND CHOP TO PREVENT THIS!" );
msg( thisFile, __LINE__, "Even so, it may be easier to find an ext than identify whether max or min." );
numPts = size(xVals,2);
fVals = funchF(xVals)

% Confirm no ptwise absMax
for n=2:numPts-1
if ( abs(fVals(n)) >= abs(fVals(n-1)) ...
  && abs(fVals(n)) >= abs(fVals(n+1))  )
	msg( thisFile, __LINE__, sprintf( ...
	  "Found a local max at %d!", n ) );
	error( "Found a local max." );
end
end

% point-wise extremum
[ ptweAbsF, ptweIndex ] = min( abs(fVals) )
assert( 2 <= ptweIndex );
assert( numPts-1 >= ptweIndex );
ptweX = xVals(ptweIndex);
%
n = ptweIndex;
if ( xVals(n-1) < xVals(n) - 10.0*( xVals(n+1) - xVals(n) ) )
	xBalanceL = xVals(n) - 5.0*( xVals(n+1) - xVals(n) )
	error( "Don't forget to use xBalanceL!" );
	return
end
if ( xVals(n+1) > xVals(n) + 10.0*( xVals(n) - xVals(n-1) ) )
	xBalanceR = xVals(n) + 5.0*( xVals(n) - xVals(n-1) )
	error( "Don't forget to use xBalanceR!" );
	return
end
%
vecXA = xVals(ptweIndex-1:ptweIndex+1)';
vecFA = fVals(ptweIndex-1:ptweIndex+1)';
matXA = [ ones(3,1), vecXA, vecXA.^2 ];
vecCA = matXA\vecFA;
funchGA = @(x)( vecCA(1) + vecCA(2)*x + vecCA(3)*x.^2 );
assert( 0.0 <= vecCA(3) );
%
xExtA = -vecCA(2)/(2.0*vecCA(3))
msg( thisFile, __LINE__, ...
  "Here we get -0.979251248665930, but test_groot1d gives -0.958662835519853?" );
if ( xExtA > ptweX )
	indexB = ptweIndex+1
else
	indexB = ptweIndex-1
end
assert( 2 <= indexB );
assert( numPts-1 >= indexB );
%
vecXB = xVals(indexB-1:indexB+1)';
vecFB = fVals(indexB-1:indexB+1)';
matXB = [ ones(3,1), vecXB, vecXB.^2 ];
vecCB = matXB\vecFB;
funchGB   = @(x)( vecCB(1) + vecCB(2)*x + vecCB(3)*x.^2 );
xExtB = -vecCB(2)/(2.0*vecCB(3))
if ( vecCB(3) < 0.0 )
	msg( thisFile, __LINE__, "Model B has wrong curvature." );
elseif ( xExtB < vecXB(1) )
	msg( thisFile, __LINE__, "extB is past left data point." );
elseif ( xExtB > vecXB(3) )
	msg( thisFile, __LINE__, "extB is past right data point." );
end
xAvg = (xExtA+xExtB)/2.0

if (0)
	vecXC = xVals(4:7)';
	vecFC = fVals(4:7)';
	matXC = [ ones(4,1), vecXC, vecXC.^2 ];
	vecCC = matXC\vecFC;
	%matWC = diag( 1.0./sqrt(sqrt(abs(vecFC))) );
	%vecCC = (matWC*matXC)\(matWC*vecFC);
	funchGC = @(x)( vecCC(1) + vecCC(2)*x + vecCC(3)*x.^2 );
	x4ptQuad = -vecCC(2)/(2.0*vecCC(3))
	%
	% two-plus-two quad, approximately
	vecXD = xVals(4:7)';
	vecFD = fVals(4:7)';
	matXD = [ ones(4,1), vecXD, vecXD.^2 ];
	matWD = diag([ 1.0, 1e8, 1e8, 1.0 ]');
	%matWD = diag([ 1.0, 1e8, 1e8, 1.0 ]'./sqrt(abs(vecFD)));
	vecCD = (matWD*matXD)\(matWD*vecFD);
	funchGD = @(x)( vecCD(1) + vecCD(2)*x + vecCD(3)*x.^2 );
	xApprox2p2Quad = -vecCD(2)/(2.0*vecCD(3))
	%
	funchF([ xAvg, x4ptQuad, xApprox2p2Quad ])
	return;
end



epsFD = sqrt(sqrt(eps));
funchDF = @(x)( (funchF(x+epsFD)-funchF(x-epsFD))/(2.0*epsFD) );
funchDDF = @(x)( (funchF(x+epsFD)+funchF(x-epsFD)-2*funchF(x))/(epsFD^2) );
funchLoggyF = @(x)( funchDF(x)./abs(funchDDF(x)) );
funchDGA = @(x)( vecCA(2) + 2.0*vecCA(3)*x );
funchDGB = @(x)( vecCB(2) + 2.0*vecCB(3)*x );
funchDDGA = @(x)( 2.0*vecCA(3)*ones(size(x)) );
funchDDGB = @(x)( 2.0*vecCB(3)*ones(size(x)) );
funchLoggyGA = @(x)( funchDGA(x)./abs(funchDDGA(x)) );
funchLoggyGB = @(x)( funchDGB(x)./abs(funchDDGB(x)) );



viz_numPts = 5000;
viz_xVals = linspace(-1.6,1.6,viz_numPts);
%viz_xVals = linspace(-1.0,-0.8,viz_numPts);
viz_xAVals = linspace(min(vecXA),max(vecXA),viz_numPts);
viz_xBVals = linspace(min(vecXB),max(vecXB),viz_numPts);
%
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, funchF(viz_xVals), 'o-', ...
  viz_xVals, funchGA(viz_xVals), '^-', ...
  viz_xVals, funchGB(viz_xVals), 'v-', ...
  xVals, fVals, 'ko', 'linewidth', 4, 'markersize', 25, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
%
capLo = -1.0;
capHi =  1.0;
xCrossVals = [ vecXB(2), vecXA(2) ];
%xCrossVals = [ (vecXB(1)+vecXB(2))/2.0, (vecXA(2)+vecXA(3))/2.0 ];
%xCrossVals = [ (vecXB(1)+vecXB(2)+vecXB(3))/2.0, (vecXA(2)+vecXA(3)+vecXA(3))/3.0 ];
numFigs++; figure(numFigs);
plot( ...
  viz_xVals, cap(funchLoggyF(viz_xVals),capLo,capHi), 'o-', ...
  viz_xAVals, funchLoggyGA(viz_xAVals), '^-', ...
  viz_xBVals, funchLoggyGB(viz_xBVals), 'v-', ...
  vecXA, funchLoggyGA(vecXA), '^', 'markersize', 25, 'linewidth', 4, ...
  vecXB, funchLoggyGB(vecXB), 'v', 'markersize', 25, 'linewidth', 4, ...
  xCrossVals, [funchLoggyGB(xCrossVals(1)),funchLoggyGA(xCrossVals(2))], 's-', 'markersize', 25, 'linewidth', 4, ...
  xAvg, 0.0, 'cx', 'markersize', 30, 'linewidth', 3, ...
  x_secret, 0.0, 'y*', 'markersize', 30, 'linewidth', 3, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f'/|f''|)" );
title( "f'/|f''| vs x" );
axis([-0.9 -0.86 -0.1 0.1]);
