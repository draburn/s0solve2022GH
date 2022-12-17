printf("\n\n\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
numFigs = 0;
msg( __FILE__, __LINE__, "Hai!" );
matX = [ datOut.vecXBest, datOut.matX ];
matG = [ datOut.vecGBest, datOut.matG ];
rvecF = [ datOut.fBest, datOut.rvecF ];
sizeX = size(matX,1)
numRecsWB = size(matX,2)
cPrm = [];
cPrm.funchFG = funchFG;
sxPrm = datOut.prm;
%
msg( __FILE__, __LINE__, "Starting from this point..." );
%%%bwdIndex = 27
bwdIndex = 28
matX = matX(:,end-bwdIndex:end);
matG = matG(:,end-bwdIndex:end);
rvecF = rvecF(end-bwdIndex:end);
%
[ fBest, indexOfBest ] = min(rvecF);
vecXBest = matX(:,indexOfBest);
vecGBest = matG(:,indexOfBest);
dBest = norm(vecXBest-vecXSecret)
gBest = norm(vecGBest)
fBest = fBest
%
msg( __FILE__, __LINE__, "In orth() subspace from zero..." );
matV = orth( matX(:,2:end) );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
matD = matX;
[ fFit, vecGFit, matHFit ] = hessfit( matV'*matD, rvecF, matV'*matG, hfPrm );
vecDScale = max( abs(matV'*matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + sxPrm.epsB * max(vecDScale) ) );
datOrthZero = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, matB, cPrm );
%
msg( __FILE__, __LINE__, "In orth() subspace..." );
matV = orth( matX(:,2:end) - vecXBest );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
matD = matX - vecXBest;
[ fFit, vecGFit, matHFit ] = hessfit( matV'*matD, rvecF, matV'*matG, hfPrm );
vecDScale = max( abs(matV'*matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + sxPrm.epsB * max(vecDScale) ) );
datOrth = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, matB, cPrm );

%
msg( __FILE__, __LINE__, "In utorthdrop() subspace..." );
matV = utorthdrop( matX(:,2:end) - vecXBest, sqrt(eps) );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
matD = matX - vecXBest;
[ fFit, vecGFit, matHFit ] = hessfit( matV'*matD, rvecF, matV'*matG, hfPrm );
vecDScale = max( abs(matV'*matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + sxPrm.epsB * max(vecDScale) ) );
datUtorth = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, matB, cPrm );
eigmin = min(eig(matHFit))
vecZCompare = [ datUtorth.lev.matZ(:,end), myhessmin( max([fFit, fBest]), matV'*vecGBest, matHFit, [], sxPrm.trDCoeff*1e8 ) ]
vecZCompareScl = [ datUtorth.levScl.matZ(:,end), myhessmin( max([fFit, fBest]), matV'*vecGBest, matHFit, matB, sxPrm.trDCoeff*1e8 ) ]
msg( __FILE__, __LINE__, "DO THE vecZCompareScl VALUES DISAGREE? WHY?" );
%
msg( __FILE__, __LINE__, "In standard basis..." );
hfPrm = [];
hfPrm.epsHRegu = 0.0;
matD = matX - vecXBest;
[ fFit, vecGFit, matHFit ] = hessfit( matX - vecXBest, rvecF, matG, hfPrm );
vecDScale = max( abs(matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + sxPrm.epsB * max(vecDScale) ) );
datFullSpace = post_sxsolve1208b__calc( vecXBest, eye(sizeX,sizeX), fBest, vecGBest, matHFit, matB, cPrm );

msg( __FILE__, __LINE__, "Using matHSecret..." );
matD = matX - vecXBest;
vecDScale = max( abs(matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + sxPrm.epsB * max(vecDScale) ) );
datHSecret = post_sxsolve1208b__calc( vecXBest, eye(sizeX,sizeX), fBest, vecGBest, matHSecret, matB, cPrm );

%
numFigs++; figure(numFigs);
plot( ...
  %datOrth.lev.rvecDeltaNorm, datOrth.lev.rvecF, 's-', ...
  datUtorth.lev.rvecDeltaNorm, datUtorth.lev.rvecF, 'x-', ...
  datFullSpace.lev.rvecDeltaNorm, datFullSpace.lev.rvecF, '+-', ...
  datHSecret.lev.rvecDeltaNorm, datHSecret.lev.rvecF, '^-' );
grid on;
title( "FActual Along Lev" );
legend( ...
  %"orth", ...
  "utorthdrop", ...
  "std basis", ...
  "h secret", ...
  "location", "NorthEast" );
%
numFigs++; figure(numFigs);
plot( ...
  %datOrth.levScl.rvecBZNorm, datOrth.levScl.rvecF, 's-', ...
  datUtorth.levScl.rvecBZNorm, datUtorth.levScl.rvecF, 'x-', ...
  datFullSpace.levScl.rvecBZNorm, datFullSpace.levScl.rvecF, '+-', ...
  datHSecret.levScl.rvecBZNorm, datHSecret.levScl.rvecF, '^-' );
grid on;
title( "FActual Along LevScl" );
legend( ...
  %"orth", ...
  "utorthdrop", ...
  "std basis", ...
  "h secret", ...
  "location", "NorthEast" );
%
numFigs++; figure(numFigs);
plot( ...
  %datOrth.lev.rvecDeltaNorm, datOrth.lev.rvecFModel, 's-', ...
  datUtorth.lev.rvecDeltaNorm, datUtorth.lev.rvecFModel, 'x-', ...
  datFullSpace.lev.rvecDeltaNorm, datFullSpace.lev.rvecFModel, '+-', ...
  datHSecret.lev.rvecDeltaNorm, datHSecret.lev.rvecFModel, '^-' );
grid on;
title( "FModel Along Lev" );
legend( ...
  %"orth", ...
  "utorthdrop", ...
  "std basis", ...
  "h secret", ...
  "location", "NorthEast" );
%
numFigs++; figure(numFigs);
plot( ...
  %datOrth.levScl.rvecBZNorm, datOrth.levScl.rvecFModel, 's-', ...
  datUtorth.levScl.rvecBZNorm, datUtorth.levScl.rvecFModel, 'x-', ...
  datFullSpace.levScl.rvecBZNorm, datFullSpace.levScl.rvecFModel, '+-', ...
  datHSecret.levScl.rvecBZNorm, datHSecret.levScl.rvecFModel, '^-' );
grid on;
title( "FModel Along LevScl" );
legend( ...
  %"orth", ...
  "utorthdrop", ...
  "std basis", ...
  "h secret", ...
  "location", "NorthEast" );
