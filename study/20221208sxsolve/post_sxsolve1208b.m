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
msg( __FILE__, __LINE__, "Setting epsB..." );
epsB = eps^0.5
%
msg( __FILE__, __LINE__, "Starting from this point..." );
bwdIndex = 35
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
[ fFit, vecGFit, matHFit ] = hessfit( matV'*matX, rvecF, matV'*matG, hfPrm );
datOrthZero = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, matB, cPrm );
%
msg( __FILE__, __LINE__, "In orth() subspace..." );
matV = orth( matX(:,2:end) - vecXBest );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
matD = matX - vecXBest;
[ fFit, vecGFit, matHFit ] = hessfit( matD, rvecF, matV'*matG, hfPrm );
vecDScale = max( abs(matV'*matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + epsB * max(vecDScale) ) );
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
matB = diag( 1.0 ./ ( vecDScale + epsB * max(vecDScale) ) );
datUtorth = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, matB, cPrm );
%
msg( __FILE__, __LINE__, "In full space..." );
hfPrm = [];
hfPrm.epsHRegu = 0.0;
matD = matX - vecXBest;
[ fFit, vecGFit, matHFit ] = hessfit( matX - vecXBest, rvecF, matG, hfPrm );
vecDScale = max( abs(matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + epsB * max(vecDScale) ) );
datFullSpace = post_sxsolve1208b__calc( vecXBest, eye(sizeX,sizeX), fBest, vecGBest, matHFit, matB, cPrm );

msg( __FILE__, __LINE__, "Using matHSecret..." );
matD = matX - vecXBest;
vecDScale = max( abs(matD), [], 2 );
matB = diag( 1.0 ./ ( vecDScale + epsB * max(vecDScale) ) );
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
  "full space", ...
  "h secret", ...
  "location", "NorthEast" );
%
numFigs++; figure(numFigs);
plot( ...
  %datOrth.levScl.rvecBDeltaNorm, datOrth.levScl.rvecF, 's-', ...
  datUtorth.levScl.rvecBDeltaNorm, datUtorth.levScl.rvecF, 'x-', ...
  datFullSpace.levScl.rvecBDeltaNorm, datFullSpace.levScl.rvecF, '+-', ...
  datHSecret.levScl.rvecBDeltaNorm, datHSecret.levScl.rvecF, '^-' );
grid on;
title( "FActual Along LevScl" );
legend( ...
  %"orth", ...
  "utorthdrop", ...
  "full space", ...
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
  "full space", ...
  "h secret", ...
  "location", "NorthEast" );
%
numFigs++; figure(numFigs);
plot( ...
  %datOrth.levScl.rvecBDeltaNorm, datOrth.levScl.rvecFModel, 's-', ...
  datUtorth.levScl.rvecBDeltaNorm, datUtorth.levScl.rvecFModel, 'x-', ...
  datFullSpace.levScl.rvecBDeltaNorm, datFullSpace.levScl.rvecFModel, '+-', ...
  datHSecret.levScl.rvecBDeltaNorm, datHSecret.levScl.rvecFModel, '^-' );
grid on;
title( "FModel Along LevScl" );
legend( ...
  %"orth", ...
  "utorthdrop", ...
  "full space", ...
  "h secret", ...
  "location", "NorthEast" );
