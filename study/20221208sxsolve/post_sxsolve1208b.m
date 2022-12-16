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
datOrthZero = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, [], cPrm );
%
msg( __FILE__, __LINE__, "In orth() subspace..." );
matV = orth( matX(:,2:end) - vecXBest );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
[ fFit, vecGFit, matHFit ] = hessfit( matV'*(matX - vecXBest), rvecF, matV'*matG, hfPrm );
datOrth = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, [], cPrm );

%
msg( __FILE__, __LINE__, "In utorthdrop() subspace..." );
matV = utorthdrop( matX(:,2:end) - vecXBest, sqrt(eps) );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
[ fFit, vecGFit, matHFit ] = hessfit( matV'*(matX - vecXBest), rvecF, matV'*matG, hfPrm );
datUtorth = post_sxsolve1208b__calc( vecXBest, matV, fBest, matV'*vecGBest, matHFit, [], cPrm );
%
msg( __FILE__, __LINE__, "In full space..." );
hfPrm = [];
hfPrm.epsHRegu = 0.0;
[ fFit, vecGFit, matHFit ] = hessfit( matX - vecXBest, rvecF, matG, hfPrm );
datFullSpace = post_sxsolve1208b__calc( vecXBest, eye(sizeX,sizeX), fBest, vecGBest, matHFit, [], cPrm );

msg( __FILE__, __LINE__, "Using matHSecret..." );
datHSecret = post_sxsolve1208b__calc( vecXBest, eye(sizeX,sizeX), fBest, vecGBest, matHSecret, [], cPrm );

%
numFigs++; figure(numFigs);
plot( ...
  datOrth.lev.rvecDeltaNorm, datOrth.lev.rvecF, 's-', ...
  datUtorth.lev.rvecDeltaNorm, datUtorth.lev.rvecF, 'x-', ...
  datFullSpace.lev.rvecDeltaNorm, datFullSpace.lev.rvecF, '+-', ...
  datHSecret.lev.rvecDeltaNorm, datHSecret.lev.rvecF, '^-' );
grid on;
legend( ...
  "orth", ...
  "utorthdrop", ...
  "full space", ...
  "h secret", ...
  "location", "NorthEast" );
  
%
numFigs++; figure(numFigs);
plot( ...
  datOrth.lev.rvecDeltaNorm, datOrth.lev.rvecFModel, 's-', ...
  datUtorth.lev.rvecDeltaNorm, datUtorth.lev.rvecFModel, 'x-', ...
  datFullSpace.lev.rvecDeltaNorm, datFullSpace.lev.rvecFModel, '+-', ...
  datHSecret.lev.rvecDeltaNorm, datHSecret.lev.rvecFModel, '^-' );
grid on;
legend( ...
  "orth", ...
  "utorthdrop", ...
  "full space", ...
  "h secret", ...
  "location", "NorthEast" );
