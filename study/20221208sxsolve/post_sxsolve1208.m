printf("\n\n\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv");
msg( __FILE__, __LINE__, "Hai!" );
matX = [ datOut.vecXBest, datOut.matX ];
matG = [ datOut.vecGBest, datOut.matG ];
rvecF = [ datOut.fBest, datOut.rvecF ];
sizeX = size(matX,1)
numRecsWB = size(matX,2)
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
vecXNext = vecXBest - matV * ( matHFit \ vecGBest );
[ fNext, vecGNext ] = funchFG( vecXNext );
dNext = norm(vecXNext-vecXSecret)
gNext = norm(vecGNext)
fNext = fNext
%
msg( __FILE__, __LINE__, "In orth() subspace..." );
matV = orth( matX(:,2:end) - vecXBest );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
[ fFit, vecGFit, matHFit ] = hessfit( matV'*(matX - vecXBest), rvecF, matV'*matG, hfPrm );
vecXNext = vecXBest - matV * ( matHFit \ vecGBest );
[ fNext, vecGNext ] = funchFG( vecXNext );
dNext = norm(vecXNext-vecXSecret)
gNext = norm(vecGNext)
fNext = fNext
%
msg( __FILE__, __LINE__, "In utorthdrop() subspace..." );
matV = utorthdrop( matX(:,2:end) - vecXBest, sqrt(eps) );
sizeK = size(matV,2)
hfPrm = [];
hfPrm.epsHRegu = 0.0;
[ fFit, vecGFit, matHFit ] = hessfit( matV'*(matX - vecXBest), rvecF, matV'*matG, hfPrm );
vecXNext = vecXBest - matV * ( matHFit \ vecGBest );
[ fNext, vecGNext ] = funchFG( vecXNext );
dNext = norm(vecXNext-vecXSecret)
gNext = norm(vecGNext)
fNext = fNext
%
msg( __FILE__, __LINE__, "In full space..." );
hfPrm = [];
hfPrm.epsHRegu = 0.0;
[ fFit, vecGFit, matHFit ] = hessfit( matX - vecXBest, rvecF, matG, hfPrm );
vecXNext = vecXBest - matHFit \ vecGBest;
[ fNext, vecGNext ] = funchFG( vecXNext );
dNext = norm(vecXNext-vecXSecret)
gNext = norm(vecGNext)
fNext = fNext
%
msg( __FILE__, __LINE__, "Using matHSecret..." );
vecXNext = vecXBest - matHSecret \ vecGBest;
[ fNext, vecGNext ] = funchFG( vecXNext );
dNext = norm(vecXNext-vecXSecret)
gNext = norm(vecGNext)
fNext = fNext
