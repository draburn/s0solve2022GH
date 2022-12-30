clear;
setprngstates(0);
sizeX = 5;
numPts = 6;
%
fCrit = 0.0;
matHFullSpace = mtm(randn(sizeX,sizeX));
vecXCrit = 0.1*randn(sizeX,1);
funchD = @(x)( x - vecXCrit );
funchGOfD = @(d)( matHFullSpace * d );
funchFOfD = @(d)( fCrit + sum( d .* (matHFullSpace*d), 1 )/2.0 );
funchG = @(x)(funchGOfD(funchD(x)));
funchF = @(x)(funchFOfD(funchD(x)));
%
%matX = randn(sizeX,numPts);
matX = [ zeros(sizeX,1), diag([1:numPts-1]) ]
matG = funchG(matX);
rvecF = funchF(matX)
%
indexAnchor = 1;
vecXAnchor = matX(:,indexAnchor);
vecGAnchor = matG(:,indexAnchor);
fAnchor = rvecF(indexAnchor);
matD = matX-vecXAnchor;
[ matV, rvecDrop ] = utorthdrop(matD);
matVTD = matV'*matD;
matVTG = matV'*matG;
[ fFit, vecGammaFit, matHFit, fitDat ] = hessfit( matVTD, rvecF, matVTG );
[ vecZNewt, levsolDat ] = levsol_eig( fFit, vecGammaFit, matHFit );

trCoeff = 0.1
vecCap = max(abs(matVTD),[],2)
matB = diag(1.0./vecCap);
bMax = trCoeff;
[ vecZ, levsolDat ] = levsol_eig( fFit, vecGammaFit, matHFit, matB, bMax );
vecDelta = matV*vecZ
abs(vecDelta)./vecCap
norm(abs(vecDelta)./vecCap)



return;
vecDelta = matV*vecZ;
vecXNew = vecXAnchor + vecDelta;
fNew = funchF(vecXNew)
