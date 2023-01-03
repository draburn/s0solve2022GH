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
matX = randn(sizeX,numPts);
matG = funchG(matX);
rvecF = funchF(matX);
%
[ fAnchor, indexAnchor ] = min( rvecF );
vecXAnchor = matX(:,indexAnchor);
vecGAnchor = matG(:,indexAnchor);
matDSans = [ matX(:,1:indexAnchor-1), matX(:,indexAnchor+1:end) ] - vecXAnchor;
matGSans = [ matG(:,1:indexAnchor-1), matG(:,indexAnchor+1:end) ];
%
vecGammaAnchor = matDSans'*vecGAnchor;
matGammaSans = matDSans'*matGSans;
%
%matDSTHDS_model = matGammaSans - vecGammaAnchor
matDSTHDS_model = symm( matGammaSans - vecGammaAnchor )
matDSTHDS_true = matDSans'*matHFullSpace*matDSans
matRes = matDSTHDS_model - matDSTHDS_true
rd = reldiff( matDSTHDS_model, matDSTHDS_true )
