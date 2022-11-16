clear
mydefs;
setprngstates( 0 );
%
sizeX = 100
numPts = 5
numUnk = (sizeX*(sizeX+1))/2 + sizeX + 1
numGvn = (1+sizeX)*numPts
%
secret_f0 = 0.0;
secret_vecX0 = randn(sizeX,1);
secret_matH = randn(sizeX,sizeX);
secret_matH = secret_matH'*secret_matH;
funchF = @(x)( secret_f0 + 0.5*sum((x-secret_vecX0).*( secret_matH * (x-secret_vecX0) ),1) );
funchG = @(x)( secret_matH*(x-secret_vecX0) );
%
%matX = randn(sizeX,numPts);
%matX = full(eye(sizeX,numPts));
matX = randn(sizeX,numPts); matX(1,:) += 100.0; matX(2,2:end) += 100.0;
rvecF = funchF(matX);
matG = funchG(matX);
%
matG += randn(size(matG));
%
[ f_minWas, pt0 ] = min( rvecF )
prm = [];
prm.useCnstF = true;
prm.f0 = rvecF(pt0); % Only relevant for hessfit()?
%prm.useCnstG = true;
%prm.vecG0 = matG(:,pt0); % Only relevant for hessfit()?
%[ vecX0, matV, fss, vecGss, matHss, datOut ] = sshessfit1115_nonorth( sizeX, numPts, matX, rvecF, matG, pt0, prm );
%[ vecX0, matV, fss, vecGss, matHss, datOut ] = sshessfit1115_orth( sizeX, numPts, matX, rvecF, matG, pt0, prm );
[ fss, vecGss, matHss, datOut ] = hessfit( sizeX, numPts, matX, rvecF, matG, prm ); vecX0 = zeros(sizeX,1); matV = eye(sizeX);
%
eigVals = eig(matHss);
eigMin = min(eigVals)
%assert( min(eigVals) > 0.0 );
vecDeltaY = -(matHss\vecGss);
fModelNext = fss + vecDeltaY'*vecGss + 0.5*vecDeltaY'*matHss*vecDeltaY
%vecGModelNext = vecGss + matHss*vecDeltaY
vecXNext = vecX0 + matV*vecDeltaY;
fNext = funchF(vecXNext)
%vecGNext = funchG(vecXNext)
%vecVTGNext = matV'*vecGNext

return;
