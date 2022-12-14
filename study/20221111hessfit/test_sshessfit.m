clear
mydefs;
setprngstates( 0 );
%
sizeX1 = 1E1;
sizeX = 2E4
numPts = 1E2
numUnk = (sizeX*(sizeX+1))/2 + sizeX + 1
numGvn = (1+sizeX)*numPts
%
secret_f0 = 0.0;
secret_vecX0 = 1000.0*randn(sizeX,1);
msg( __FILE__, __LINE__, "Generating Hessian..." ); tic();
%%%secret_matH = randn(sizeX,sizeX); secret_matH = secret_matH'*secret_matH;
secret_matH = randn(sizeX1,sizeX); secret_matH = secret_matH'*secret_matH;
msgnnl( __FILE__, __LINE__, "Finished generating Hessian. " ); toc();
funchF = @(x)( secret_f0 + 0.5*sum((x-secret_vecX0).*( secret_matH * (x-secret_vecX0) ),1) );
funchG = @(x)( secret_matH*(x-secret_vecX0) );
%
matX = randn(sizeX,numPts);
%matX = full(eye(sizeX,numPts));
%matX = randn(sizeX,numPts); matX(1,:) += 100.0; matX(2,2:end) += 100.0;
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
msg( __FILE__, __LINE__, "Generating fit..." ); tic();
[ vecX0, matV, fss, vecGss, matHss, datOut ] = sshessfit1115_nonorth( sizeX, numPts, matX, rvecF, matG, pt0, prm );
%[ vecX0, matV, fss, vecGss, matHss, datOut ] = sshessfit1115_orth( sizeX, numPts, matX, rvecF, matG, pt0, prm );
%[ fss, vecGss, matHss, datOut ] = hessfit( sizeX, numPts, matX, rvecF, matG, prm ); vecX0 = zeros(sizeX,1); matV = eye(sizeX);
msgnnl( __FILE__, __LINE__, "Finished generating fit. " ); toc();
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
fallFrac = 1.0 - fNext/f_minWas

return;
