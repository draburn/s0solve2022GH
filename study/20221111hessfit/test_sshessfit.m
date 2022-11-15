clear
mydefs;
setprngstates( 0 );
%
sizeX = 10
numPts = 5
%
secret_f0 = 0.0;
secret_vecX0 = randn(sizeX,1);
secret_matH = randn(sizeX,sizeX);
secret_matH = secret_matH'*secret_matH;
funchF = @(x)( secret_f0 + 0.5*sum((x-secret_vecX0).*( secret_matH * (x-secret_vecX0) ),1) );
funchG = @(x)( secret_matH*(x-secret_vecX0) );
%
matX = randn(sizeX,numPts);
rvecF = funchF(matX);
matG = funchG(matX);
%
[ f_minWas, pt0 ] = min( rvecF )
prm = [];
[ vecX0, matV, fss, vecGss, matHss, datOut ] = sshessfit1115_nonorth( sizeX, numPts, matX, rvecF, matG, pt0, prm );
%
eigVals = eig(matHss)
assert( min(eigVals) > 0.0 );
vecDeltaY = -(matHss\vecGss);
fModelNext = fss + vecDeltaY'*vecGss + 0.5*vecDeltaY'*matHss*vecDeltaY
vecGModelNext = vecGss + matHss*vecDeltaY
vecXNext = vecX0 + matV*vecDeltaY;
fNext = funchF(vecXNext)
vecGNext = funchG(vecXNext)
vecVTGNext = matV'*vecGNext

return;




numUnk = (sizeX*(sizeX+1))/2 + sizeX + 1
numGvn = (1+sizeX)*numPts
%
%
matX(:,1) = 0.0;
secret_matDX = matX - secret_vecX0;
rvecF = secret_f0 + sum( secret_matDX .* ( secret_matH * secret_matDX ), 1 )/2.0;
matG = secret_matH * secret_matDX;
%
%matG += 0.1*randn(size(matG));
%matG = randn(size(matG));
%
msg( __FILE__, __LINE__, "Begin test main..." );
sizeFGL = 1 + sizeX + ((sizeX*(sizeX+1))/2)
%
prm = [];
prm.f0 = rvecF(1);
prm.vecG0 = matG(:,1);
prm.useCnstF = true;
msg( __FILE__, __LINE__, "Calling hessfit()..." );
tic();
[ calc_f, calc_vecG, calc_matH, datOut ] = hessfit( sizeX, numPts, matX, rvecF, matG, prm );
toc();
calc_vecX0 = -(calc_matH\calc_vecG);
secret_vecDX0 = calc_vecX0 - secret_vecX0;
f_minWas = min(rvecF)
f_wouldBe = secret_f0 + (secret_vecDX0'*secret_matH*secret_vecDX0)/2.0
%
msg( __FILE__, __LINE__, "End of test." ); return;
