commondefs;
vecFM = funchF( vecXM );
vecFMHat = vecFM/norm(vecFM);
matP = eye(sizeF,sizeF) - (vecFMHat*(vecFMHat'));
%
%fdjacoPrm.vecF0 = vecFM;
fdjacoPrm.epsFD = sqrt(eps);
fdjacoPrm.fdOrder = 2;
funchJ_accurate = @(x) fdjaco( funchF, x, fdjacoPrm );
matJM = funchJ_accurate( vecXM );
[ matUM, matSM, matVM ] = svd( matJM );
vecVNull = matVM(:,end);
%
stepD = 2E-2;
%
vecXPrev = vecXM;
vecX0 = vecXM + stepD * vecVNull;
%
%funchDistFromPrev(repmat(vecX0,[1,2]))
%funchFPRC(repmat(vecX0,[1,2]))
%return;
startTime = time;
for n=1:500
	funchDistFromPrev = @(x)( sqrt(sum( (x-repmat(vecXPrev,[1,size(x,2)])).^2, 1 )) );
	funchFPRC = @(x)( (matP * funchF( x )) + (vecFMHat*((funchDistFromPrev(x)/stepD)-1.0)) );
	groot_prm.funchJ = @(x) fdjaco( funchFPRC, x );
	groot_prm.verbLev = VERBLEV__WARN;
	[ groot_vecX, groot_retCode, goot_datOut ] = groot( funchFPRC, vecX0, groot_prm );
	assert( 0 == groot_retCode );
	matXGroot(:,n) = groot_vecX;
	vecDelta = groot_vecX - vecXPrev;
	[time-startTime n norm(vecDelta) norm(groot_vecX - vecXM) groot_vecX(3)]
	%
	vecXPrev = groot_vecX;
	vecX0 = groot_vecX + (stepD*vecDelta/norm(vecDelta));	;
end
semilogy( sum( (funchF(matXGroot)).^2,1),'o-');grid on
