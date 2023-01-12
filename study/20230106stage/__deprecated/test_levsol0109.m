	clear;
	setprngstates(0);
	sizeX = 5;
	%%%matH = mtm(randn(sizeX,sizeX));
	matH = mtm(randn(sizeX-2,sizeX)) - mtm(randn(2,sizeX));
	%%%matA = randn(sizeX,sizeX);
	%%%matH = (matA'+matA)/2.0
	fCrit = abs(randn())+1E1
	vecXCrit = randn(sizeX,1)
	matB = diag(abs(randn(sizeX,1)));
	%
	funchD = @(x)( x - vecXCrit );
	funchFOfD = @(d)( fCrit + (d'*matH*d)/2.0 );
	funchGOfD = @(d)( matH*d );
	funchF = @(x) funchFOfD(funchD(x));
	funchG = @(x) funchGOfD(funchD(x));
	%
	vecX0 = zeros(sizeX,1);
	f0 = funchF(vecX0)
	vecG = funchG(vecX0);
	%bInf = norm(matB*(vecXCrit-vecX0))
	%
	vecDeltaN = levsol0109( f0, vecG, matH, matB );
	fN = funchF( vecX0 + vecDeltaN )
	vecGN = funchG( vecX0 + vecDeltaN )
	f0OC = f0 - fCrit
	fNOC = fN - fCrit
	%return
	bN = norm(matB*vecDeltaN)
	xN = norm(vecDeltaN);
	if ( min(eig(matH)) > 0.0 && fCrit > 0.0 )
		assert( reldiff( vecDeltaN, vecXCrit -vecX0 ) < sqrt(eps) );
	endif
	%
	vecDeltaNP = levsol0109( f0, vecG, matH, matB, bN+0.01 );
	vecDeltaN0 = levsol0109( f0, vecG, matH, matB, bN );
	vecDeltaNM = levsol0109( f0, vecG, matH, matB, bN-0.01 );
	mat_nearEnd = [ vecDeltaN, vecDeltaNP, vecDeltaN0, vecDeltaNM ]
	%
	vecDelta06 = levsol0109( f0, vecG, matH, matB, bN*0.6 );
	vecDelta05 = levsol0109( f0, vecG, matH, matB, bN*0.5 );
	vecDelta04 = levsol0109( f0, vecG, matH, matB, bN*0.4 );
	mat_nearMid = [ vecDelta06, vecDelta05, vecDelta04 ]
	%
	vecDelta002 = levsol0109( f0, vecG, matH, matB, bN*0.02 );
	vecDelta001 = levsol0109( f0, vecG, matH, matB, bN*0.01 );
	vecDelta000 = levsol0109( f0, vecG, matH, matB, bN*0.00 );
	mat_nearZer = [ vecDelta002, vecDelta001, vecDelta000 ]
	%
	bMax = []
	xMax = []
	vecDelta = levsol0109( f0, vecG, matH, matB, bMax, xMax );
	b = norm(matB*vecDelta)
	x = norm(vecDelta)
	%
	bMax = []
	xMax = 2.0
	vecDelta = levsol0109( f0, vecG, matH, matB, bMax, xMax );
	b = norm(matB*vecDelta)
	x = norm(vecDelta)
	%
	bMax = 2.0
	xMax = []
	vecDelta = levsol0109( f0, vecG, matH, matB, bMax, xMax );
	b = norm(matB*vecDelta)
	x = norm(vecDelta)
	%
	bMax = 2.0
	xMax = 2.0
	vecDelta = levsol0109( f0, vecG, matH, matB, bMax, xMax );
	b = norm(matB*vecDelta)
	x = norm(vecDelta)
	%
	bMax = 2.0
	xMax = 1.5
	vecDelta = levsol0109( f0, vecG, matH, matB, bMax, xMax );
	b = norm(matB*vecDelta)
	x = norm(vecDelta)
	%
	msg( __FILE__, __LINE__, "vvvvv" );
	matH(:,:) = 0.0;
	vecDelta = levsol0109( f0, vecG, matH, matB, bMax, xMax );
	msg( __FILE__, __LINE__, "^^^^^" );
	vecDelta
	%
	msg( __FILE__, __LINE__, "Please check the above results for reasonableness." );
