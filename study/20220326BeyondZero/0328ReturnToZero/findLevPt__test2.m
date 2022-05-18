	clear;
	numFigs = 0;
	setprngstates();
	sizeX = 500;
	sizeF = 100;
	sizeB = sizeX;
	%
	vecX = randn(sizeX,1) .* exp(randn(sizeX,1)) .* exp(5.0*randn());
	matJ = randn(sizeF,sizeX) .* exp(3.0*randn(sizeF,sizeX)) .* exp(5.0*randn());
	matB = randn(sizeB,sizeX) .* exp(3.0*randn(sizeB,sizeX)) .* exp(5.0*randn());
	vecF = matJ*vecX;
	%
	vecG = matJ'*vecF;
	matH = matJ'*matJ;
	matC = matB'*matB;
	matH += sqrt(eps)*diag(diag(matH));
	%
	vecYN = matH\vecG;
	bN = norm(matB*vecYN);
	%
	prm = [];
	prm.matBTB = matC;
	bTrgt0 = 0.001*bN
	[ vecY0, vecYPrime0, b0, bPrime0 ] = findLevPt( vecG, matH, bTrgt0, matB, prm );
	%
	prm = [];
	prm.matBTB = matC;
	bTrgt5 = 0.5*bN
	[ vecY5, vecYPrime5, b5, bPrime5 ] = findLevPt( vecG, matH, bTrgt5, matB, prm );
	%
	prm = [];
	prm.matBTB = matC;
	bTrgt9 = 0.999*bN
	[ vecY9, vecYPrime9, b9, bPrime9 ] = findLevPt( vecG, matH, bTrgt9, matB, prm );
	%
	bAct0 = norm(matB*vecY0);
	bAct5 = norm(matB*vecY5);
	bAct9 = norm(matB*vecY9);
	[ bAct0, bTrgt0, bAct0 - bTrgt0 ]
	[ bAct5, bTrgt5, bAct5 - bTrgt5 ]
	[ bAct9, bTrgt9, bAct9 - bTrgt9 ]
	msg( __FILE__, __LINE__, sprintf( "Rel res: %0.3e, %0.3e, %0.3e.", bAct0/bTrgt0 - 1.0, bAct5/bTrgt5 - 1.0, bAct9/bTrgt9 - 1.0 ) );
	%
	%
	numVals = 101;
	foo = linspace( 1.0, 0.0, numVals );
	tVals = ( 1.0 - (foo.^4) ).^4;
	%
	matRegu = sqrt(eps)*matC;
	for n=1:numVals
		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, tVals(n), matC, matRegu );
		vecYVals(:,n) = vecY;
		vecYPrimeVals(:,n) = vecYPrime;
		sVals(n) = s;
		sPrimeVals(n) = sPrime;
	endfor
	numFigs++; figure(numFigs);
	plot( ...
	  tVals, sVals, 'o-' );
	grid on;
	numFigs++; figure(numFigs);
	plot( ...
	  tVals, sPrimeVals, 'o-' );
	grid on;
