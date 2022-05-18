	clear;
	numFigs = 0;
	%setprngstates(86372480); sizeX = 500; sizeF = 100; sizeB = sizeX; % Poorly scaled on both ends.
	%setprngstates(92548928); sizeX = 5; sizeF = 5; sizeB = sizeX; % Both hDom and cDom are bad.
	%setprngstates(55143024); sizeX = 10; sizeF = 10; sizeB = sizeX; % Order unity?
	% Switch to not using exp(randn)...
	setprngstates(46560016); sizeX = 100; sizeF = 100; sizeB = sizeX; % Yeah, okay, the problem is clear now.
	%
	if (0)
	vecX = randn(sizeX,1) .* exp(randn(sizeX,1)) .* exp(5.0*randn());
	matJ = randn(sizeF,sizeX) .* exp(3.0*randn(sizeF,sizeX)) .* exp(5.0*randn());
	matB = randn(sizeB,sizeX) .* exp(3.0*randn(sizeB,sizeX)) .* exp(5.0*randn());
	else
	vecX = randn(sizeX,1);
	matJ = randn(sizeF,sizeX);
	matB = randn(sizeB,sizeX);
	endif
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
	vecCYN = matC*vecYN;
	alpha_hDom = (3.0/8.0)*sumsq(matB*vecYN)/(vecCYN'*(matH\vecCYN))
	alpha_cDom = 2.0*norm(matB*(matC\vecG))/norm(matB*vecYN)
	[ norm(matH), alpha_hDom*norm(matC), alpha_cDom*norm(matC) ]
	%
	if (0)
	prm = [];
	prm.matBTB = matC;
	bTrgt0 = 0.001*bN
	[ vecY0, vecYPrime0, b0, bPrime0, n0 ] = findLevPt( vecG, matH, bTrgt0, matB, prm );
	%
	prm = [];
	prm.matBTB = matC;
	bTrgt5 = 0.5*bN
	[ vecY5, vecYPrime5, b5, bPrime5, n5 ] = findLevPt( vecG, matH, bTrgt5, matB, prm );
	%
	prm = [];
	prm.matBTB = matC;
	bTrgt9 = 0.9*bN
	[ vecY9, vecYPrime9, b9, bPrime9, n9 ] = findLevPt( vecG, matH, bTrgt9, matB, prm );
	%
	bAct0 = norm(matB*vecY0);
	bAct5 = norm(matB*vecY5);
	bAct9 = norm(matB*vecY9);
	[ bAct0, bTrgt0, bAct0 - bTrgt0 ]
	[ bAct5, bTrgt5, bAct5 - bTrgt5 ]
	[ bAct9, bTrgt9, bAct9 - bTrgt9 ]
	msg( __FILE__, __LINE__, sprintf( "Num iter: %d, %d, %d.", n0, n5, n9 ) );
	msg( __FILE__, __LINE__, sprintf( "Rel res: %0.3e, %0.3e, %0.3e.", bAct0/bTrgt0 - 1.0, bAct5/bTrgt5 - 1.0, bAct9/bTrgt9 - 1.0 ) );
	endif
	%
	%
	numVals = 101;
	foo = linspace( 1.0, 0.0, numVals );
	tVals = ( 1.0 - (foo.^4) ).^4;
	%
	matRegu = sqrt(eps)*matC;
	vecCIG = matC\vecG;
	vecCIHCIG = matC\(matH*vecCIG);
	for n=1:numVals
		t = tVals(n);
		%
		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, t, matC, matRegu );
		vecYVals(:,n) = vecY;
		vecYPrimeVals(:,n) = vecYPrime;
		sVals(n) = s;
		sPrimeVals(n) = sPrime;
		%
		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, t, matC * alpha_cDom, matRegu );
		sVals_cDom(n) = s / sqrt(alpha_cDom);
		sPrimeVals_cDom(n) = sPrime / sqrt(alpha_cDom);
		%
		[ vecY, vecYPrime, s, sPrime ] = calcLevPt( vecG, matH, t, matC * alpha_hDom, matRegu );
		sVals_hDom(n) = s / sqrt(alpha_hDom);
		sPrimeVals_hDom(n) = sPrime / sqrt(alpha_hDom);
		%
		%sVals_modelLo(n) = norm(matB*( t*vecCIG - (t^2)*vecCIHCIG ));
		sVals_modelLo(n) = t*norm(matB*vecCIG);
		sVals_check(n) = t*norm(matB*( (t*matH+(1.0-t)*matC) \ vecG ) );
	endfor
	%[ sVals_cDom(end)/sVals(end), sVals_hDom(end)/sVals(end) ]
	numFigs++; figure(numFigs);
	plot( ...
	  tVals, cap(sVals_modelLo,0.0,sVals(end)), 'o-', 'markersize', 30, 'linewidth', 2, ...
	  tVals, sVals, 'o-', 'markersize', 25, 'linewidth', 2, ...
	  tVals, sVals_cDom, '^-', 'markersize', 20, 'linewidth', 2, ...
	  tVals, sVals_hDom, 'v-', 'markersize', 15, 'linewidth', 2 );
	grid on;
	numFigs++; figure(numFigs);
	semilogy( ...
	  tVals(2:end), sVals(2:end), 'o-', 'markersize', 25, 'linewidth', 2, ...
	  tVals(2:end), sVals_cDom(2:end), '^-', 'markersize', 20, 'linewidth', 2, ...
	  tVals(2:end), sVals_hDom(2:end), 'v-', 'markersize', 15, 'linewidth', 2 );
	grid on;
	numFigs++; figure(numFigs);
	semilogy( ...
	  tVals, sPrimeVals, 'o-', 'markersize', 25, 'linewidth', 2, ...
	  tVals, sPrimeVals_cDom, '^-', 'markersize', 20, 'linewidth', 2, ...
	  tVals, sPrimeVals_hDom, 'v-', 'markersize', 15, 'linewidth', 2 );
	grid on;
