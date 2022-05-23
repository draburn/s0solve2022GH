	clear;
	numFigs = 0;
	%setprngstates(86372480); sizeX = 500; sizeF = 100; sizeB = sizeX; % Poorly scaled on both ends.
	%setprngstates(92548928); sizeX = 5; sizeF = 5; sizeB = sizeX; % Both hDom and cDom are bad.
	%setprngstates(55143024); sizeX = 10; sizeF = 10; sizeB = sizeX; % Order unity?
	% Switch to not using exp(randn)...
	%setprngstates(46560016); sizeX = 100; sizeF = 100; sizeB = sizeX; % Yeah, okay, the problem is clear now.
	% Switc to using exp(randn) but then normalize...
	%setprngstates(48947952); sizeX = 500; sizeF = 500; sizeB = sizeX; % 14, 9, 6 with old.
	%setprngstates(81975888); sizeX = 500; sizeF = 500; sizeB = sizeX; % 5, 12, 4 with old.
	%setprngstates(25336720); sizeX = 500; sizeF = 500; sizeB = sizeX; % 8, 3, 17 with _basic(); 21, 6, 3 with old.
	%setprngstates(43606048); sizeX = 500; sizeF = 400; sizeB = sizeX; % 21, 4, 3 with old; 13, 5, 3 with old+useLoop0518 orig.
	%setprngstates(51107136); sizeX = 1000; sizeF = 500; sizeB = sizeX; % 19, 4, 2 with old; 12, 6, 4 with old+useLoop0518 orig, 11, 6, 4 with revised.
	%setprngstates(65757664); sizeX = 1000; sizeF = 900; sizeB = sizeX; %0519 wass horrid.
	%setprngstates(64609504); sizeX = 1000; sizeF = 900; sizeB = sizeX; %0519 meh.
	%setprngstates(92034128); sizeX = 1000; sizeF = 900; sizeB = sizeX; %0519 better.
	%setprngstates(12629168); sizeX = 500; sizeF = 200; sizeB = sizeX; %0519 is 9x13x13; old (null) just fails.
	%setprngstates(91450512); sizeX = 1000; sizeF = 900; sizeB = sizeX; %0519 is 7x10x10; old (null) just fails.
	%setprngstates(31746976); sizeX = 1000; sizeF = 900; sizeB = sizeX;
	%setprngstates(98014544); sizeX = 1000; sizeF = 500; sizeB = sizeX;
	setprngstates(65738160); sizeX = 1000; sizeF = 500; sizeB = sizeX; %0521 is 10x5x5
	%
	useBasic = false;
	use0522 = true;
	use0521 = true;
	use0519 = false;
	if (0)
	vecX = randn(sizeX,1) .* exp(randn(sizeX,1)) .* exp(5.0*randn());
	matJ = randn(sizeF,sizeX) .* exp(3.0*randn(sizeF,sizeX)) .* exp(5.0*randn());
	matB = randn(sizeB,sizeX) .* exp(3.0*randn(sizeB,sizeX)) .* exp(5.0*randn());
	elseif (0)
	vecX = randn(sizeX,1);
	matJ = randn(sizeF,sizeX);
	matB = randn(sizeB,sizeX);
	else
	vecX = randn(sizeX,1) .* exp(randn(sizeX,1)) .* exp(5.0*randn());
	matJ = randn(sizeF,sizeX) .* exp(3.0*randn(sizeF,sizeX)) .* exp(5.0*randn());
	matB = randn(sizeB,sizeX) .* exp(3.0*randn(sizeB,sizeX)) .* exp(5.0*randn());
	matJ /= norm(diag(matJ));
	matB /= norm(diag(matB));
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
	alpha_hDom = (3.0/8.0)*sumsq(matB*vecYN)/(vecCYN'*(matH\vecCYN));
	alpha_cDom = 2.0*norm(matB*(matC\vecG))/norm(matB*vecYN);
	%[ norm(matH), alpha_hDom*norm(matC), alpha_cDom*norm(matC) ]
	%alpha = alpha_cDom; matC *= alpha; matB *= sqrt(alpha); bN *= sqrt(alpha);
	%
	msg( __FILE__, __LINE__, "--- Please ignore any warnings above this line! ---" );
	if (1)
	tic();
	prm = [];
	prm.matBTB = matC;
	dat.matC = matC;
	%bTrgt0 = 0.001*bN;
	bTrgt0 = 1.0e-4*bN;
	if ( useBasic )
	[ vecY0, vecYPrime0, b0, bPrime0, n0 ] = findLevPt_basic( vecG, matH, bTrgt0, matB, prm );
	elseif (use0522)
	[ vecY0, datOut ] = findLevPt_0522( vecG, matH, bTrgt0, matB, prm, dat );
	n0 = datOut.iterCount;
	elseif (use0521)
	[ vecY0, datOut ] = findLevPt_0521( vecG, matH, bTrgt0, matB, prm, dat );
	n0 = datOut.iterCount;
	elseif (use0519)
	[ vecY0, datOut ] = findLevPt_0519( vecG, matH, bTrgt0, matB, prm, dat );
	n0 = datOut.iterCount;
	else
	[ vecY0, vecYPrime0, b0, bPrime0, n0 ] = findLevPt( vecG, matH, bTrgt0, matB, prm );
	endif
	%
	prm = [];
	prm.matBTB = matC;
	bTrgt5 = 0.5*bN;
	if ( useBasic )
	[ vecY5, vecYPrime5, b5, bPrime5, n5 ] = findLevPt_basic( vecG, matH, bTrgt5, matB, prm );
	elseif (use0522)
	[ vecY5, datOut ] = findLevPt_0522( vecG, matH, bTrgt5, matB, prm, dat );
	n5 = datOut.iterCount;
	elseif (use0521)
	[ vecY5, datOut ] = findLevPt_0521( vecG, matH, bTrgt5, matB, prm, dat );
	n5 = datOut.iterCount;
	elseif (use0519)
	mydefs;
	%prm.verbLev = VERBLEV__COPIOUS;
	[ vecY5, datOut ] = findLevPt_0519( vecG, matH, bTrgt5, matB, prm, dat );
	n5 = datOut.iterCount;
	%msg( __FILE__, __LINE__, "Goodbye!" ); return;
	else
	[ vecY5, vecYPrime5, b5, bPrime5, n5 ] = findLevPt( vecG, matH, bTrgt5, matB, prm );
	endif
	%
	prm = [];
	prm.matBTB = matC;
	bTrgt9 = 0.9*bN;
	if ( useBasic )
	[ vecY9, vecYPrime9, b9, bPrime9, n9 ] = findLevPt_basic( vecG, matH, bTrgt9, matB, prm );
	elseif (use0522)
	[ vecY9, datOut ] = findLevPt_0522( vecG, matH, bTrgt9, matB, prm, dat );
	n9 = datOut.iterCount;
	elseif (use0521)
	[ vecY9, datOut ] = findLevPt_0521( vecG, matH, bTrgt9, matB, prm, dat );
	n9 = datOut.iterCount;
	elseif (use0519)
	[ vecY9, datOut ] = findLevPt_0519( vecG, matH, bTrgt9, matB, prm, dat );
	n9 = datOut.iterCount;
	else
	[ vecY9, vecYPrime9, b9, bPrime9, n9 ] = findLevPt( vecG, matH, bTrgt9, matB, prm );
	endif
	%
	bAct0 = norm(matB*vecY0);
	bAct5 = norm(matB*vecY5);
	bAct9 = norm(matB*vecY9);
	[ bAct0, bTrgt0, bAct0 - bTrgt0 ]
	[ bAct5, bTrgt5, bAct5 - bTrgt5 ]
	[ bAct9, bTrgt9, bAct9 - bTrgt9 ]
	msg( __FILE__, __LINE__, sprintf( "Num iter: %d, %d, %d.", n0, n5, n9 ) );
	msg( __FILE__, __LINE__, sprintf( "Rel res: %0.3e, %0.3e, %0.3e.", bAct0/bTrgt0 - 1.0, bAct5/bTrgt5 - 1.0, bAct9/bTrgt9 - 1.0 ) );
	toc();
	return;
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
