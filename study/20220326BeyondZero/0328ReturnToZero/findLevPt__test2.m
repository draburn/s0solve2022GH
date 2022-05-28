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
	%setprngstates(65738160); sizeX = 1000; sizeF = 500; sizeB = sizeX; %0521 is 10x5x5
	%
	%%%setprngstates(57428304); sizeX = 100; sizeF = 99; sizeB = sizeX;
	%setprngstates(11053568); sizeX = 500; sizeF = 500; sizeB = sizeX;
	%setprngstates(22617792); sizeX = 500; sizeF = 500; sizeB = sizeX;
	%setprngstates(17078064); sizeX = 500; sizeF = 400; sizeB = sizeX;
	%setprngstates(1255056); sizeX = 1000; sizeF = 900; sizeB = sizeX;
	setprngstates(0); sizeX = 100; sizeF = 100; sizeB = sizeX;
	%
	%
	useBasic = false;
	use0527 = true;
	use0522 = true;
	use0521 = false;
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
	bN = norm(matB*vecYN)
	%
	vecCYN = matC*vecYN;
	alpha_hDom = (3.0/8.0)*sumsq(matB*vecYN)/(vecCYN'*(matH\vecCYN));
	alpha_cDom = 2.0*norm(matB*(matC\vecG))/norm(matB*vecYN);
	%[ norm(matH), alpha_hDom*norm(matC), alpha_cDom*norm(matC) ]
	%alpha = alpha_cDom; matC *= alpha; matB *= sqrt(alpha); bN *= sqrt(alpha);
	%
	%
	if (0)
		numPts = 1001;
		%muPts = 100.0*linspace(1000.0,2000.0,numPts)';
		muPts = 1.0e8*(0.0+linspace(0.0,1.0,numPts).^4);
		for n=1:numPts
			matR = chol( matH + muPts(n)*matC );
			bPts(n) = norm(matB*(matR\(matR'\vecG)));
		endfor
	endif
	if (0)
		n0 = 100;
		mu0 = muPts(n0)
		%b0 = bPts(n0)
		%bPrime0 = (bPts(n0+1)-bPts(n0-1))/(muPts(n0+1)-muPts(n0-1));
		%mu_plus  = muPts(n0+1)-muPts(n0);
		%bP_plus  = (bPts(n0+1)-bPts(n0))/(muPts(n0+1)-muPts(n0));
		%mu_minus = muPts(n0)-muPts(n0-1);
		%bP_minus = (bPts(n0)-bPts(n0-1))/(muPts(n0)-muPts(n0-1));
		%bPrimePrime0 = (bP_plus-bP_minus)/(mu_plus-mu_minus)
		%
		matM0 = matH + mu0*matC;
		vecY0 = matM0\(-vecG);
		vecY1 = matM0\(matC*vecY0);
		vecY2 = matM0\(matC*vecY1);
		alpha_crude = norm(vecY1)/norm(vecY0)
		c2 = norm(matB*vecY0);
		c1 = -2.0*norm(matM0*vecY1);
		c0 = norm(matB*vecY1);
		alpha_crit = -0.5*c1/c2
		discrim = (c1^2)-(4.0*c0*c2);
		if ( discrim >= 0.0 )
			alpha_plus = ( (-c1) + sqrt(discrim) ) / (2.0*c2)
			alpha_minus = ( (-c1) - sqrt(discrim) ) / (2.0*c2)
		endif
		%alpha = 3.0e-5
		alpha = alpha_minus
		%
		dPts = muPts-mu0;
		xPts = alpha*(muPts-mu0);
		for n=1:numPts
			x = xPts(n);
			d = dPts(n);
			vecYPts(:,n) = vecY0/(1.0+x) + (alpha*vecY0-vecY1)*d/((1.0+x)^2) + ( (alpha^2)*vecY0 - 2.0*alpha*vecY1 + vecY2 )*(d^2)/((1.0+x)^3);
			bModelPts(n) = norm(matB*vecYPts(:,n));
		endfor
		%
		numFigs++; figure(numFigs);
		funchVizF = @(f)( f );
		funchVizX = @(x)( x );
		plot( funchVizX(muPts), funchVizF(bPts), 'o-', 'linewidth', 2, funchVizX(muPts), funchVizF(bModelPts), 'x-' );
		grid on;
		axis([0,1e8,0,1e-3])
		%
		numFigs++; figure(numFigs);
		funchVizF = @(f)( f );
		funchVizX = @(x)( x );
		plot( funchVizX(muPts), funchVizF(bPts), 'o-', 'linewidth', 2, funchVizX(muPts), funchVizF(bModelPts), 'x-' );
		grid on;
		axis([0,1e5,0,1e0])
		return
	elseif (0)
		n0 = 100;
		mu0 = muPts(n0)
		b0 = bPts(n0)
		bPrime0 = (bPts(n0+1)-bPts(n0-1))/(muPts(n0+1)-muPts(n0-1));
		mu_plus  = muPts(n0+1)-muPts(n0);
		bP_plus  = (bPts(n0+1)-bPts(n0))/(muPts(n0+1)-muPts(n0));
		mu_minus = muPts(n0)-muPts(n0-1);
		bP_minus = (bPts(n0)-bPts(n0-1))/(muPts(n0)-muPts(n0-1));
		bPrimePrime0 = (bP_plus-bP_minus)/(mu_plus-mu_minus)
		%
		denom = b0 * bPrimePrime0 - (bPrime0^2)
		assert( denom > 0)
		P = bPrime0^2 / denom
		Q = b0
		R = -bPrime0 / ( P * Q )
		bModelPts = Q * ( 1.0 + (muPts-mu0) ).^(-P);
		%
		numFigs++; figure(numFigs);
		funchVizF = @(f)( f );
		funchVizX = @(x)( x );
		plot( funchVizX(muPts), funchVizF(bPts), 'o-', 'linewidth', 2, funchVizX(muPts), funchVizF(bModelPts), 'x-' );
		grid on;
		return
	elseif (0)
		n0 = 100;
		mu0 = muPts(n0)
		b0 = bPts(n0)
		bPrime0 = (bPts(n0+1)-bPts(n0-1))/(muPts(n0+1)-muPts(n0-1));
		mu_plus  = muPts(n0+1)-muPts(n0);
		bP_plus  = (bPts(n0+1)-bPts(n0))/(muPts(n0+1)-muPts(n0));
		mu_minus = muPts(n0)-muPts(n0-1);
		bP_minus = (bPts(n0)-bPts(n0-1))/(muPts(n0)-muPts(n0-1));
		bPrimePrime0 = (bP_plus-bP_minus)/(mu_plus-mu_minus)
		P = b0;
		assert( bPrimePrime0 > 0.0 );
		%Q = sqrt( 0.5 * bPrimePrime0 / b0 );
		%R = bPrime0 + sqrt( 0.5 * bPrimePrime0 * b0 );
		c0 = bPrimePrime0;
		c1 = 4.0*bPrime0;
		c2 = 2.0*b0;
		discrim = (c1^2) - (4.0*c0*c2);
		if ( discrim < 0 )
			msg( __FILE__, __LINE__, "Forcing discrim to zero." );
			discrim = 0.0;
		endif
		Q = ( (-c1) + sqrt(discrim) )/(2.0*c2);
		R = bPrime0 + Q*P;
		bModelPts = P./(1.0+Q*(muPts-mu0)) + R*muPts./((1.0+Q*(muPts-mu0)).^2);
		2.0*(Q^2)*P - 4.0*Q*R
		%
		numFigs++; figure(numFigs);
		funchVizF = @(f)( f );
		funchVizX = @(x)( x );
		plot( funchVizX(muPts), funchVizF(bPts), 'o-', 'linewidth', 2, funchVizX(muPts), funchVizF(bModelPts), 'x-' );
		grid on;
		return
	elseif (0)
		n0 = 100;
		mu0 = muPts(n0)+sqrt(eps)
		b0 = bPts(n0);
		bPrime0 = (bPts(n0+1)-bPts(n0-1))/(muPts(n0+1)-muPts(n0-1));
		bModel0Pts = b0 ./ ( 1.0 - (bPrime0/b0)*(muPts-mu0) );
		msk0 = (bModel0Pts<max(bPts))&(bModel0Pts>0.0);
		%
		n1 = 200;
		mu1 = muPts(n1);
		b1 = bPts(n1);
		bPrime1 = (bPts(n1+1)-bPts(n1-1))/(muPts(n1+1)-muPts(n1-1));
		Q = (b0/b1-1.0)/(mu1-mu0);
		bModel01Pts = b0 ./ ( 1.0 + Q*(muPts-mu0) );
		msk01 = (bModel01Pts<max(bPts))&(bModel01Pts>0.0);
		%
		tPts = 1.0./(1.0+muPts);
		vuPts = 1+muPts;
		numFigs++; figure(numFigs);
		%funchViz = @(x)( 1.0./x - 1.0/max(bPts) );
		funchViz = @(f)( f );
		plot( vuPts, funchViz(bPts), 'o-', vuPts(msk0), funchViz(bModel0Pts(msk0)), '-', 'linewidth', 2, vuPts(msk01), funchViz(bModel01Pts(msk01)), '-' );
		grid on;
		return;
	endif
	%
	%
	msg( __FILE__, __LINE__, "--- Please ignore any warnings above this line! ---" );
	if (1)
	tic();
	prm = [];
	prm.matBTB = matC;
	dat.matC = matC;
	bTrgt0 = 0.001*bN;
	if ( useBasic )
	[ vecY0, vecYPrime0, b0, bPrime0, n0 ] = findLevPt_basic( vecG, matH, bTrgt0, matB, prm );
	elseif (use0527)
	[ vecY0, pt0, n0, retCode0, datOut ] = findLevPt_0527( vecG, matH, matC, matB, bTrgt0, prm, dat );
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
	elseif (use0527)
	[ vecY5, pt5, n5, retCode5, datOut ] = findLevPt_0527( vecG, matH, matC, matB, bTrgt5, prm, dat );
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
	elseif (use0527)
	[ vecY9, pt9, n9, retCode9, datOut ] = findLevPt_0527( vecG, matH, matC, matB, bTrgt9, prm, dat );
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
