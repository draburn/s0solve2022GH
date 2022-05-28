function [ vecY, pt_best, iterCount, retCode, datOut ] = findLevPt_0527( vecG, matH, matC, matB, bTrgt, prmIn=[] )
	% Allow C ~= B^T * B. Consequently, ||B*y|| may not be monotonic.
	mydefs;
	datOut = []; % DRaburn 2022-05-27: Not currently supported.
	prm = __init( vecG, matH, matC, matB, bTrgt, prmIn );
	%
	iterCount = 1;
	pt_best = __calcPt( 0.0, vecG, matH, matC, matB, prm );
	if ( isempty(matB) || max(max(abs(matB))) == 0.0 || isempty(bTrgt) )
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Full step was requested." );
		vecY = pt_best.vecY;
		retCode = RETCODE__SUCCESS;
		return;
	endif
	%
	bTol = bTrgt*prm.bRelTol;
	muScale = prm.hScale/prm.cScale;
	pt0 = pt_best;
	havePt1 = false;
	while (~havePt1)
		if ( abs(pt_best.b-bTrgt) <= bTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
			vecY = pt_best.vecY;
			retCode = RETCODE__SUCCESS;
			return;
		elseif ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "IMPOSED STOP: Reached iterMax." );
			vecY = pt_best.vecY;
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		iterCount++;
		%
		mu0 = pt0.mu;
		b0 = pt0.b;
		%mu = muScale * ( (b0/bTrgt) - 1.0 ) ; % This would be correct, but, let's overshoot...
		mu = 10.0*mu0 + muScale*b0/bTrgt;
		pt = __calcPt( mu, vecG, matH, matC, matB, prm );
		%
		if ( abs(pt.b-bTrgt) < abs(pt_best.b-bTrgt) )
			pt_best = pt;
			% We could also check for success here, but, meh.
		endif
		if ( pt.b < bTrgt )
			pt1 = pt;
			havePt1 = true;
		else
			pt0 = pt;
			% We possibly overshot a solution, but we have no way to be sure.
		endif
	endwhile
	%
	applyConstraints = false;
	while (1)
		if ( abs(pt_best.b-bTrgt) <= bTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
			vecY = pt_best.vecY;
			retCode = RETCODE__SUCCESS;
			return;
		elseif ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "IMPOSED STOP: Reached iterMax." );
			vecY = pt_best.vecY;
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		iterCount++;
		%
		mu0 = pt0.mu;
		b0 = pt0.b;
		bPrime0 = pt0.bPrime;
		mu1 = pt1.mu;
		b1 = pt1.b;
		bPrime1 = pt1.bPrime;
		%
		haveSensible0 = false;
		if ( bPrime0 < 0.0 )
			mu_from0 = mu0 + ( (b0/bTrgt) - 1.0 ) * b0 / (-bPrime0);
			if ( mu_from0 < mu1 )
				haveSensible0 = true;
			endif
		endif
		%
		haveSensible1 = false;
		if ( bPrime1 < 0.0 )
			mu_from1 = mu1 + ( (b1/bTrgt) - 1.0 ) * b1 / (-bPrime1);
			if ( mu_from1 > mu0 )
				haveSensible1 = true;
			endif
		endif
		%
		if ( haveSensible1 && haveSensible0 )
			if ( abs(b0-bTrgt) < abs(b1-bTrgt) )
				mu = mu_from0;
			else
				mu = mu_from1;
			endif
		elseif ( haveSensible0 )
			mu = mu_from0;
		elseif ( haveSensible1 )
			mu = mu_from1;
		else
			mu = mu0 + ( (b0/bTrgt) - 1.0 ) * ( mu1 - mu0 ) / ( (b0/b1) - 1.0 );
		endif
		%
		if ( applyConstraints )
			mu = median([ mu0+0.1*(mu1-mu0), mu, mu1-0.1*(mu1-mu0) ]);
		endif
		pt = __calcPt( mu, vecG, matH, matC, matB, prm );
		if ( ~applyConstraints && abs(pt.b-bTrgt) > 0.5 * abs(pt_best.b-bTrgt) )
			applyConstraints = true;
		else
			applyConstraints = false;
		endif
		%
		if ( abs(pt.b-bTrgt) < abs(pt_best.b-bTrgt) )
			pt_best = pt;
		endif
		if ( pt.b < bTrgt )
			pt1 = pt;
		else
			pt0 = pt;
		endif
	endwhile
	return;
endfunction


function prm = __init( vecG, matH, matC, matB, bTrgt, prmIn )
	mydefs;
	%
	%prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__ZERO; % Optimization.
	prm.verbLev = VERBLEV__WARNING; prm.valdLev = VALDLEV__MEDIUM; % Production.
	%prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__HIGH; % Integration.
	%prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__VERY_HIGH; % Dev.
	%prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Debug.
	%prm.bRelTol = 100.0*eps;
	%prm.bRelTol = sqrt(eps);
	prm.bRelTol = 1.0e-2;
	prm.cholRelTol = sqrt(eps);
	prm.epsReguRel = sqrt(eps);
	prm.iterMax = 100;
	prm = overwritefields( prm, prmIn );
	prm.cScale = norm(diag(matC));
	prm.hScale = norm(diag(matH));
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		sz = size( vecG, 1 );
		assert( isrealarray(vecG,[sz,1]) );
		%
		assert( isrealarray(matH,[sz,sz]) );
		assert( issymmetric(matH) );
		assert( min(diag(matH)) >= 0.0 );
		assert( max(diag(matH)) >= max(max(matH))*(1.0-sqrt(eps)) );
		%
		assert( isrealarray(matC,[sz,sz]) );
		assert( issymmetric(matC) );
		assert( min(diag(matC)) > 0.0 );
		assert( max(diag(matC)) >= max(max(matC))*(1.0-sqrt(eps)) );
		%
		if (~isempty(matB))
			szb = size(matB,1);
			assert( isrealarray(matB,[szb,sz]) );
			if ( szb > sz )
				msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, ...
				  "NOTE: Leading dimesion of boundary matrix is greater than problem size." );
				msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, ...
				  "  Replacing the input matB with chol(matB'*matB) might provide a speed-up." );
			endif
		endif
		%
		if (~isempty(bTrgt))
			assert( isscalar(bTrgt) );
			assert( 0.0 < bTrgt );
		endif
		%
		assert( isrealscalar(prm.verbLev) );
		assert( isrealscalar(prm.valdLev) );
		assert( isrealscalar(prm.bRelTol) );
		assert( 0.0 < prm.bRelTol );
		assert( prm.bRelTol < 1.0 );
		assert( isrealscalar(prm.cholRelTol) );
		assert( 0.0 <= prm.cholRelTol );
		assert( prm.cholRelTol <= 1.0 );
		assert( 0.0 < prm.epsReguRel );
		assert( prm.epsReguRel <= 1.0 );
		assert( isrealscalar(prm.iterMax) );
		assert( abs(prm.iterMax-round(prm.iterMax)) < sqrt(eps) );
		assert( 1 <= prm.iterMax );
	endif
	if ( prm.valdLev >= VALDLEV__UNLIMITED )
		eigH = eig(matH);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matH): %0.3e ~ %0.3e", min(eigH), max(eigH) ) );
		eigC = eig(matC);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matC): %0.3e ~ %0.3e", min(eigC), max(eigC) ) );
		%
		if ( min(eigH) < -sqrt(eps)*max(abs(eigH)) )
			error( "Hessian matrix has a clearly negative eigenvalue." );
		elseif ( min(eigC) < -sqrt(eps)*max(abs(eigC)) )
			error( "Curve scaling matrix has a clearly negative eigenvalue." );
		endif
		if ( min(eigC) <= 0.0 )
			msgif ( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, ...
			  "WARNING: Curve scaling matrix appears to be non-positive-definite." );
		endif
	endif
	return;
endfunction

function ptDat = __calcPt( mu, vecG, matH, matC, matB, prm )
	matM = matH + mu*matC;
	[ matR, cholFlag ] = chol( matM);
	if ( 0 == cholFlag && min(diag(matR)) > prm.cholRelTol * max(abs(diag(matR))) )
		[ vecY, vecYPrime, b, bPrime ] = __calcFromChol( matR, vecG, matH, matC, matB );
		ptDat.mu = mu;
		ptDat.vecY = vecY;
		ptDat.vecYPrime = vecYPrime;
		ptDat.b = b;
		ptDat.bPrime = bPrime;
		return;
	endif
	matRegu = (prm.epsReguRel*prm.hScale/prm.cScale);
	matR1 = chol( matM + matRegu );
	matR2 = chol( matM + 2.0*matRegu );
	[ vecY1, vecYPrime1, b1, bPrime1 ] = __calcFromChol( matR1, vecG, matH, matC, matB );
	[ vecY2, vecYPrime2, b2, bPrime2 ] = __calcFromChol( matR2, vecG, matH, matC, matB );
	ptDat.mu = mu;
	ptDat.vecY = 2.0*vecY1 - vecY2;
	ptDat.vecYPrime = 2.0*vecYPrime1 - vecYPrime2;
	ptDat.b = 2.0*b1 - b2;
	ptDat.bPrime = 2.0*bPrime - bPrime2;
	return;
endfunction
function [ vecY, vecYPrime, b, bPrime ] = __calcFromChol( matR, vecG, matH, matC, matB );
	vecY = matR \ ( matR' \ (-vecG) );
	vecYPrime = matR \ ( matR'\(-(matC*vecY)) );
	vecBeta = matB*vecY;
	b = norm( vecBeta );
	if ( 0.0 == b )
		bPrime = 0.0;
	else
		bPrime = (vecBeta'*(matB*vecYPrime)) / b;
	endif
	return;
endfunction

