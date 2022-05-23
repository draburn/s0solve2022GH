function [ vecY, datOut ] = findLevPt_0522( vecG, matH, bTrgt=[], matB=[], prmIn=[], datIn=[] )
	mydefs;
	[ matB, prm, dat ] = __init( vecG, matH, bTrgt, matB, prmIn, datIn );
	if ( 0.0 == bTrgt )
		vecY = zeros(size(vecG));
		datOut.mu = +Inf;
		datOut.b = 0.0;
		datOut.bPrime = 0.0;
		datOut.vecYPrime = zeros(size(vecG));
		datOut.iterCount = 0;
		msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Zero step was requested." );
		return;
	endif
	%
	levDat_newt = __calcLev_newt( vecG, matH, matB, prm, dat );
	datOut.levDat_newt = levDat_newt;
	if ( isempty(bTrgt) || abs(bTrgt-levDat_newt.b) <= bTrgt*prm.bRelTol )
		vecY = levDat_newt.vecY;
		datOut.levDat = levDat_newt;
		datOut.iterCount = 0;
		if ( isempty(bTrgt) )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Full step was requested." );
		else
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Full step is within tolerance." );
		endif
		return;
	endif
	%
	[ levDat, retCode, iterCount ] = __findLev( levDat_newt, vecG, matH, bTrgt, matB, prm, dat );
	vecY = levDat.vecY;
	datOut.levDat = levDat;
	datOut.retCode = retCode;
	datOut.iterCount = iterCount;
	return;
endfunction


function [ matB, prm, dat ] = __init( vecG, matH, bTrgt=[], matB=[], prmIn=[], datIn=[] )
	mydefs;
	%
	if ( isempty(matB) )
		matB = eye(size(matH));
	endif
	%
	prm.verbLev = VERBLEV__FLAGGED; prm.valdLev = VALDLEV__ZERO; % Production.
	prm.verbLev = VERBLEV__MAIN; prm.valdLev = VALDLEV__MEDIUM; % Integration.
	prm.verbLev = VERBLEV__DETAILS; prm.valdLev = VALDLEV__HIGH; % Performance testing.
	prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Dev.
	prm.bRelTol = sqrt(eps);
	%prm.bRelTol = 1.0e-4;
	prm.cholRelTol = sqrt(eps);
	prm.epsReguRel = sqrt(eps);
	prm.iterMax = 100;
	prm = overwritefields( prm, prmIn );
	%
	matC = mygetfield( datIn, "matC", [] );
	if (isempty(matC))
		matC = matB'*matB;
	endif
	dat.matC = matC;
	dat.hScale = norm(diag(matH));
	if ( 0 == dat.hScale )
		error( "Diagonal of Hessian matrix is all zeros." );
	endif
	dat.cScale = norm(diag(matC));
	if ( 0 == dat.cScale )
		error( "Constraint matrix is singular." );
	endif
	dat.epsReguScaled = prm.epsReguRel * dat.hScale / dat.cScale;
	dat = overwritefields( dat, prmIn );
	%
	if ( prm.valdLev >= VALDLEV__LOW )
		sz = size( vecG, 1 );
		assert( isrealarray(vecG,[sz,1]) );
		assert( isrealarray(matH,[sz,sz]) );
		assert( issymmetric(matH) );
		assert( min(diag(matH)) >= 0.0 );
		assert( max(diag(matH)) >= max(max(matH)) );
		if (~isempty(bTrgt))
			assert( isscalar(bTrgt) );
			assert( 0.0 < bTrgt );
		endif
		szb = size(matB,1);
		assert( isrealarray(matB,[szb,sz]) );
		if ( szb > sz )
			msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, ...
			  "NOTE: Leading dimesion of boundary matrix is greater than problem size." );
			msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, ...
			  "  Replacing the input matB with chol(matB'*matB) might provide a speed-up." );
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
		%
		assert( isrealarray(dat.matC,[sz,sz]) );
		assert( issymmetric(dat.matC) );
		assert( min(diag(dat.matC)) >= 0.0 );
		assert( max(diag(dat.matC)) >= max(max(dat.matC)) );
		assert( isrealscalar(dat.epsReguScaled) );
		assert( 0.0 < dat.epsReguScaled );
	endif
	if ( prm.valdLev >= VALDLEV__MEDIUM );
		assert( reldiff(dat.matC,matB'*matB) < sqrt(eps) );
	endif
	if ( prm.valdLev >= VALDLEV__UNLIMITED )
		eigH = eig(matH);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matH): %0.3e ~ %0.3e", min(eigH), max(eigH) ) );
		eigC = eig(dat.matC);
		msgif( prm.verbLev >= VERBLEV__NOTE, __FILE__, __LINE__, sprintf( "eig(matC): %0.3e ~ %0.3e", min(eigC), max(eigC) ) );
		%
		if ( min(eigH) < -sqrt(eps)*max(abs(eigH)) )
			error( "Hessian matrix has a clearly negative eigenvalue." );
		elseif ( min(eigC) < -sqrt(eps)*max(abs(eigC)) )
			error( "Constraint matrix has a clearly negative eigenvalue." );
		endif
		if ( min(eigC) <= 0.0 )
			msgif ( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, ...
			  "WARNING: Constraint matrix appears to be non-positive-definite." );
		endif
	endif
	return;
endfunction


% If chol() failes, use linear extrapolation.
% Since calculation of dy/dmu would require an additional backsub, return merely a function to allow its calculation.
function levDat = __calcLev_newt( vecG, matH, matB, prm, dat )
	mydefs;
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag )
	if ( min(diag(matR)) > prm.cholRelTol * max(abs(diag(matR))) )
		[ b, bPrime, vecY, vecRho ] = __calcFromChol( vecG, matR, matB, dat.matC );
		levDat.mu = 0.0;
		levDat.b = b;
		levDat.bPrime = bPrime;
		levDat.vecY = vecY;
		levDat.funchYPrime = @()( matR \ vecRho );
		return;
	endif
	endif
	clear matR;
	%
	msgif( prm.verbLev >= VERBLEV__INFO, __FILE__, __LINE__, "Using extrapolation for Newton step." );
	% Note that epsReguScaled should be scaled for matC.
	[ matR1, cholFlag1 ] = chol( matH + dat.epsReguScaled * dat.matC );
	if ( 0 ~= cholFlag1 )
		error( "Cholesky factorization failed even with regularization matrix." );
	endif
	[ matR2, cholFlag2 ] = chol( matH + (2.0*dat.epsReguScaled) * dat.matC );
	if ( 0 ~= cholFlag2 )
		error( "Cholesky factorization (somehow) failed with regularization matrix second time." );
	endif
	%
	[ b1, bPrime1, vecY1, vecRho1 ] = __calcFromChol( vecG, matR1, matB, dat.matC );
	[ b2, bPrime2, vecY2, vecRho2 ] = __calcFromChol( vecG, matR2, matB, dat.matC );
	levDat.mu = 0.0;
	levDat.b = (2.0*b1) - b2;
	levDat.bPrime = (2.0*bPrime1) - bPrime2;
	levDat.vecY = (2.0*vecY1) - vecY2;
	levDat.funchYPrime = @()(  (2.0*( matR1 \ vecRho1 )) - ( matR2 \ vecRho2 )  );
	return;
endfunction


% Math...
%  C = B^T * B
%  M = H + mu*C
%  y = M^-1 * (-g)
%  b = || B *y ||
%  dM/dmu = - M^-1 * C * M^-1
%  dy/dmu = - M^-1 * C * M^-1 * (-g) = - M^-1 * C * y
%  d(b^2)/dmu = -2 * y^T * C * dy/dmu = -2 * || M^-(1/2) * C * y ||^2
%  db/dmu = (d(b^2)/dmu) / ( 2 * b ) = || M^-(1/2) * C * y ||^2 / b
% Use rho = M^(-1/2) * C * y. (I forget why it's called "rho".)
% Use M = R^T * R.
% Note that dy/dmu is not calculated here, but can be calculated from rho and R.
function [ b, bPrime, vecY, vecRho ] = __calcFromChol( vecG, matR, matB, matC )
	vecY = matR \ ( matR' \ (-vecG) );
	b = norm( matB * vecY );
	vecRho = matR' \ ( matC * vecY );
	bPrime = -sumsq( vecRho ) / b;
	return;
endfunction



function levDat = __calcLev( mu, vecG, matH, matB, prm, dat )
	matR = chol( matH + mu*dat.matC );
	vecY = matR \ ( matR' \ (-vecG) );
	b = norm( matB * vecY );
	vecRho = matR' \ ( dat.matC * vecY );	
	levDat.mu = mu;
	levDat.b = b;
	levDat.bPrime = -sumsq( vecRho ) / b;
	levDat.vecY = vecY;
	levDat.funchYPrime = @()( matR \ vecRho );
endfunction


function [ levDat, retCode, iterCount ] = __findLev( levDat0, vecG, matH, bTrgt, matB, prm, dat )
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < bTrgt );
		assert( bTrgt < levDat0.b );
		assert( levDat0.bPrime < 0.0 );
	endif
	levDat_best = levDat0;
	haveFinite1 = false;
	iterCount = 0;
	%
	if ( bTrgt < 0.1*levDat0.b )
	%if ( 1 )
		% Try large-mu form.
		[ matR, cholFlag ] = chol( dat.matC );
		if ( 0~=cholFlag )
			error( "Cholesky factorization of constraint matrix failed." );
		endif
		mu = norm(matB*(matR\(matR'\vecG)))/bTrgt;
		if ( mu > levDat0.mu )
			levDat = __calcLev( mu, vecG, matH, matB, prm, dat );
			if ( abs(levDat.b-bTrgt) < abs(levDat_best.b-bTrgt) )
				levDat_best = levDat;
				if ( abs(levDat_best.b-bTrgt) < prm.bRelTol*bTrgt )
					msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "SUCCESS: Large-mu form was within tolerance." );
					retCode = RETCODE__SUCCESS;
					return;
				endif
			endif
			if ( levDat.b >= levDat0.b || levDat.bPrime >= 0.0 )
				msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Function became non-monotonic." );
				retCode = RETCODE__NUMERICAL_ISSUE;
				return;
			endif
			if ( levDat.b < bTrgt )
				levDat1 = levDat;
				haveFinite1 = true;
			else
				levDat0 = levDat;
			endif
		endif
	endif
	%
	while (~haveFinite1)
		if ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, sprintf( "IMPOSED STOP: Reached iterMax (%d).", prm.iterMax ) );
			levDat.retCode = RETCODE__IMPOSED_STOP;
			break;
		endif
		iterCount++;
		%
		mu0 = levDat0.mu;
		b0 = levDat0.b;
		bPrime0 = levDat0.bPrime;
		mu1 = +Inf;
		b1 = 0.0;
		bPrime1 = 0.0;
		mu = mu0 - ( (b0/bTrgt) - 1.0 ) * b0 / bPrime0;
		stepTypeStr = "0";
		if ( mu <= levDat0.mu )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Guess went out of bounds." );
			retCode = RETCODE__NUMERICAL_ISSUE;
		endif
		levDat = __calcLev( mu, vecG, matH, matB, prm, dat );
		msgif( prm.verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
		  " %3d:  mu: %9.3e ~ %9.3e (%8.3e);  b: %9.3e ~ %10.3e;  trial: %s, %9.3e, %10.3e.", ...
		  iterCount, mu0, mu1, mu1-mu0, b0-bTrgt, b1-bTrgt, stepTypeStr, mu, levDat.b-bTrgt ) );
		if ( abs(levDat.b-bTrgt) < abs(levDat_best.b-bTrgt) )
			levDat_best = levDat;
			if ( abs(levDat_best.b-bTrgt) < prm.bRelTol*bTrgt )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
				retCode = RETCODE__SUCCESS;
				return;
			endif
		endif
		if ( levDat.b >= levDat0.b || levDat.bPrime >= 0.0 )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Function became non-monotonic." );
			retCode = RETCODE__NUMERICAL_ISSUE;
			break;
		endif
		%
		if ( levDat.b < bTrgt )
			haveFinite1 = true;
			levDat1 = levDat;
		else
			levDat0 = levDat;
		endif
	endwhile
	%
	while (1)
		if ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "IMPOSED STOP: Reached iterMax (%d).", prm.iterMax ) );
			datOut.retCode = RETCODE__IMPOSED_STOP;
			break;
		endif
		iterCount++;
		
		%%%
		% This part needs work!
		mu0 = levDat0.mu;
		b0 = levDat0.b;
		bPrime0 = levDat0.bPrime;
		mu1 = levDat1.mu;
		b1 = levDat1.b;
		bPrime1 = levDat1.bPrime;
		mu_from0 = mu0 - ( (b0/bTrgt) - 1.0 ) * b0 / bPrime0;
		mu_from1 = mu1 - ( (b1/bTrgt) - 1.0 ) * b1 / bPrime1;
		mu_fromX = mu0 + ( (b0/bTrgt) - 1.0 ) * ( mu1 - mu0 ) / ( (b0/b1) - 1.0 );
		mu_for0 = mu0 + 0.3*(mu1-mu0);
		mu_for1 = mu0 + 0.7*(mu1-mu0);
		%if ( mu_from0 < mu_for0 && mu_from1 < mu_for0 && mu_fromX < mu_for0 )
		if ( abs(b0-bTrgt)*sqrt(sqrt(mu_from0-mu0)) < abs(b1-bTrgt)*sqrt(sqrt(mu1-mu_from1)) && mu_from0 < 0.9*mu1 )
			mu = mu_from0;
			stepTypeStr = "0";
		%elseif ( mu_from0 > mu_for1 && mu_from1 > mu_for1 && mu_fromX > mu_for1 )
		elseif ( abs(b1-bTrgt)*sqrt(sqrt(mu1-mu_from1)) < abs(b0-bTrgt)*sqrt(sqrt(mu_from0-mu0)) && mu_from1 > 1.1*mu0 )
			mu = mu_from1;
			stepTypeStr = "1";
		else
			mu = mu_fromX;
			stepTypeStr = "X";
		endif
		%%%
		
		if ( mu <= levDat0.mu || mu >= levDat1.mu )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Guess went out of bounds." );
			retCode = RETCODE__NUMERICAL_ISSUE;
		endif
		levDat = __calcLev( mu, vecG, matH, matB, prm, dat );
		msgif( prm.verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
		  " %3d:  mu: %9.3e ~ %9.3e (%8.3e);  b: %9.3e ~ %10.3e;  trial: %s, %9.3e, %10.3e.", ...
		  iterCount, mu0, mu1, mu1-mu0, b0-bTrgt, b1-bTrgt, stepTypeStr, mu, levDat.b-bTrgt ) );
		if ( abs(levDat.b-bTrgt) < abs(levDat_best.b-bTrgt) )
			levDat_best = levDat;
			if ( abs(levDat_best.b-bTrgt) < prm.bRelTol*bTrgt )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
				retCode = RETCODE__SUCCESS;
				return;
			endif
		endif
		if ( levDat.b >= levDat0.b || levDat.bPrime >= 0.0 || levDat.b <= levDat1.b )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Function became non-monotonic." );
			retCode = RETCODE__NUMERICAL_ISSUE;
			break;
		endif
		%
		if ( levDat.b < bTrgt )
			levDat1 = levDat;
		else
			levDat0 = levDat;
		endif
	endwhile
endfunction
