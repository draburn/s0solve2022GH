function [ vecY, datOut ] = findLevPt_0523( vecG, matH, bTrgt=[], matB=[], prmIn=[], datIn=[] )
	msg( __FILE__, __LINE__, "This is a proof of principle hack of findLevPt_0522." )
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
	if ( isempty(bTrgt) || levDat_newt.b <= bTrgt + bTrgt*prm.bRelTol )
		%[ bTrgt, norm(matB*vecY), levDat_newt.b ]
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
	%prm.verbLev = VERBLEV__UNLIMITED; prm.valdLev = VALDLEV__UNLIMITED; % Dev.
	%prm.bRelTol = 100.0*eps;
	prm.bRelTol = sqrt(eps);
	%prm.bRelTol = 1.0e-4;
	prm.cholRelTol = sqrt(eps);
	prm.epsReguRel = sqrt(eps);
	prm.iterMax = 10;
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
	levDat.funchYPrime = @()( matR \ vecRho ); % Is this right?
endfunction


function [ levDat_best, retCode, iterCount ] = __findLev( levDat0, vecG, matH, bTrgt, matB, prm, dat )
	mydefs;
	if ( prm.valdLev >= VALDLEV__LOW )
		assert( 0.0 < bTrgt );
		assert( bTrgt < levDat0.b );
		assert( levDat0.bPrime < 0.0 );
	endif
	levDat_best = levDat0;
	%
	%
	haveFinite1 = false;
	iterCount = 0;
	
	
	USE_3SEG = true;
	if (USE_3SEG)
	while (1)
		if ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, sprintf( "IMPOSED STOP: Reached iterMax (%d).", prm.iterMax ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		iterCount++;
		
		% Rem: mu0 < mu1, b0 > b1.
		if (~haveFinite1)
			mu1 = +Inf;
			b1 = 0.0;
			vecY1 = zeros(size(vecG));
			vecV1 = -(dat.matC\vecG);
		else
			mu1 = levDat1.mu;
			b1 = levDat1.b;
			vecY1 = levDat1.vecY;
			vecV1 = levDat1.funchYPrime();
		endif
		vecV1 /= norm(vecV1);
		% Consider ray: y = y1 + p1 * v1 for p1 >= 0.
		% Cap based on...
		%  p1 >= 0;
		%  ||B*y|| <= bTrgt;
		%  d/dp1 (omega) <= 0.0.
		%
		assert( b1 < bTrgt );
		% For cap B...
		%  || B*y1 + p1*B*v1 ||^2 = bTgt^2
		%  b1^2 + 2*p1*y1^T*C*v1 + p1^2*||B*v1||^2 = bTrgt^2
		c0 = b1^2 - bTrgt^2;
		c1 = 2*(vecY1' * dat.matC * vecV1 );
		c2 = sumsq(matB*vecV1);
		if ( 0.0 == c2 )
			msg( __FILE__, __LINE__, "EXCEPTION: Constraint matrix is singualr." );
			retCode = RETCODE__BAD_INPUT;
			return;
		endif
		discrim = (c1^2) - (4.0*c0*c2);
		assert( discrim >= 0.0 );
		p1_capB = ( -c1 + sqrt(discrim) ) / (2.0*c2);
		%
		% Conider deltaOmega = g^T*y + 0.5*(y^T*H*y),
		% Critical point is at: (g+H*y)^T * v1 = 0
		% So... (g+H*y1+p1*H*v)^T * v1 = 0
		%  p1*v1^T*H*v1 =-(g+H*y1)^T*v1
		v1thv1 = vecV1' * matH * vecV1;
		if ( 0.0 >= v1thv1 )
			msg( __FILE__, __LINE__, "EXCEPTION: Hessian is singualr." );
			retCode = RETCODE__BAD_INPUT;
			return;
		endif
		p1_capH = -(vecG+(matH*vecY1))'*vecV1/v1thv1;
		%
		%%%p1Min = 0.1*min([ p1_capB, p1_capH ]);
		p1Min = 0.0;
		p1Max = max([ p1Min, min([ p1_capB, p1_capH ]) ]);
		%
		%
		%
		mu0 = levDat0.mu;
		b0 = levDat0.b;
		vecY0 = levDat0.vecY;
		vecV0 = levDat0.funchYPrime();
		vecV0 /= norm(vecV0);
		assert( b0 > bTrgt );
		% Find p0 which minimizes distance between y0 + p0*v0 and y1+p1*v1,
		%  for p1 between p1Min and p1Max, and p0 >= 0.0.
		% Let's solver for the unconstrained case, then cap sensibly --
		%  getting the truly optimal solution might(?) require explict checking
		%  of all possible combinationsof the various bound, but, meh,
		%  even if that were the case, just finding the unconstrained solution
		%  and constraining it is probably good enough for us:
		%  this is jus a model anyway.
		% So...
		if (0)
			% Minimize || y0 - y1 + p0*v0 - p1*v1 ||^2 w.r.t p0 and p1...
			%  d/dp0:  v0'*( y0-y1+p0*v0-p1*v1 ) = 0;
			%  d/dp1: -v1'*( y0-y1+p0*v0-p1*v1 ) = 0;
			% Ergo...
			%   ( v0'*v0) * p0 + (-v1'*v0) * p1 =  v0'*(y1-y0);
			%   (-v1'*v0) * p0 + ( v1'*v1) * p1 = -v1'*(y1-y0);
			%[ sumsq(vecV0), -(vecV1'*vecV0), sumsq(vecV1) ]
			p01 = [ sumsq(vecV0), -(vecV1'*vecV0); -(vecV1'*vecV0), sumsq(vecV1) ] \ [ vecV0'*(vecY1-vecY0); -(vecV1'*(vecY1-vecY0)) ];
			p0 = p01(1);
			p1 = p01(2);
		else
			% Minimize || B* ( y0 - y1 + p0*v0 - p1*v1 ) ||^2 w.r.t p0 and p1...
			%  d/dp0:  v0'*C*( y0-y1+p0*v0-p1*v1 ) = 0;
			%  d/dp1: -v1'*C*( y0-y1+p0*v0-p1*v1 ) = 0;
			% Ergo...
			%   ( v0'*C*v0) * p0 + (-v1'*C*v0) * p1 =  v0'*C*(y1-y0);
			%   (-v1'*C*v0) * p0 + ( v1'*C*v1) * p1 = -v1'*C*(y1-y0);
			%[ sumsq(matB*vecV0), -(vecV1'*dat.matC*vecV0), sumsq(matB*vecV1) ]
			p01 = [ sumsq(matB*vecV0), -(vecV1'*dat.matC*vecV0); -(vecV1'*vecV0), sumsq(vecV1'*dat.matC*vecV0) ] ...
			  \ [ vecV0'*dat.matC*(vecY1-vecY0); -(vecV1'*dat.matC*(vecY1-vecY0)) ];
			p0 = p01(1);
			p1 = p01(2);
		endif
		p1 = median([ p1Min, p1, p1Max ]);
		%%%p0Min = 0.1*p1Min;
		p0Min = 0.0;
		p0 = max([ p0, p0Min ]);
		%
		% We know have the four points of our 3-segment model:
		%  vecY1, vecYA, vecYB, vecY0.
		vecYA = vecY1 + p1*vecV1;
		vecYB = vecY0 + p0*vecV0;
		%
		% Start with the first segment, find where we cross the bounday...
		assert( norm(matB*vecY1) < bTrgt );
		if ( norm(matB*vecYA) >= bTrgt )
			% Solve || matB*(vecY1 + p1 * vecV1) || = bTrgt.
			% We already did this.
			assert( p1_capB <= p1 );
			assert( p1 >= 0.0 );
			vecY3Seg = vecY1 + p1_capB*vecV1;
		elseif ( norm(matB*vecYB) >= bTrgt )
			assert( norm(matB*vecYA) <= bTrgt );
			assert( norm(matB*vecYB) >= bTrgt );
			% Solve || matB*( vecYA + q * ( vecYB - vecYA ) ) || = bTrgt.
			%  ||matB*vecYA||^2 + 2.0*q*(vecYB-vecYA)'*matC*vecYA + q^2*||matB*(vecYB-vecYA)||^2 = bTrgt^2
			c0 = sumsq(matB*vecYA) - bTrgt^2;
			c1 = 2.0*(vecYB-vecYA)'*dat.matC*vecYA;
			c2 = sumsq(matB*(vecYB-vecYA));
			discrim = (c1^2) - (4.0*c0*c2);
			assert( discrim >= 0.0 );
			assert( c2 >= 0.0 );
			q = ( (-c1) + sqrt(discrim) ) / (2.0*c2);
			%qM = ( (-c1) - sqrt(discrim) ) / (2.0*c2)
			assert( q >= 0.0 );
			assert( q <= 1.0 );
			vecY3Seg = vecYA + q*(vecYB-vecYA);
		else
			assert( norm(matB*vecY0) > bTrgt );
			% Solve || matB*( vecYB + q * ( vecY0 - vecYB ) ) || = bTrgt.
			%  ||matB*vecYB||^2 + 2.0*q*(vecY0-vecYB)'*matC*vecYB + q^2*||matB*(vecY0-vecYB)||^2 = bTrgt^2
			c0 = sumsq(matB*vecYB) - bTrgt^2;
			c1 = 2.0*(vecY0-vecYB)'*dat.matC*vecYB;
			c2 = sumsq(matB*(vecY0-vecYB));
			discrim = (c1^2) - (4.0*c0*c2);
			assert( discrim >= 0.0 );
			assert( c2 >= 0.0 );
			q = ( (-c1) + sqrt(discrim) ) / (2.0*c2);
			%qM = ( (-c1) - sqrt(discrim) ) / (2.0*c2)
			vecY3Seg = vecYB + q*(vecY0-vecYB);
		endif
		%
		% Now find pt on Lev curve that closest(ish) to vecY3Seg...
		%  Idealy, we should ahve y3Seg = -(H+mu*C)\g for some mu,
		%   or, equivalently: mu*C*y3Seg = -(g+H*y3Seg) for some mu.
		%  But, since it's not exact, we could min || y3Seg + (H+mu*C)\g ||,
		%   which, if zero, would be the same as min'zing || (H+mu*C)*y3Seg + g ||,
		%  With scaling, this would be || B * ( y3Seg + (H+mu*C)\g ) ||,
		%   but, this is no easier than our original problem!
		%  So, instead, we'll min || B * M \ ( (H+mu*C)*y3Seg + g ) ||
		%   for some sensible M.
		% We'll use M = M0 for now...
		%  || B*M\(H*y3Seg+g) + mu*B\M*(C*y3Seg) ||^2...
		matM0 = matH + mu0*dat.matC;
		vecW1 = matB*( matM0\(matH*vecY3Seg+vecG) );
		vecW2 = matB*( matM0\(dat.matC*vecY3Seg) );
		%  d/dmu || w1 + mu*w2 ||^2 = 0...
		mu = (-vecW1'*vecW2)/(vecW2'*vecW2);
		%[ norm(vecY3Seg), norm((matH+mu*dat.matC)\vecG) ]
		%[ norm(matB*vecY3Seg), norm(matB*((matH+mu*dat.matC)\vecG)) ]
		%[ norm(vecW1), mu*norm(vecW2) ]
		norm(matB*((matH+mu*dat.matC)\vecG));
		stepTypeStr = "3";
		levDat = __calcLev( mu, vecG, matH, matB, prm, dat );
		msgif( prm.verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
		  " %3d:  mu: %9.3e ~ %9.3e (%9.3e);  b: %9.3e ~ %10.3e;  trial: %2s, %9.3e, %10.3e.", ...
		  iterCount, mu0, mu1, mu1-mu0, b0-bTrgt, b1-bTrgt, stepTypeStr, mu, levDat.b-bTrgt ) );
		%
		if ( levDat.b >= b0 || levDat.bPrime >= 0.0 )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Function became non-monotonic." );
			retCode = RETCODE__NUMERICAL_ISSUE;
			return;
		endif
		%
		if ( abs(levDat.b-bTrgt) < abs(levDat_best.b-bTrgt) )
			levDat_best = levDat;
			if ( abs(levDat_best.b-bTrgt) < prm.bRelTol*bTrgt )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
				retCode = RETCODE__SUCCESS;
				return;
			endif
		endif
		%
		if ( levDat.b < bTrgt )
			levDat1 = levDat;
			haveFinite1 = true;
		else
			levDat0 = levDat;
		endif
	endwhile
	endif
	
	
	while (~haveFinite1)
		if ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, sprintf( "IMPOSED STOP: Reached iterMax (%d).", prm.iterMax ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		iterCount++;
		%
		mu0 = levDat0.mu;
		b0 = levDat0.b;
		bPrime0 = levDat0.bPrime;
		mu1 = +Inf;
		b1 = 0.0;
		bPrime1 = 0.0;
		%
		mu = mu0 + 10.0*( (b0/bTrgt) - 1.0 ) * b0 / (-bPrime0); % Intentional overshoot.
		stepTypeStr = "O";
		% Draburn 2022-05-22: We could also consider a step like...
		%bGrad = norm(matB*(dat.matC\vecG));
		%mu = mu0 + bGrad*( (1.0/bTrgt) - (1.0/b0) );
		% but, POITROME.
		%
		if ( mu <= mu0 )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE:Guess went out of bounds." );
			retCode = RETCODE__NUMERICAL_ISSUE;
			return;
		endif
		levDat = __calcLev( mu, vecG, matH, matB, prm, dat );
		msgif( prm.verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
		  " %3d:  mu: %9.3e ~ %9.3e (%9.3e);  b: %9.3e ~ %10.3e;  trial: %2s, %9.3e, %10.3e.", ...
		  iterCount, mu0, mu1, mu1-mu0, b0-bTrgt, b1-bTrgt, stepTypeStr, mu, levDat.b-bTrgt ) );
		if ( levDat.b >= b0 || levDat.bPrime >= 0.0 )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Function became non-monotonic." );
			retCode = RETCODE__NUMERICAL_ISSUE;
			return;
		endif
		%
		if ( abs(levDat.b-bTrgt) < abs(levDat_best.b-bTrgt) )
			levDat_best = levDat;
			if ( abs(levDat_best.b-bTrgt) < prm.bRelTol*bTrgt )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
				retCode = RETCODE__SUCCESS;
				return;
			endif
		endif
		%
		if ( levDat.b < bTrgt )
			levDat1 = levDat;
			haveFinite1 = true;
		else
			levDat0 = levDat;
		endif
	endwhile
	%
	applyConstraints = false;
	while (1)
		if ( iterCount >= prm.iterMax )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, sprintf( "IMPOSED STOP: Reached iterMax (%d).", prm.iterMax ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		endif
		iterCount++;
		%
		mu0 = levDat0.mu;
		b0 = levDat0.b;
		bPrime0 = levDat0.bPrime;
		mu1 = levDat1.mu;
		b1 = levDat1.b;
		bPrime1 = levDat1.bPrime;
		mu_from0 = mu0 + ( (b0/bTrgt) - 1.0 ) * b0 / (-bPrime0);
		mu_from1 = mu1 + ( (b1/bTrgt) - 1.0 ) * b1 / (-bPrime1);
		mu_fromX = mu0 + ( (b0/bTrgt) - 1.0 ) * ( mu1 - mu0 ) / ( (b0/b1) - 1.0 );
		%if ( abs(b0-bTrgt)^2*abs(bPrime1)^0.5 < abs(b1-bTrgt)^2*abs(bPrime0)^0.5 && mu_from0 < mu1 )
		if ( abs(b0-bTrgt) < abs(b1-bTrgt) && mu_from0 < mu1 )
			mu = mu_from0;
			stepTypeStr = "0";
		%elseif ( abs(b0-bTrgt)^2*abs(bPrime1)^0.5 > abs(b1-bTrgt)^2*abs(bPrime0)^0.5 && mu_from1 > mu0 )
		elseif ( abs(b0-bTrgt) > abs(b1-bTrgt) && mu_from1 > mu0 )
			mu = mu_from1;
			stepTypeStr = "1";
		else
			mu = mu_fromX;
			stepTypeStr = "X";
		endif
		if ( applyConstraints )
			mu = median([ mu0+0.1*(mu1-mu0), mu, mu1-0.1*(mu1-mu0) ]);
			stepTypeStr = [ stepTypeStr, "c" ];
		endif
		%
		if ( mu <= mu0 || mu >= mu1 )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE:Guess went out of bounds." );
			retCode = RETCODE__NUMERICAL_ISSUE;
			return;
		endif
		levDat = __calcLev( mu, vecG, matH, matB, prm, dat );
		msgif( prm.verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( ...
		  " %3d:  mu: %9.3e ~ %9.3e (%9.3e);  b: %9.3e ~ %10.3e;  trial: %2s, %9.3e, %10.3e.", ...
		  iterCount, mu0, mu1, mu1-mu0, b0-bTrgt, b1-bTrgt, stepTypeStr, mu, levDat.b-bTrgt ) );
		if ( levDat.b >= b0 || levDat.bPrime >= 0.0 || levDat.b <= b1 )
			msgif( prm.verbLev >= VERBLEV__WARNING, __FILE__, __LINE__, "NUMERICAL ISSUE: Function became non-monotonic." );
			retCode = RETCODE__NUMERICAL_ISSUE;
			return;
		endif
		%
		if ( ~applyConstraints && abs(levDat.b-bTrgt) > 0.5 * abs(levDat_best.b-bTrgt) )
			applyConstraints = true;
		else
			applyConstraints = false;
		endif
		%
		if ( abs(levDat.b-bTrgt) < abs(levDat_best.b-bTrgt) )
			levDat_best = levDat;
			if ( abs(levDat_best.b-bTrgt) < prm.bRelTol*bTrgt )
				msgif( prm.verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, sprintf( "SUCCESS: Converged in %d iterations.", iterCount ) );
				retCode = RETCODE__SUCCESS;
				return;
			endif
		endif
		%
		if ( levDat.b < bTrgt )
			levDat1 = levDat;
		else
			levDat0 = levDat;
		endif
	endwhile
endfunction
