function [ vecDelta, datOut ] = levsol_eig( f0, vecG, matH, matB=[], bMax=[], prm=[] )
	datOut = [];
	%
	% Initialize stuff.
	[ bMax, bTol, fMin, fTol, matRB, prm, initDat ] = __init( f0, vecG, matH, matB, bMax, prm );
	datOut.initDat = initDat;
	datOut.prm = prm;
	%
	% Apply scaling.
	vecG_scaled = matRB \ ( matRB' \ vecG );
	matH_scaled = ( (matRB\(matRB'\matH)) / matRB ) / (matRB');
	%
	% Do main work.
	[ vecDelta_scaled, solveDat ] = __solve( f0, vecG_scaled, matH_scaled, bMax, bTol, fMin, fTol, prm );
	datOut.solveDat = solveDat;
	%
	% Properly de-scale result.
	vecDelta = matRB \ ( matRB' \ vecDelta_scaled );
return;
endfunction


function [ bMax, bTol, fMin, fTol, matRB, prm, initDat ] = __init( f0, vecG, matH, matB, bMax, prmIn )
	initDat = [];
	%
	% Validate input arguments.
	sizeX = size(vecG,1);
	assert( 1 <= sizeX );
	assert( isrealscalar(f0) );
	assert( isrealarray(vecG,[sizeX,1]) );
	assert( isrealarray(matH,[sizeX,sizeX]) );
	assert( issymmetric(matH) );
	if ( ~isempty(matB) )
		assert( isrealarray(matB) );
		assert( issymmetric(matB) );
		if ( ~isdiag(matB) )
			msg( __FILE__, __LINE__, "WARNING: Boundary matrix is not diagonal. Performance may be very slow." );
		endif
	endif
	if ( ~isempty(bMax) )
		assert( isrealscalar(bMax) );
		assert( 0.0 <= bMax );
	endif
	%
	% Set default parameters.
	prm.bTol = myternary( isempty(bMax), [], 0.1*bMax );
	prm.fMin = 0.0;
	prm.fTol = myternary( isempty(mygetfield(prm,"fMin",[])), [], 0.1*f0 );
	prm.matRB = []; % This is just a "suggestion" for matRB.
	prm.epsEig = sqrt(eps);
	%
	% Apply input parameters.
	prm = overwritefields( prm, prmIn );
	%
	% Extract params that get converted to args and process them.
	bTol = prm.bTol;
	fMin = prm.fMin;
	fTol = prm.fTol;
	if ( ~isempty(bMax) )
		assert( isrealscalar(bTol) );
		assert( bTol >= 0.0 );
	else
		assert( isempty(bTol) );
	endif
	if ( ~isempty(fMin) )
		assert( isrealscalar(fMin) );
		assert( fMin <= f0 );
		assert( isrealscalar(fTol) );
		assert( fTol >= 0.0 );
	else
		assert( isempty(fTol) );
	endif
	if ( isempty(matB) )
		assert( isempty(prm.matRB) );
		matRB = eye(sizeX,sizeX);
	elseif ( isempty(prm.matRB) )
		matRB = chol(matB);
	else
		matRB = prm.matRB;
		assert( isrealarray(matRB,[sizeX,sizeX]) );
		assert( istriu(matRB) );
		rd = reldiff( matB, matRB'*matRB );
		if ( rd > sqrt(eps) )
			msg( __FILE__, __LINE__, "WARNING: Provided Cholesky factorization of barrier matrix may be incorrect." );
		endif
	endif
	%
	% Validate other parameters.
	assert( isrealscalar(prm.epsEig) );
	assert( prm.epsEig >= 0.0 );
return;
endfunction


function [ vecDelta, datOut ] = __solve( f0, vecG, matH, bMax, bTol, fMin, fTol, prm )
	datOut = [];
	sizeX = size(vecG,1);
	gScl = norm(vecG);
	hScl = sqrt(sum(sum(matH.^2)));
	%
	% Handle degenerate cases.
	if ( ~isempty(bMax) )
	if ( 0.0 == bMax )
		msg( __FILE__, __LINE__, "WARNING: Max step size is zero." );
		vecDelta = zeros(sizeX,1);
		return;
	endif
	endif
	if ( ~isempty(fMin) )
	if ( f0 == fMin )
		msg( __FILE__, __LINE__, "WARNING: Objective function minimum is initial value." );
		% We technically could still take a non-zero step in this case,
		% allowing f0 to go down to fMin - fTol. But, meh.
		vecDelta = zeros(sizeX,1);
		return;
	endif
	endif
	if ( 0.0 == gScl )
		msg( __FILE__, __LINE__, "WARNING: Initial gradient is zero." );
		vecDelta = zeros(sizeX,1);
		return;
	endif
	if ( 0.0 == hScl )
		msg( __FILE__, __LINE__, "WARNING: Hessian matrix is zero." );
		if ( ~isempty(bMax) && ~isempty(fMin) )
			s = bMax / gScl;
			vecDelta = -s*vecG;
			f = f0 + vecDelta'*vecG;
			if ( f <= fMin - fTol )
				s = ( f0 - fMin ) / (gScl^2);
				vecDelta = -s*vecG;
			endif
			return;
		elseif ( ~isempty(fMin) )
			s = ( f0 - fMin ) / (gScl^2);
			vecDelta = -s * vecG;
			return;
		elseif ( ~isempty(bMax) )
			s = bMax / gScl;
			vecDelta = -s * vecG;
			return;
		else
			error( "Hessian is zero, gradient is not zero, and no constraints have been specified. There is no solution." );
		endif
	endif
	%
	% Get eigenfactorization.
	% This approach is simple, though it may be unaccetpable slow.
	% Note that we can probably calculate just the eigenvalues more quickly,
	%  and repeated Cholesky factorization would likely be faster still.
	% But, KISSxPOITROME.
	[ matPsi, matLambda ] = eig( matH );
	vecLambda = diag( matLambda );
	vecPsiTNG = matPsi' * (-vecG);
	%
	% Itendify maximum (Newton) step size.
	if ( min(vecLambda) < prm.epsEig * max(abs(vecLambda)) )
		if ( min(vecLambda) < 0.0 )
			muMin = abs(min(vecLambda)) + prm.epsEig * max(abs(vecLambda));
		else
			muMin = prm.epsEig * max(abs(vecLambda));
		endif
		sMax = hScl / ( muMin + hScl );
	else
		sMax = 1.0;
	endif
	%
	[ vecDelta_trial, b, f ] = __getDeltaBFOfS( sMax, hScl, f0, vecG, matH, matPsi, vecLambda, vecPsiTNG );
	bSatisfied = isempty(bMax) || ( b <= bMax + bTol );
	fSatisfied = isempty(fMin) || ( f >= fMin - fTol );
	% Note, if both are satisfied, if sMax < 1.0, we could consider "extrapolation" to sMax = 1.0.
	% We would only perform this if it looks like vecDelta is well-behaved as s -> 1.0.
	% Actually, with the eigenfactorization, we could perturb just the least-positive eigenvalues...
	% But, it's not clear the accuracy would even be much better.
	%
	% Backtrack, as needed, to ensure both constrants are satisfied.
	% Note that fzero() may be suboptimal, but has the benefit of simplicity.
	% Potential benefits of a custom zero-finder:
	%  1. Use derivative information.
	%  2. Intelligently bisect when main algorithm is stalling.
	%  3. Consider both goal and constraint at the same time.
	%  4. Store additional informaion (such as calculated vecDelta) to avoid need to re-calculate.
	% See "AMUSING_myfzerod.m" for a decent code which hits points 1, 2, and 4.
	% See "AMUSING_myfzero1216.m" for a not-so-great code which hits point 1 and 3 while attempting 2.
	% 2022-12-25: NOTE: Have not verified that setting TolX works correctly.
	if ( ~bSatisfied )
		fzero_options = optimset();
		if ( ~isempty(bTol) )
			fzero_options = optimset( fzero_options, "TolX", bTol );
		endif
		bTrgt = bMax;
		sMax = fzero( @(s) __getBResOfS( s, hScl, bTrgt, matPsi, vecLambda, vecPsiTNG ), [ 0.0, sMax ] );
		[ vecDelta_trial, b, f ] = __getDeltaBFOfS( sMax, hScl, f0, vecG, matH, matPsi, vecLambda, vecPsiTNG );
		bSatisfied = isempty(bMax) || ( b <= bMax + bTol );
		fSatisfied = isempty(fMin) || ( f >= fMin - fTol );
		assert( bSatisfied );
	endif
	if ( ~fSatisfied )
		fzero_options = optimset();
		if ( ~isempty(fTol) )
			fzero_options = optimset( fzero_options, "TolX", fTol );
		endif
		fTrgt = fMin;
		sMax = fzero( @(s) __getFResOfS( s, hScl, fTrgt, f0, vecG, matH, matPsi, vecLambda, vecPsiTNG ), [ 0.0, sMax ] );
		[ vecDelta_trial, b, f ] = __getDeltaBFOfS( sMax, hScl, f0, vecG, matH, matPsi, vecLambda, vecPsiTNG );
		bSatisfied = isempty(bMax) || ( b <= bMax + bTol );
		fSatisfied = isempty(fMin) || ( f >= fMin - fTol );
		assert( bSatisfied );
		assert( fSatisfied );
	endif
	vecDelta = vecDelta_trial;
return;
endfunction


function [ vecDelta ] = __getDeltaOfMu( mu, matPsi, vecLambda, vecPsiTNG )
	vecDelta = matPsi * ( vecPsiTNG ./ ( vecLambda + mu ) );
return;
endfunction
function [ vecDelta ] = __getDeltaOfS( s, hScl, matPsi, vecLambda, vecPsiTNG )
	if ( 0.0 == s )
		vecDelta = zeros(size(vecPsiTNG));
	else
		vecDelta = __getDeltaOfMu( hScl*((1.0/s)-1.0), matPsi, vecLambda, vecPsiTNG );
	endif
return;
endfunction
function [ vecDelta, b, f ] = __getDeltaBFOfS( s, hScl, f0, vecG, matH, matPsi, vecLambda, vecPsiTNG )
	vecDelta = __getDeltaOfS( s, hScl, matPsi, vecLambda, vecPsiTNG );
	b = norm(vecDelta);
	f = f0 + vecDelta'*vecG + (vecDelta'*matH*vecDelta)/2.0;
return;
endfunction
function [ bRes ] = __getBResOfS( s, hScl, bTrgt, matPsi, vecLambda, vecPsiTNG )
	vecDelta = __getDeltaOfS( s, hScl, matPsi, vecLambda, vecPsiTNG );
	bRes = norm(vecDelta) - bTrgt;
return;
endfunction
function [ fRes ] = __getFResOfS( s, hScl, fTrgt, f0, vecG, matH, matPsi, vecLambda, vecPsiTNG )
	vecDelta = __getDeltaOfS( s, hScl, matPsi, vecLambda, vecPsiTNG );
	fRes = f0 + vecDelta'*vecG + (vecDelta'*matH*vecDelta)/2.0 - fTrgt;
return;
endfunction
