% function vecX = levsol0109( f0, vecG, matH, matB=[], bMax=[], xMax=[], prm=[] )
% Uses "lambdaFloor" for positive-defness and fMin.
% Allows for constrants on both ||B*x|| and ||x||.

function vecX = levsol0109( f0, vecG, matH, matB=[], bMax=[], xMax=[], prm=[] )
	% Validate input.
	sz = size(vecG,1);
	assert( isrealscalar(f0) );
	assert( f0 > 0.0 );
	assert( isrealarray(vecG,[sz,1]) );
	assert( norm(vecG) > 0.0 );
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	assert( sum(sum(abs(matH))) > 0.0 );
	if ( isempty(matB) )
		vecB = ones(sz,1);
	else
		assert( isrealarray(matB,[sz,sz]) );
		assert( isdiag(matB) );
		vecB = diag(matB);
		clear matB;
	endif
	assert( min(vecB) > 0.0 );
	%
	% Apply scaling.
	vecBInv = 1.0./vecB;
	vecGScl = vecBInv .* vecG;
	matHScl = vecBInv .* matH .* (vecBInv');
	matHScl = (matHScl'+matHScl)/2.0;
	clear vecG;
	clear matH;
	clear vecB;
	%
	% Calculate eigenfactorization.
	% ( This may be too slow for some applications. )
	[ matPsi, matLambdaScl ] = eig(matHScl);
	vecGamma = -(matPsi'*vecGScl);
	vecLambdaScl = diag(matLambdaScl);
	clear matLambdaScl;
	clear vecGScl;
	clear matHScl;
	%
	vecLambdaMod = __findLambdaMod( f0, vecGamma, vecLambdaScl, prm );
	clear vecLamdbaScl;
	vecX = __findVecX( vecGamma, vecLambdaMod, matPsi, vecBInv, bMax, xMax );
return;
endfunction

function vecX = __findVecX( vecGamma, vecLambda, matPsi, vecBInv, xSclMax, xMax );
	s1 = 1.0;
	if ( ~isempty(xSclMax) )
		xSclNorm1 = __xSclNormOfS( s1, vecGamma, vecLambda );
		if ( xSclNorm1 > xSclMax )
			s1 = fzero( @(s) __xSclNormOfS( s, vecGamma, vecLambda ) - xSclMax, [ 0.0, s1 ] );
		endif
		clear xSclNorm1;
	endif
	if ( ~isempty(xMax) )
		xNorm1 = __xNormOfS( s1, vecGamma, vecLambda, vecBInv, matPsi );
		if ( xNorm1 > xMax )
			s1 = fzero( @(s) __xNormOfS( s, vecGamma, vecLambda, vecBInv, matPsi ) - xMax, [ 0.0, s1 ] );
		endif
		clear xNorm1;
	endif
	vecX = __vecXOfS( s1, vecGamma, vecLambda, vecBInv, matPsi );
return;
endfunction

function vecLambdaMod = __findLambdaMod( f0, vecGamma, vecLambda, prm );
	vecLambdaMod = vecLambda; % But, may be modified.
	% Do we even need to modify lambda?
	if ( min(vecLambda) > 0.0 )
		fCrit = f0 - (( (vecGamma.^2)./vecLambda )/2.0);
		if ( fCrit > -0.1*f0 )
			return;
		endif
	endif
	%
	% Let's get some bounds.
	lambdaHi = 2.0 * sumsq(vecGamma) / f0;
	if ( lambdaHi > max(vecLambda) )
		vecLambdaMod(:) = lambdaHi;
		assert( abs(f0 - (( (vecGamma.^2)./vecLambdaMod )/2.0)) < sqrt(eps)*f0 );
		return;
	endif
	lambdaLo = sqrt(eps) * max(abs(vecLambda));
	assert( lambdaLo > 0.0 );
	if ( __fCritOfLambdaFloor( lambdaLo, f0, vecGamma, vecLambda ) >= 0.0 )
		vecLambdaMod( vecLambda < lambdaLo ) = lambdaLo;
		assert( abs(f0 - (( (vecGamma.^2)./vecLambdaMod )/2.0)) < sqrt(eps)*f0 );
		return;
	endif
	%
	% Let's do a 1D solve.
	lambdaFloor = fzero( @(lamf) __fCritOfLambdaFloor(lamf,f0,vecGamma,vecLamba), [ lambdaLo, lambdaHi ] );
	vecLambdaMod( vecLambda < lambdaFloor ) = lambdaFloor;
	assert( abs(f0 - (( (vecGamma.^2)./vecLambdaMod )/2.0)) < sqrt(eps)*f0 );
return;
endfunction

function fCrit = __fCritOfLambdaFloor( lambdaFloor, f0, vecGamma, vecLambda )
	assert( lambdaFloor > 0.0 );
	vecLambda( vecLambda < lambdaFloor) = lambdaFloor;
	fCrit = f0 - (( (vecGamma.^2)./vecLambda )/2.0);
return;
endfunction

function xSclNorm = __xSclNormOfS( s, vecGamma, vecLambda )
	if ( 0.0 == s )
		xSclNorm = 0.0;
	else
		mu = max(vecLambda) * ( (1.0/s) - 1.0 );
		xSclNorm = sqrt(sumsq( vecGamma ./ ( vecLambda + mu ) ));
	endif
return;
endfunction

function xNorm = __xNormOfS( s, vecGamma, vecLambda, vecBInv, matPsi );
	xNorm = sqrt(sumsq( __vecXOfS( s, vecGamma, vecLambda, vecBInv, matPsi ) ));
return;
endfunction

function vecX = __vecXOfS( s, vecGamma, vecLambda, vecBInv, matPsi );
	if ( 0.0 == s )
		vecX = zeros(size(vecBInv));
	else
		mu = max(vecLambda) * ( (1.0/s) - 1.0 );
		vecX = vecBInv .* ( matPsi * ( vecGamma ./ ( vecLambda + mu ) ) );
	endif
return;
endfunction
