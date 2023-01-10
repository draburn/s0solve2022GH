% function vecX = levsol0109( f0, vecG, matH, matB=[], bMax=[], xMax=[], prm=[] )
% Uses "lambdaFloor" for positive-defness and fMin.
% Allows for constrants on both ||B*x|| and ||x||.
% ... It would probably feel more natural to specify D for curve and some set of Bs for boundaries.
% Also, for simplicity, should REALLY just do everything in terms of delta(mu).

function vecX = levsol0109( f0, vecG, matH, matB=[], bMax=[], xMax=[], prm=[] )
	% Validate input.
	sz = size(vecG,1);
	assert( isrealscalar(f0) );
	assert( f0 > 0.0 );
	assert( isrealarray(vecG,[sz,1]) );
	assert( norm(vecG) > 0.0 );
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	%assert( sum(sum(abs(matH))) > 0.0 );
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
	if (0)
		msg( __FILE__, __LINE__, "Infodump..." );
		matHMod = diag(1.0./vecBInv) * matPsi * diag(vecLambdaMod) * (matPsi') * diag(1.0./vecBInv)
		[ vecLambdaScl, vecLambdaMod ]
		fCrit_modmod = f0 - sum( (vecLambdaMod-(vecLambdaScl/2.0)) .* ((vecGamma./vecLambdaMod).^2)  )
		fCrit = f0 - (sum( (vecGamma.^2)./vecLambdaMod )/2.0)
		vecX = [];
		return
	endif
	vecX = __findVecX( vecGamma, vecLambdaMod, matPsi, vecBInv, bMax, xMax, f0, vecLambdaScl, prm );
return;
endfunction

function vecX = __findVecX( vecGamma, vecLambdaC, matPsi, vecBInv, xSclMax, xMax, f0, vecLambdaF, prm );
	s1 = 1.0;
	f1 = __fOfS( s1, vecGamma, vecLambdaC, f0, vecLambdaF );
	if ( f1 < -0.1*f0 )
		% This constraint is typically the last important;
		% but, we only want to trigger it if f is *sufficiently* negative,
		% with "sufficiently" not depending on any other criteria.
		%msgnnl( __FILE__, __LINE__, "Reducing s1 from " );
		%printf( "%g to ", s1 );
		s1 = __fzeroWrapper( @(s) __fOfS( s, vecGamma, vecLambdaC, f0, vecLambdaF ), [ 0.0, s1 ], prm );
		%printf( "%g.\n", s1 );
		%f1 = __fOfS( s1, vecGamma, vecLambdaC, matPsi, f0, vecLambdaF )
	endif
	clear f1;
	if ( ~isempty(xSclMax) )
		xSclNorm1 = __xSclNormOfS( s1, vecGamma, vecLambdaC );
		if ( xSclNorm1 > xSclMax )
			%s1 = fzero( @(s) __xSclNormOfS( s, vecGamma, vecLambda ) - xSclMax, [ 0.0, s1 ] );
			%msgnnl( __FILE__, __LINE__, "Reducing s1 from " );
			%printf( "%g to ", s1 );
			s1 = __fzeroWrapper( @(s) __xSclNormOfS( s, vecGamma, vecLambdaC ) - xSclMax, [ 0.0, s1 ], prm );
			%printf( "%g.\n", s1 );
			%xSclNorm1 = __xSclNormOfS( s1, vecGamma, vecLambdaC )
		endif
		clear xSclNorm1;
	endif
	if ( ~isempty(xMax) )
		xNorm1 = __xNormOfS( s1, vecGamma, vecLambdaC, vecBInv, matPsi );
		if ( xNorm1 > xMax )
			%msgnnl( __FILE__, __LINE__, "Reducing s1 from " );
			%printf( "%g to ", s1 );
			%s1 = fzero( @(s) __xNormOfS( s, vecGamma, vecLambda, vecBInv, matPsi ) - xMax, [ 0.0, s1 ] );
			s1 = __fzeroWrapper( @(s) __xNormOfS( s, vecGamma, vecLambdaC, vecBInv, matPsi ) - xMax, [ 0.0, s1 ], prm );
			%printf( "%g.\n", s1 );
		endif
		clear xNorm1;
	endif
	vecX = __vecXOfS( s1, vecGamma, vecLambdaC, vecBInv, matPsi );
return;
endfunction

function vecLambdaMod = __findLambdaMod( f0, vecGamma, vecLambda, prm );
	vecLambdaMod = vecLambda; % But, may be modified.
	% Do we even need to modify lambda?
	if ( min(vecLambda) > 0.0 )
		fCrit = f0 - (sum( (vecGamma.^2)./vecLambdaMod )/2.0);
		if ( fCrit > -0.1*f0 )
			return;
		endif
	endif
	%
	% Let's get some bounds.
	% This hi bound should handle the case H = zero, FWIW.
	lambdaHi = sumsq(vecGamma) / ( 2.0 * f0 );
	if ( lambdaHi > max(vecLambda) )
		vecLambdaMod(:) = lambdaHi;
		fCrit = f0 - (sum( (vecGamma.^2)./vecLambdaMod )/2.0);
		assert( abs(fCrit) < sqrt(eps)*f0 );
		return;
	endif
	% This lo bound may push fCrit positive,
	%  but is about as safely as we can go while making H positive definite,
	%  which is good, because then we don't have to consider a positive-semi-definite case.
	lambdaLo = sqrt(eps) * max(abs(vecLambda));
	assert( lambdaLo > 0.0 );
	if ( __fCritOfLambdaFloor( lambdaLo, f0, vecGamma, vecLambda ) >= 0.0 )
		vecLambdaMod( vecLambda < lambdaLo ) = lambdaLo;
		fCrit = f0 - (sum( (vecGamma.^2)./vecLambdaMod )/2.0);
		assert( fCrit > -sqrt(eps)*f0 );
		return;
	endif
	%
	% Let's do a 1D solve.
	%lambdaFloor = fzero( @(lamf) __fCritOfLambdaFloor(lamf,f0,vecGamma,vecLambda), [ lambdaLo, lambdaHi ] );
	lambdaFloor = __fzeroWrapper( @(lamf) __fCritOfLambdaFloor(lamf,f0,vecGamma,vecLambda), [ lambdaLo, lambdaHi ], prm );
	vecLambdaMod( vecLambda < lambdaFloor ) = lambdaFloor;
	fCrit = f0 - (sum( (vecGamma.^2)./vecLambdaMod )/2.0);
	assert( abs(fCrit) < sqrt(eps)*f0 );
return;
endfunction

function x = __fzeroWrapper( funchF, xVals, prm=[] )
	% Built-in fzero sometimes returns worse of two bounds...
	[ fzero_x, fzero_fval, fzero_info, fzero_output ] = fzero( funchF, xVals );
	xOut = fzero_output.bracketx;
	fOut = fzero_output.brackety;
	if ( abs(fOut(2)-fOut(1)) > sqrt(eps)*(abs(fOut(2))+abs(fOut(1))) )
		x = ( xOut(1)*fOut(2) - xOut(2)*fOut(1) ) / ( fOut(2) - fOut(1) );
	else
		x = fzero_x;
	endif
return;
endfunction

function fCrit = __fCritOfLambdaFloor( lambdaFloor, f0, vecGamma, vecLambda )
	assert( lambdaFloor > 0.0 );
	vecLambda( vecLambda < lambdaFloor) = lambdaFloor;
	fCrit = f0 - (sum( (vecGamma.^2)./vecLambda )/2.0);
	%%%vecLambdaMod = vecLambda;
	%%%vecLambdaMod( vecLambda < lambdaFloor) = lambdaFloor;
	%%%fCrit = f0 - sum( (vecLambdaMod-(vecLambda/2.0)) .* ((vecGamma./vecLambdaMod).^2)  );
return;
endfunction

function xSclNorm = __xSclNormOfS( s, vecGamma, vecLambda )
	if ( 0.0 == s )
		xSclNorm = 0.0;
	else
		%mu = max(vecLambda) * ( (1.0/s) - 1.0 );
		mu = sqrt(mean(vecLambda.^2)) * ( (1.0/s) - 1.0 );
		xSclNorm = sqrt(sumsq( vecGamma ./ ( vecLambda + mu ) ));
	endif
return;
endfunction

function xNorm = __xNormOfS( s, vecGamma, vecLambda, vecBInv, matPsi );
	xNorm = sqrt(sumsq( __vecXOfS( s, vecGamma, vecLambda, vecBInv, matPsi ) ));
return;
endfunction

function f = __fOfS( s, vecGamma, vecLambdaC, f0, vecLambdaF );
	if ( 0.0 == s )
		f = f0;
	else
		%mu = max(vecLambdaC) * ( (1.0/s) - 1.0 );
		mu = sqrt(mean(vecLambdaC.^2)) * ( (1.0/s) - 1.0 );
		vecLCPMI = vecLambdaC + mu;
		vecGLInv = vecGamma./vecLCPMI;
		f = f0 - sum(vecGamma.*vecGLInv) + sum(vecGLInv.*vecLambdaF.*vecGLInv)/2.0;
	endif
return;
endfunction

function vecX = __vecXOfS( s, vecGamma, vecLambda, vecBInv, matPsi );
	if ( 0.0 == s )
		vecX = zeros(size(vecBInv));
	else
		%mu = max(vecLambda) * ( (1.0/s) - 1.0 );
		mu = sqrt(mean(vecLambda.^2)) * ( (1.0/s) - 1.0 );
		vecX = vecBInv .* ( matPsi * ( vecGamma ./ ( vecLambda + mu ) ) );
	endif
return;
endfunction
