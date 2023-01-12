% function [ vecDelta, datOut ] = levsol0110( f0, vecG, matH, vecS=[], prm=[] )
%
% Going with simplicity, here.

function [ vecDelta, datOut ] = levsol0110( f0, vecG, matH, vecS=[], sMax=[], dMax=[], prm=[] )
	datOut = [];
	__validateInput( f0, vecG, matH, vecS, sMax, dMax, prm );
	%
	if ( isempty(vecS) )
		vecGScl = vecG;
		matHScl = matH;
	else
		vecGScl = vecG ./ vecS;
		matHScl = symm( (matH ./ vecS) ./ (vecS') ); % Autobroadcast.
	endif
	[ matPsi, matLambdaOrig ] = eig( matHScl );
	vecLambdaOrig = diag( matLambdaOrig );
	vecGamma = matPsi' * (-vecGScl);
	%
	vecLambdaMod = __findVecLambdaMod( f0, vecGamma, vecLambdaOrig, prm );
	assert( __fOfP( 1.0, f0, vecGamma, vecLambdaMod, vecLambdaMod ) >= -sqrt(eps)*f0 );
	vecDelta = __findVecDelta( f0, vecGamma, vecLambdaMod, vecLambdaOrig, matPsi, vecS, sMax, dMax, prm );
return;
endfunction


function __validateInput( f0, vecG, matH, vecS, sMax, dMax, prm )
	sz = size(vecG,1);
	assert( isrealscalar(f0) );
	assert( f0 > 0.0 );
	assert( isrealarray(vecG,[sz,1]) );
	assert( norm(vecG) > 0.0 );
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	if (~isempty(vecS))
		assert( isrealarray(vecS,[sz,1]) );
		assert( min(vecS) > 0.0 );
	endif
	if (~isempty(sMax))
		assert( isrealscalar(sMax) );
		assert( sMax >= 0.0 );
	endif
	if (~isempty(dMax))
		assert( isrealscalar(dMax) );
		assert( dMax >= 0.0 );
	endif
return;
endfunction


function vecLambdaMod = __findVecLambdaMod( f0, vecGamma, vecLambdaOrig, prm )
	if ( min(vecLambdaOrig) > 0.0 )
	if ( __fModCritOfLambdaFloor( 0.0, f0, vecGamma, vecLambdaOrig ) >= 0.0 )
		vecLambdaMod = vecLambdaOrig;
		return;
	endif
	endif
	%
	lambdaHi = sumsq(vecGamma) / ( 2.0 * f0 );
	if ( lambdaHi > max(vecLambdaOrig) )
		vecLambdaMod = vecLambdaOrig;
		vecLambdaMod(:) = lambdaHi;
		return;
	endif
	%
	lambdaLo = sqrt(eps) * max(abs(vecLambdaOrig));
	assert( lambdaLo > 0.0 );
	if ( __fModCritOfLambdaFloor( lambdaLo, f0, vecGamma, vecLambdaOrig ) >= 0.0 )
		vecLambdaMod = vecLambdaOrig;
		vecLambdaMod( vecLambdaOrig < lambdaLo ) = lambdaLo;
		return;
	endif	
	%
	lambdaFloor = fzerowrap( @(lamf) __fModCritOfLambdaFloor( lamf, f0, vecGamma, vecLambdaOrig ), [ lambdaLo, lambdaHi ] );
	vecLambdaMod = vecLambdaOrig;
	vecLambdaMod( vecLambdaOrig < lambdaFloor ) = lambdaFloor;
return;
endfunction


function vecDelta = __findVecDelta( f0, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS, sMax, dMax, prm )
	p1 = 1.0;
	if ( norm(__vecPhiOfP( p1, vecGamma, vecLambdaCurve )) > sMax )
		p1 = fzerowrap( @(p) (norm(__vecPhiOfP( p, vecGamma, vecLambdaCurve )) - sMax), [ 0.0, p1 ] );
	endif
	if ( norm(__vecDeltaOfP( p1, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS )) > dMax )
		p1 = fzerowrap( @(p) (norm(__vecDeltaOfP( p, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS )) - dMax), [ 0.0, p1 ] );
	endif
	if ( __fOfP( 1.0, f0, vecGamma, vecLambdaCurve, vecLambdaFunc ) < -0.01 * f0 ...
	  && __fOfP( p1, f0, vecGamma, vecLambdaCurve, vecLambdaFunc ) < 0.0 )
		p1 = fzerowrap( @(p) __fOfP( p, f0, vecGamma, vecLambdaCurve, vecLambdaFunc ), [ 0.0, p1 ] );
	endif
	vecDelta = __vecDeltaOfP( p1, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS );
return
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecPhi = __vecPhiOfP( p, vecGamma, vecLambdaCurve )
	if ( 0.0 == p )
		vecPhi = zeros(size(vecGamma));
	else
		assert( min(vecLambdaCurve) > 0.0 );
		mu = mean(vecLambdaCurve) * ( (1.0/p) - 1.0 );
		vecPhi = vecGamma ./ ( vecLambdaCurve + mu );
		% Alterlative scaling of mu might make things easier for 1D solver.
	endif
return;
endfunction

function f = __fOfP( p, f0, vecGamma, vecLambdaCurve, vecLambdaFunc )
	vecPhi = __vecPhiOfP( p, vecGamma, vecLambdaCurve );
	f = f0 - ( vecGamma' * vecPhi ) + (( vecPhi' * (vecPhi.*vecLambdaFunc) )/2.0);
return;
endfunction

function vecDelta = __vecDeltaOfP( p, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS )
	vecDelta = ( matPsi * __vecPhiOfP( p, vecGamma, vecLambdaCurve ) ) ./ vecS;
return;
endfunction

function fModCrit = __fModCritOfLambdaFloor( lambdaFloor, f0, vecGamma, vecLambdaOrig )
	vecLambdaMod = vecLambdaOrig;
	vecLambdaMod( vecLambdaOrig < lambdaFloor ) = lambdaFloor;
	fModCrit = __fOfP( 1.0, f0, vecGamma, vecLambdaMod, vecLambdaMod );
	% Note: We use Mod in place of Orig in the calculation of f because we want "f *MOD* Crit".
return;
endfunction
