%function [ vecDelta, datOut ] = levsol0111( f0, vecG, matH, vecS=[], sMax=[], dMax=[], prm=[] )
%
% DRaburn 2023-01-11
%  This version includes fMin and fMinRegu.

function [ vecDelta, datOut ] = levsol0111( f0, vecG, matH, vecS=[], sMax=[], dMax=[], prm=[] )
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
	%
	[ matPsi, matLambdaOrig ] = eig( matHScl );
	vecLambdaOrig = diag( matLambdaOrig );
	vecGamma = matPsi' * (-vecGScl);
	%
	vecLambdaMod = __findVecLambdaMod( f0, vecGamma, vecLambdaOrig, prm );
	vecDelta = __findVecDelta( f0, vecGamma, vecLambdaMod, vecLambdaOrig, matPsi, vecS, sMax, dMax, prm );
return;
endfunction


function __validateInput( f0, vecG, matH, vecS, sMax, dMax, prm )
	assert( 7 == nargin );
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
	%
	if ( isfield(prm,"fModMin") )
		msg( __FILE__, __LINE__, "WARNING: prm has field 'fModMin'; did you mean 'fMinRegu'?" );
	endif
return;
endfunction


function vecLambdaMod = __findVecLambdaMod( f0, vecGamma, vecLambdaOrig, prm )
	assert( 4 == nargin );
	fMinRegu = mygetfield( prm, "fMinRegu", 0.0 );
	if ( min(vecLambdaOrig) > 0.0 )
	if ( isempty(fMinRegu) || __fModCritOfLambdaFloor( 0.0, f0, vecGamma, vecLambdaOrig ) >= fMinRegu ) % Short-circuit.
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
	if ( isempty(fMinRegu) || __fModCritOfLambdaFloor( lambdaLo, f0, vecGamma, vecLambdaOrig ) >= fMinRegu ) % Short-circuit.
		vecLambdaMod = vecLambdaOrig;
		vecLambdaMod( vecLambdaOrig < lambdaLo ) = lambdaLo;
		return;
	endif
	assert( fMinRegu < f0 );
	%
	lambdaFloor = fzerowrap( @(lamf) __fModCritOfLambdaFloor( lamf, f0, vecGamma, vecLambdaOrig ) - fMinRegu, [ lambdaLo, lambdaHi ] );
	vecLambdaMod = vecLambdaOrig;
	vecLambdaMod( vecLambdaOrig < lambdaFloor ) = lambdaFloor;
return;
endfunction


function vecDelta = __findVecDelta( f0, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS, sMax, dMax, prm )
	assert( 9 == nargin );
	p1 = 1.0;
	if ( ~isempty( sMax ) )
	if ( norm(__vecPhiOfP( p1, vecGamma, vecLambdaCurve )) > sMax )
		p1 = fzerowrap( @(p) (norm(__vecPhiOfP( p, vecGamma, vecLambdaCurve )) - sMax), [ 0.0, p1 ] );
	endif
	endif
	if ( ~isempty(dMax) )
	if ( norm(__vecDeltaOfP( p1, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS )) > dMax )
		p1 = fzerowrap( @(p) (norm(__vecDeltaOfP( p, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS )) - dMax), [ 0.0, p1 ] );
	endif
	endif
	fMin = mygetfield( prm, "fMin", 0.0 );
	if ( ~isempty(fMin) )
	assert( fMin < f0 );
	if ( __fOfP( p1, f0, vecGamma, vecLambdaCurve, vecLambdaFunc ) < fMin )
		p1 = fzerowrap( @(p) __fOfP( p, f0, vecGamma, vecLambdaCurve, vecLambdaFunc ), [ 0.0, p1 ] );
	endif
	endif
	vecDelta = __vecDeltaOfP( p1, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS );
return
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecPhi = __vecPhiOfP( p, vecGamma, vecLambdaCurve )
	assert( 3 == nargin );
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
	assert( 5 == nargin );
	vecPhi = __vecPhiOfP( p, vecGamma, vecLambdaCurve );
	f = f0 - ( vecGamma' * vecPhi ) + (( vecPhi' * (vecPhi.*vecLambdaFunc) )/2.0);
return;
endfunction

function vecDelta = __vecDeltaOfP( p, vecGamma, vecLambdaCurve, vecLambdaFunc, matPsi, vecS )
	assert( 6 == nargin );
	vecDelta = ( matPsi * __vecPhiOfP( p, vecGamma, vecLambdaCurve ) ) ./ vecS;
return;
endfunction

function fModCrit = __fModCritOfLambdaFloor( lambdaFloor, f0, vecGamma, vecLambdaOrig )
	assert( 4 == nargin );
	vecLambdaMod = vecLambdaOrig;
	vecLambdaMod( vecLambdaOrig < lambdaFloor ) = lambdaFloor;
	fModCrit = __fOfP( 1.0, f0, vecGamma, vecLambdaMod, vecLambdaMod );
	% Note: We use Mod in place of Orig in the calculation of f because we want "f *MOD* Crit".
return;
endfunction
