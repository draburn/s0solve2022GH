% function [ vecDelta, datOut ] = eigfloorsol0111( f0, vecG, matH, vecS=[], sMax=[], dMax=[], prm=[] )
%
% DRaburn 2023-01-11

function [ vecDelta, datOut ] = eigfloorsol0111( f0, vecG, matH, vecS=[], sMax=[], dMax=[], prm=[] )
	datOut = [];
	[ fMin, fModMin ] = __init( f0, vecG, matH, vecS, sMax, dMax, prm );
	%
	if ( isempty(vecS) )
		vecGScl = vecG;
		matHScl = matH;
		vecS = ones(size(vecG));
	else
		vecGScl = vecG ./ vecS;
		matHScl = symm( (matH ./ vecS) ./ (vecS') ); % Autobroadcast.
	endif
	%
	[ matPsi, matLambda ] = eig( matHScl );
	vecLambda = diag( matLambda );
	vecGamma = matPsi' * (-vecGScl);
	assert( norm(vecGamma) > 0.0 );
	%
	if ( min(vecLambda) > 0.0 )
		mu0 = min(vecLambda);
		[ vecDelta, mu ] = __findVecDeltaPD( f0, vecGamma, vecLambda, matPsi, vecS, mu0, sMax, dMax, fMin, fModMin, prm );
	else
		if ( isempty(sMax) && isempty(dMax) && isempty(fMin) && isempty(fModMin) )
			msg( __FILE__, __LINE__, "ERROR: matH is not posdef and there are no constraitns." );
			msg( __FILE__, __LINE__, "  There is no (supported) solution in this case." );
			msg( __FILE__, __LINE__, "  Note that 'perfectly balanced' solutions are not currently supported." );
			vecDelta = [];
			mu = [];
			return;
		endif
		if ( max(vecLambda) <= 0.0 )
			[ vecDelta, mu ] = __findVecDeltaNSD( f0, vecGamma, vecLambda, matPsi, vecS, sMax, dMax, fMin, fModMin, prm );
		else
			% This case is hackish, but should typically be fine.
			mu0 = sqrt(eps) * max(abs(vecLambda));
			[ vecDelta, mu ] = __findVecDeltaPD( f0, vecGamma, vecLambda, matPsi, vecS, mu0, sMax, dMax, fMin, fModMin, prm );
		endif
	endif
	%
	vecLambdaMod = vecLambda;
	vecLambdaMod( vecLambdaMod < mu ) = mu;
	matHMod = symm( diag(vecS) * matPsi * diag(vecLambdaMod) * (matPsi') * (diag(vecS)') );
	datOut.fOrigPred = f0 + vecDelta'*vecG + (vecDelta'*matH*vecDelta)/2.0;
	datOut.fModPred = f0 + vecDelta'*vecG + (vecDelta'*matHMod*vecDelta)/2.0;
	datOut.vecGOrigPred = vecG + matH * vecDelta;
	datOut.vecGModPred = vecG + matHMod * vecDelta;
	datOut.matHMod = matHMod;
	
	if (0)
		msg( __FILE__, __LINE__, "BEGIN INFODUMP..." );
		vecLambda'
		mu
		d = norm(vecDelta)
		s = norm(vecS.*vecDelta)
		f0
		fOrigPred = f0 + vecDelta'*vecG + (vecDelta'*matH*vecDelta)/2.0
		fModPred = f0 + vecDelta'*vecG + (vecDelta'*matHMod*vecDelta)/2.0
		g0 = norm( vecG )
		gOrigPred = norm( vecG + matH * vecDelta )
		gModPred = norm( vecG + matHMod * vecDelta )
		gos0 = norm( vecG./vecS )
		gosOrigPred = norm( (vecG + matH * vecDelta)./vecS )
		gosModPred = norm( (vecG + matHMod * vecDelta)./vecS )
		%assert( min(vecLambda) > 0.0 )
		msg( __FILE__, __LINE__, "END INFODUMP." );
	endif
	
return;
endfunction


function [ fMin, fModMin ] = __init( f0, vecG, matH, vecS, sMax, dMax, prm )
	assert( 7 == nargin );
	sz = size(vecG,1);
	assert( isrealscalar(f0) );
	assert( f0 > 0.0 );
	assert( isrealarray(vecG,[sz,1]) );
	assert( norm(vecG) > 0.0 );
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	if ( ~isempty(vecS) )
		assert( isrealarray(vecS,[sz,1]) );
		assert( min(vecS) > 0.0 );
	endif
	if ( ~isempty(sMax) )
		assert( isrealscalar(sMax) );
		assert( sMax >= 0.0 );
	endif
	if ( ~isempty(dMax) )
		assert( isrealscalar(dMax) );
		assert( dMax >= 0.0 );
	endif
	fMin = mygetfield( prm, "fMin", -1.0E-4*f0 );
	if ( ~isempty(fMin) )
		assert( isrealscalar(fMin) );
		assert( fMin < f0 );
	endif
	fModMin = mygetfield( prm, "fModMin", -1.0E-4*f0 );
	if ( ~isempty(fModMin) )
		assert( isrealscalar(fModMin) );
		assert( fModMin < f0 );
	endif
	%
	if ( isfield(prm,"fMinRegu") )
		msg( __FILE__, __LINE__, "WARNING: prm has field 'fMinRegu'; did you mean 'fModMin'?" );
	endif
return;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ vecDelta, mu ] = __findVecDeltaNSD( f0, vecGamma, vecLambda, matPsi, vecS, sMax, dMax, fMin, fModMin, prm )
	assert( 10 == nargin );
	mu = [];
	if ( ~isempty(sMax) )
		% vecPhi = vecGamma / mu.
		mu = max([ mu, norm(vecGamma)/sMax ]);
	endif
	if ( ~isempty(dMax) )
		% vecDelta = ( matPsi * vecGamma ) / mu.
		mu = max([ mu, norm((matPsi*vecGamma)./vecS)/dMax ]);
	endif
	if ( ~isempty(fMin) )
		% muCrit = ( vecGamma' * matLambda * vecGamma ) / ( vecGamma' * vecGamma ) <= 0.0 b/c NSD.
		% f = f0 - ( vecGamma' * vecGamma ) / mu + ( vecGamma' * matLambda * vecGamma ) / ( 2.0 * mu^2 ).
		discrim = ( vecGamma' * vecGamma )^2 - (2.0 * ( vecGamma' * ( vecGamma .* vecLambda ) ) * ( f0 - fMin ));
		assert( discrim >= 0.0 );
		mu = max([ mu, ( (vecGamma'*vecGamma) + sqrt(discrim) ) / ( 2.0 * ( f0 - fMin ) ) ]);
	endif
	if ( ~isempty(fModMin) )
		% fMod = f0 - sumsq(vecGamma) / ( 2.0 * mu ).
		mu = max([ mu, sumsq(vecGamma) / ( 2.0 * ( f0 - fModMin ) ) ]);
	endif
	assert( ~isempty(mu) );
	vecDelta = __vecDeltaOfMu( mu, vecGamma, vecLambda, matPsi, vecS );
return;
endfunction


function [ vecDelta, mu ] = __findVecDeltaPD( f0, vecGamma, vecLambda, matPsi, vecS, mu0, sMax, dMax, fMin, fModMin, prm );
	assert( 11 == nargin );
	p1 = 1.0;
	if ( ~isempty(sMax) )
	if ( norm(__vecPhiOfP( p1, mu0, vecGamma, vecLambda)) > sMax )
		p1 = fzerowrap( @(p)norm(__vecPhiOfP( p, mu0, vecGamma, vecLambda)) - sMax, [ 0.0, p1 ] );
	endif
	endif
	if ( ~isempty(dMax) )
	if ( norm(__vecDeltaOfP( p1, mu0, vecGamma, vecLambda, matPsi, vecS )) > dMax )
		p1 = fzerowrap( @(p)norm(__vecDeltaOfP( p, mu0, vecGamma, vecLambda, matPsi, vecS )) - dMax, [ 0.0, p1 ] );
	endif
	endif
	if ( ~isempty(fMin) )
	if ( __fOfP( p1, mu0, f0, vecGamma, vecLambda ) < fMin )
		p1 = fzerowrap( @(p)__fOfP( p, mu0, f0, vecGamma, vecLambda ) - fMin, [ 0.0, p1 ] );
	endif
	endif
	if ( ~isempty(fModMin) )
	if ( __fModOfP( p1, mu0, f0, vecGamma, vecLambda ) < fModMin )
		p1 = fzerowrap( @(p)__fModOfP( p, mu0, f0, vecGamma, vecLambda ) - fModMin, [ 0.0, p1 ] );
	endif
	endif
	vecDelta = __vecDeltaOfP( p1, mu0, vecGamma, vecLambda, matPsi, vecS );
	mu = mu0/p1;
return;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecPhi =__vecPhiOfMu( mu, vecGamma, vecLambda )
	assert( 3 == nargin );
	vecLambdaMod = vecLambda;
	vecLambdaMod( vecLambda < mu ) = mu;
	vecPhi = vecGamma ./ vecLambdaMod;
return;
endfunction
function vecDelta =__vecDeltaOfMu( mu, vecGamma, vecLambda, matPsi, vecS )
	assert( 5 == nargin );
	vecPhi = __vecPhiOfMu( mu, vecGamma, vecLambda );
	vecDelta = ( matPsi * vecPhi ) ./ vecS;
return;
endfunction
function f = __fOfMu( mu, f0, vecGamma, vecLambda )
	assert( 4 == nargin );
	vecPhi = __vecPhiOfMu( mu, vecGamma, vecLambda );
	f = f0 - ( vecPhi' * vecGamma ) + (( vecPhi' * ( vecPhi.*vecLambda ) )/2.0);
return;
endfunction
function fMod = __fModOfMu( mu, f0, vecGamma, vecLambda )
	assert( 4 == nargin );
	vecLambdaMod = vecLambda;
	vecLambdaMod( vecLambda < mu ) = mu;
	vecPhi = vecGamma ./ vecLambdaMod;
	fMod = f0 - ( vecPhi' * vecGamma ) + (( vecPhi' * ( vecPhi.*vecLambdaMod ) )/2.0);
return;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecPhi =__vecPhiOfP( p, mu0, vecGamma, vecLambda )
	assert( 4 == nargin );
	if ( 0.0 == p )
		vecPhi = zeros(size(vecGamma));
	else
		vecPhi = __vecPhiOfMu( mu0/p, vecGamma, vecLambda );
	endif
return;
endfunction
function vecDelta =__vecDeltaOfP( p, mu0, vecGamma, vecLambda, matPsi, vecS )
	assert( 6 == nargin );
	if ( 0.0 == p )
		vecDelta = zeros(size(vecGamma));
	else
		vecDelta = __vecDeltaOfMu( mu0/p, vecGamma, vecLambda, matPsi, vecS );
	endif
return;
endfunction
function f = __fOfP( p, mu0, f0, vecGamma, vecLambda )
	assert( 5 == nargin );
	if ( 0.0 == p )
		f = f0;
	else
		f = __fOfMu( mu0/p, f0, vecGamma, vecLambda );
	endif
return;
endfunction
function fMod = __fModOfP( p, mu0, f0, vecGamma, vecLambda )
	assert( 5 == nargin );
	if ( 0.0 == p )
		fMod = f0;
	else
		fMod = __fModOfMu( mu0/p, f0, vecGamma, vecLambda );
	endif
return;
endfunction
