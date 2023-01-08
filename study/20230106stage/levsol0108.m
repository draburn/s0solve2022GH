function vecX = levsol0108( vecG, matH, matB, bMax=[], prm=[] )
	sz = size(vecG,1);
	assert( isrealarray(vecG,[sz,1]) );
	assert( isrealarray(matH,[sz,sz]) );
	assert( issymmetric(matH) );
	assert( isrealarray(matB,[sz,sz]) );
	assert( isdiag(matB) );
	vecB = diag(matB);
	assert( min(vecB) > 0.0 );
	vecBInv = 1.0./vecB;
	vecGScl = vecBInv .* vecG;
	matHScl = vecBInv .* matH .* (vecBInv');
	vecXScl = __solve( vecGScl, matHScl, bMax, prm );
	if ( isempty(vecXScl) )
		vecX = [];
	else
		vecX = vecBInv .* vecXScl;
	endif
return;
endfunction

function vecX = __solve( vecG, matH, bMax, prm )
	sz = size(vecG,1);
	[ matPsi, matLambda_original ] = eig(matH);
	vecLambda = diag(matLambda_original); % But, we'll modify this below.
	eigMax = max(vecLambda);
	if ( eigMax <= 0.0 )
		% We won't bother with this atypical case.
		vecX = [];
		return;
	endif
	epsEig = mygetfield( prm, "epsEig", sqrt(eps) );
	assert( isrealscalar(epsEig) );
	assert( epsEig > 0.0 );
	assert( epsEig < 1.0 );
	vecLambda( vecLambda < epsEig*eigMax ) = epsEig*eigMax;
	vecGamma = matPsi' * vecG;
	vecXNewt = matPsi * ( -vecGamma./vecLambda );
	%
	if ( isempty(bMax) )
		% We're not optimizing for this case, but, here you go.
		vecX = vecXNewt;
		return;
	endif
	assert( isrealscalar(bMax) );
	if ( 0.0 == bMax )
		% This is stupid, but, meh.
		vecX = zeros(sz,1);
		return;
	endif
	assert( bMax > 0.0 );
	if ( norm(vecXNewt) <= bMax )
		vecX = vecXNewt;
		return;
	endif
	%
	% Set up 1D solver.
	hScl = sqrt( sum(vecLambda.^2) / sz );
	funch_resOfS = @(s)(  __bOfS( s, hScl, vecGamma, vecLambda ) - bMax );
	assert( funch_resOfS(0.0) < 0.0 );
	assert( funch_resOfS(1.0) > 0.0 );
	s = fzero( funch_resOfS, [ 0.0, 1.0 ] );
	if ( s < 0.0 )
		s
		error( "s is negative" );
	elseif ( 0.0 == s )
		vecX = zeros(sz,1);
		warning( "s is zero." );
	elseif ( s < 1.0 )
		mu = hScl * ( (1.0/s) - 1.0 );
		vecX = matPsi * ( -vecGamma ./ ( vecLambda + mu ) );
	elseif ( 1.0 == s )
		vecX = vecXNewt;
		warning( "s is unity." );
	else
		s
		error( "s is greater than unity" );
	endif
	
	assert( s > 0.0 );
	assert( s < 1.0 );
	% Constraints on s might fail if problem is extremely poorly poosed.
return;
endfunction


function b = __bOfS( s, hScl, vecGamma, vecLambda )
	if ( 0.0 == s )
		b = 0.0;
	else
		mu = hScl * ( (1.0/s) - 1.0 );
		b = sqrt(sumsq( vecGamma ./ ( vecLambda + mu ) ));
	endif
return;
endfunction

%!test
%!	clear;
%!	setprngstates(0);
%!	sizeX = 5;
%!	matH = mtm(randn(sizeX,sizeX));
%!	vecXCrit = randn(sizeX,1)
%!	matB = diag(abs(randn(sizeX,1)));
%!	%
%!	vecX0 = zeros(sizeX,1);
%!	vecG = matH*(vecX0-vecXCrit);
%!	bInf = norm(matB*(vecXCrit-vecX0))
%!	%
%!	vecDeltaUndef = levsol0108( vecG, matH, matB, [] );
%!	gCheck = norm( vecG + matH*(vecX0+vecDeltaUndef) )
%!	vecDeltaInfP = levsol0108( vecG, matH, matB, bInf+0.01 );
%!	vecDeltaInf = levsol0108( vecG, matH, matB, bInf );
%!	vecDeltaInfM = levsol0108( vecG, matH, matB, bInf-0.01 );
%!	%
%!	mat_nearEnd = [ vecDeltaUndef, vecDeltaInfP, vecDeltaInf, vecDeltaInfM ]
%!	vecDelta06 = levsol0108( vecG, matH, matB, bInf*0.6 );
%!	vecDelta05 = levsol0108( vecG, matH, matB, bInf*0.5 );
%!	vecDelta04 = levsol0108( vecG, matH, matB, bInf*0.4 );
%!	mat_nearMid = [ vecDelta06, vecDelta05, vecDelta04 ]
%!	%
%!	vecDelta002 = levsol0108( vecG, matH, matB, bInf*0.02 );
%!	vecDelta001 = levsol0108( vecG, matH, matB, bInf*0.01 );
%!	vecDelta000 = levsol0108( vecG, matH, matB, bInf*0.00 );
%!	mat_nearZer = [ vecDelta002, vecDelta001, vecDelta000 ]
%!	%
%!	msg( __FILE__, __LINE__, "Please check the above results for reasonableness." );
