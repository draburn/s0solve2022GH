% Function ...
function [ vecX, datOut ] = findLocMin_alytJ_mk2( vecX0, funchFJ, prm=[] )
	% Use retcode?
	%
	%
	% Parse input.
	sizeX = size(vecX0,1);
	%
	debugMode = mygetfield( prm, "debugMode", true );
	omegaMin = mygetfield( prm, "omegaMin", 0.0 );
	gNormTol = mygetfield( prm, "gNormTol", 0.0 );
	iterLimit = mygetfield( prm, "iterLimit", 100 );
	omegaFallTol = mygetfield( prm, "omegaFallTol", sqrt(eps) );
	stepSizeTol = mygetfield( prm, "stepSizeTol", sqrt(eps) );
	hessianType = mygetfield( prm, "hessianType", 0 );
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealscalar(omegaMin) );
		assert( isrealscalar(gNormTol) );
		assert( isrealscalar(iterLimit) );
		assert( isrealscalar(omegaFallTol) );
		assert( isrealscalar(stepSizeTol) );
		assert( isrealscalar(hessianType) );
	endif
	%
	[ vecF0, matJ0 ] = funchFJ( vecX0 );
	sizeF = size(vecF0,1);
	if (debugMode)
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		vecDelta_test = 0.001*ones(sizeX,1);
		vecJF_testA = ( funchFJ( vecX0 + vecDelta_test ) - funchFJ( vecX0 - vecDelta_test ) ) / 2.0;
		vecJF_testB = matJ0*vecDelta_test;
		if ( sumsq(vecJF_testA-vecJF_testB) > sqrt(eps)*(sumsq(vecJF_testA)+sumsq(vecJF_testA)) )
			msg( __FILE__, __LINE__, "*** WARNING: Jacobian calculated by funchFJ appears to be incorrect. ***" );
		endif
		clear vecDelta_test;
		clear vecJF_testA;
		clear vecJF_testB;
	endif
	%
	omega0 = sumsq(vecF0,1)/2.0;
	vecG0 = matJ0'*vecF0;
	matIX = eye(sizeX,sizeX);
	if ( nargout >= 2 )
		datOut = [];
	end
	%
	%
	% Prep main loop.
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	omega = omega0;
	vecG = vecG0;
	iterCount = 0;
	while (1)
		% Check pre-iter success conditions.
		if ( omega <= omegaMin )
			msgif( debugMode, __FILE__, __LINE__, "Reached omegaMin." );
			return;
		elseif ( norm(vecG) <= gNormTol )
			msgif( debugMode, __FILE__, __LINE__, "Reached gNormTol." );
			return;
		endif
		%
		%
		% Check pre-iter imposed limits.
		iterCount++;
		if ( iterCount > iterLimit )
			msgif( debugMode, __FILE__, __LINE__, "Reached iterLimit." );
			return;
		endif
		%
		%
		% Calculate Hessian.
		matH_jtj = matJ'*matJ;
		switch ( hessianType )
		case 0
			matH = matH_jtj;
		otherwise
			error( "Invalid hessianType." );
		endswitch
		%
		%
		%
		% Pick step.
		% REPLACE THIS WITH BACKTRACKING.
		mu = 0.0;
		[ matR, cholFlag ] = chol( matH + (mu*matIX) );
		if ( 0 ~= cholFlag )
			msg( __FILE__, __LINE__, "Cholesky factorization of Hessian failed." );
			return;
		end
		vecDelta = -( matR \ ( matR' \ vecG ) );
		%
		vecX_next = vecX + vecDelta;
		[ vecF_next, matJ_next ] = funchFJ( vecX_next );
		omega_next = sumsq(vecF_next)/2.0;
		if ( omega_next >= omega )
			msg( __FILE__, __LINE__, "Omega increased." );
			return;
		end
		vecG_next = matJ_next'*vecF_next;
		%
		%
		% Take step.
		vecX_prev = vecX;
		vecF_prev = vecF;
		matJ_prev = matJ;
		omega_prev = omega;
		vecG_prev = vecG;
		vecX = vecX_next;
		vecF = vecF_next;
		matJ = matJ_next;
		omega = omega_next;
		vecG = vecG_next;
		%
		%
		% Check post-iter imposed limits.
		if ( abs(omega-omega_prev) <= omegaFallTol )
			msg( __FILE__, __LINE__, "Reached omegaFallTol." );
			return;
		elseif ( norm(vecX-vecX_prev) <= stepSizeTol )
			msg( __FILE__, __LINE__, "Reached stepSizeTol." );
			return;
		endif
	endwhile
return;
end


%!function [ vecF, matJ ] = funcFJ_trivial( vecX )
%!	sizeX = size(vecX,1);
%!	vecF = vecX-1.0;
%!	matJ = eye(sizeX,sizeX);
%!endfunction
%!function [ vecF, matJ ] = funcFJ_easy( vecX )
%!	sizeX = size(vecX,1);
%!	matJ = diag((1:sizeX));
%!	vecF = matJ*( vecX - (1:sizeX)' );
%!	%vecF = vecX - (1:sizeX)';
%!endfunction


%!test
%!	caseNum = 1;
%!	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
%!	switch (caseNum)
%!	case 0
%!		sizeX = 2;
%!		sizeF = 2;
%!		vecX0 = zeros(sizeX,1);
%!		funchFJ = @(dummyX)( funcFJ_trivial(dummyX) );
%!	case 1
%!		sizeX = 10;
%!		sizeF = 10;
%!		vecX0 = zeros(sizeX,1);
%!		funchFJ = @(dummyX)( funcFJ_easy(dummyX) );
%!	otherwise
%!		error( "Invalid caseNum." );
%!	endswitch
%!	%
%!	[ vecXF, datOut ] = findLocMin_alytJ_mk2( vecX0, funchFJ );
