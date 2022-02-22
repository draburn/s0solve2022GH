% Function ...
%  2022-02-18 Let's study a new approach to "backtracking".
%   But, don't worry about re-using information between "nonlinear iterations".
%   And, assume we have convenient analytic J.
function [ vecX, datOut ] = findLocMin_alytJ_mk3( vecX0, funchFJ, prm=[] )
	% Use retcode?
	%
	%
	% Parse input.
	sizeX = size(vecX0,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
	endif
	%
	omegaMin = mygetfield( prm, "omegaMin", 0.0 );
	gNormTol = mygetfield( prm, "gNormTol", 0.0 );
	iterLimit = mygetfield( prm, "iterLimit", 5 );
	%stepType = mygetfield( prm, "stepType", 102 );
	%stepType = mygetfield( prm, "stepType", 200 );
	stepType = mygetfield( prm, "stepType", 201 );
	c0_trustRegion = mygetfield( prm, "c0_trustRegion", 100.0 );
	stepSizeTol = mygetfield( prm, "stepSizeTol", sqrt(eps) );
	omegaFallAbsTol = mygetfield( prm, "omegaFallAbsTol", eps );
	omegaFallRelTol = mygetfield( prm, "omegaFallRelTol", sqrt(eps) );
	if (debugMode)
		assert( isrealscalar(omegaMin) );
		assert( isrealscalar(gNormTol) );
		assert( isrealscalar(iterLimit) );
		assert( isrealscalar(c0_trustRegion) );
		assert( isrealscalar(stepSizeTol) );
		assert( isrealscalar(omegaFallAbsTol) );
		assert( isrealscalar(omegaFallRelTol) );
	endif
	%
	%
	% Set dummy return values.
	vecX = vecX0;
	if ( nargout >= 2 )
		datOut = [];
	endif
	%
	%
	% Do pre-loop analysis.
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
	omega0 = sumsq(vecF0,1)/2.0;
	if ( omega0 <= omegaMin )
		msg( __FILE__, __LINE__, "Initial omega is below target min." );
		return;
	endif
	vecG0 = matJ0'*vecF0;
	if ( norm(vecG0) <= gNormTol )
		msg( __FILE__, __LINE__, "Initial gradient is below tolerance." );
		return;
	endif
	matJTJ0 = matJ0'*matJ0;
	dTreg0 = c0_trustRegion * norm(vecF0) / sqrt(max(abs(diag(matJTJ0))));
	%%%dTreg0 = c0_trustRegion * norm(vecG0) / max(abs(diag(matJTJ0))); % This might make more sense.
	%
	trackHJBroyd = false;
	if ( trackHJBroyd )
		matJ_HJBroyd = matJ0
		matH_HJBroyd = matJ0'*matJ0
	endif
	%
	winters_dTreg = [];
	winters_matK = [];
	%
	% Prep main loop.
	vecX = vecX0; % Probably repeated, but, clarifying.
	vecF = vecF0;
	matJ = matJ0;
	omega = omega0;
	vecG = vecG0;
	matJTJ = matJTJ0;
	dTreg = dTreg0;
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
		msg( __FILE__, __LINE__, sprintf( "~~~ ITERATION %d ~~~", iterCount ) );
		if ( iterCount > iterLimit )
			msgif( debugMode, __FILE__, __LINE__, "Reached iterLimit." );
			return;
		endif
		%
		%
		% Do work!
		% Note that dTreg may be change internally.
		switch (stepType)
		case 0
			% Simple grad.
			vecX_next = vecX - 0.0001*vecG;
		case 10
			% Blind Newton.
			vecX_next = vecX - matJ\vecG;
		case 11
			% Blind Newton with enforced regularization.
			hScale = max(abs(diag(matJTJ)));
			vecX_next = vecX - ( matJTJ + 1e-4*hScale*eye(sizeX) ) \ vecG;
		case 20
			% Single enforced patch then blind Newton.
			epsX = 1e-4;
			[ matPsi, matLambda ] = eig(matJTJ);
			[ lambdaAbsMin, nOfAbsMin ] = min(abs(diag(matLambda)));
			vecPhiHat = matPsi(:,nOfAbsMin);
			h = vecF'*( funchFJ( vecX + epsX*vecPhiHat ) + funchFJ( vecX - epsX*vecPhiHat ) - 2.0*vecF )/(epsX*epsX);
			matH = matJTJ + abs(h)*(vecPhiHat*(vecPhiHat'));
			vecX_next = vecX - matH \ vecG;
		case 21
			% Single enforced patch then Newton with enforced regularization.
			epsX = 1e-4;
			[ matPsi, matLambda ] = eig(matJTJ);
			[ lambdaAbsMin, nOfAbsMin ] = min(abs(diag(matLambda)));
			vecPhiHat = matPsi(:,nOfAbsMin);
			h = vecF'*( funchFJ( vecX + epsX*vecPhiHat ) + funchFJ( vecX - epsX*vecPhiHat ) - 2.0*vecF )/(epsX*epsX);
			matH = matJTJ + abs(h)*(vecPhiHat*(vecPhiHat'));
			hScale = max(abs(diag(matH)));
			vecX_next = vecX - ( matH + 1e-4*hScale*eye(sizeX) ) \ vecG;
		case 100
			findLocMin_alytJ_mk3__findNext_hupd;
		case 102
			fnprm.dTreg = dTreg;
			vecX_next = findLocMin_alytJ_mk3__findNext_mk2( vecX, vecF, matJ, funchFJ, fnprm );
		case 200
			fnprm = [];
			vecX_next = findLocMin_alytJ_mk3__findNext_winters( vecX, vecF, matJ, funchFJ, fnprm );
		case 201
			fnprm = [];
			if ( isrealscalar(winters_dTreg) )
				fnprm.dTreg0 = winters_dTreg;
			endif
			if ( isrealarray(winters_matK,[sizeX,sizeX]) )
				fnprm.matK0 = winters_matK;
			endif
			[ vecX_next, fndatOut ] = findLocMin_alytJ_mk3__findNext_winters( vecX, vecF, matJ, funchFJ, fnprm );
			winters_matK = fndatOut.matK;
			winters_dTreg = fndatOut.dTreg;
		otherwise
			error( "Invalid stepType." );
		endswitch
		%
		%
		% Validate step.
		if ( ~isrealarray(vecX_next,[sizeX,1]) );
			msg( __FILE__, __LINE__, "vecX_next is not a valid vector." );
			return;
		elseif ( norm(vecX_next-vecX) > dTreg )
			msg( __FILE__, __LINE__, "Step moves outside trust region!" );
		endif
		%echo__vecX_next = vecX_next
		[ vecF_next, matJ_next ] = funchFJ( vecX_next );
		omega_next = sumsq(vecF_next)/2.0
		if ( omega_next >= omega )
			msg( __FILE__, __LINE__, "Omega did not decrease." );
			echo__omega_next = omega_next;
			return;
		endif
		vecG_next = matJ_next'*vecF_next;
		matJTJ_next = matJ_next'*matJ_next;
		%
		%
		% Take step.
		vecX_prev = vecX;
		vecF_prev = vecF;
		matJ_prev = matJ;
		omega_prev = omega;
		vecG_prev = vecG;
		matJTJ_prev = matJTJ;
		%
		vecX = vecX_next;
		vecF = vecF_next;
		matJ = matJ_next;
		omega = omega_next;
		vecG = vecG_next;
		matJTJ = matJTJ_next;
		%
		%
		useWintersK = false;
		if ( useWintersK )
			if ( isrealarray(winters_matK,[sizeX,sizeX]) )
				% Impose deltaK * deltaX = (deltaJ)'*J*deltaX.
				fooX = vecX - vecX_prev;
				fooY = (matJ-matJ_prev)' * (matJ*(fooX));
				fooS = sumsq(fooX);
				winters_matK += ( fooY*(fooX') + fooX*(fooY') - (fooX*(fooX'))*(fooX'*fooY)/fooS )/fooS;
			else
				msg( __FILE__, __LINE__, "Don't have winters K." );
			endif
		endif
		%
		%
		% Check post-iter imposed limits.
		if ( norm(vecX-vecX_prev) <= stepSizeTol )
			msg( __FILE__, __LINE__, "Reached stepSizeTol." );
			return;
		elseif ( abs(omega-omega_prev) <= omegaFallAbsTol )
			msg( __FILE__, __LINE__, "Reached omegaFallAbsTol." );
			return;
		elseif ( abs(omega-omega_prev) <= omegaFallRelTol*abs(omega) )
			msg( __FILE__, __LINE__, "Reached omegaFallRelTol." );
			return;
		endif
	endwhile
return;
endfunction


%!function [ vecF, matJ ] = funcFJ_trivial( vecX )
%!	sizeX = size(vecX,1);
%!	vecF = vecX-1.0;
%!	matJ = eye(sizeX,sizeX);
%!endfunction
%!function [ vecF, matJ ] = funcFJ_easy( vecX )
%!	sizeX = size(vecX,1);
%!	matJ = diag((1:sizeX));
%!	vecF = matJ*( vecX - (1:sizeX)' );
%!endfunction
%!function [ vecF, matJ ] = funcFJ_nonlin( vecX )
%!	sizeX = size(vecX,1);
%!	vecXE = (1:sizeX)';
%!	matA = diag((1:sizeX));
%!	matA(2,1) = (sizeX+1);
%!	vecF = matA*((vecX-vecXE)).^2;
%!	matJ = diag((1:sizeX));
%!endfunction


%!test
%!	caseNum = 104;
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
%!	case 10
%!		sizeX = 2;
%!		sizeF = 2;
%!		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0); % Calls setprngstates.
%!		%echo__vecFE = testFuncPrm.vecFE
%!		%echo__omegaE = sumsq(testFuncPrm.vecFE)/2.0
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = zeros(sizeX,1);
%!	case 11
%!		sizeX = 2;
%!		sizeF = 2;
%!		tfpPrm.matJPreMod = [ 1.0, 1.0; 0.0, 0.0 ];
%!		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0,true,true,true,tfpPrm); % Calls setprngstates.
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = zeros(sizeX,1);
%!	case 20
%!		sizeX = 2;
%!		sizeF = 2;
%!		tfpPrm.matJPreMod = ones(sizeF,sizeX);
%!		tfpPrm.matJPreMod(1,1) = 100.0;
%!		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0,true,true,true,tfpPrm); % Calls setprngstates.
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = zeros(sizeX,1);
%!	case 100
%!		msg( __FILE__, __LINE__, "*** WARNING: This is a 'perfectly balanced' case! ***" );
%!		sizeX = 2;
%!		sizeF = 2;
%!		testFuncPrm.sizeX = 2;
%!		testFuncPrm.sizeF = 2;
%!		testFuncPrm.vecXE = [ 0.0; 0.0 ];
%!		testFuncPrm.vecFE = [ 0.0; 1.0 ];
%!		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = ones(sizeX,1);
%!	case 101
%!		msg( __FILE__, __LINE__, "*** WARNING: This is a 'perfectly balanced' case! ***" );
%!		sizeX = 2;
%!		sizeF = 2;
%!		testFuncPrm.sizeX = 2;
%!		testFuncPrm.sizeF = 2;
%!		testFuncPrm.vecXE = [ 0.0; 0.0 ];
%!		testFuncPrm.vecFE = [ 0.0; 1.0 ];
%!		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = 0.01*ones(sizeX,1);
%!	case 102
%!		sizeX = 2;
%!		sizeF = 2;
%!		testFuncPrm.sizeX = 2;
%!		testFuncPrm.sizeF = 2;
%!		testFuncPrm.vecXE = [ 0.0; 0.0 ];
%!		testFuncPrm.vecFE = [ 0.0; 1.0 ];
%!		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = [ 1.0; 0.9 ];
%!	case 103
%!		sizeX = 2;
%!		sizeF = 2;
%!		testFuncPrm.sizeX = 2;
%!		testFuncPrm.sizeF = 2;
%!		testFuncPrm.vecXE = [ 0.0; 0.0 ];
%!		testFuncPrm.vecFE = [ 0.0; 1.0 ];
%!		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = 0.01*[ 1.0; 0.9 ];
%!	case 104
%!		sizeX = 2;
%!		sizeF = 2;
%!		testFuncPrm.sizeX = 2;
%!		testFuncPrm.sizeF = 2;
%!		testFuncPrm.vecXE = [ 0.0; 0.0 ];
%!		testFuncPrm.vecFE = [ 0.0; 1.0 ];
%!		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
%!		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
%!		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
%!		vecX0 = [ 1.0; 2.0 ];
%!	otherwise
%!		error( "Invalid caseNum." );
%!	endswitch
%!	%
%!	[ vecXF, datOut ] = findLocMin_alytJ_mk3( vecX0, funchFJ );
