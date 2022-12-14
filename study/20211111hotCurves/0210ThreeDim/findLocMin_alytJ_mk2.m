% Function ...
function [ vecX, datOut ] = findLocMin_alytJ_mk2( vecX0, funchFJ, prm=[] )
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
	iterLimit = mygetfield( prm, "iterLimit", 10 );
	omegaFallTol = mygetfield( prm, "omegaFallTol", sqrt(eps) );
	stepSizeTol = mygetfield( prm, "stepSizeTol", sqrt(eps) );
	%patchDimMax = mygetfield( prm, "patchDimMax", 0 ); patchThresh = mygetfield( prm, "patchThresh", 0.0 );
	%patchDimMax = mygetfield( prm, "patchDimMax", 1 ); patchThresh = mygetfield( prm, "patchThresh", 1e-2 );
	patchDimMax = mygetfield( prm, "patchDimMax", 1 ); patchThresh = mygetfield( prm, "patchThresh", 2.0 );
	%patchDimMax = mygetfield( prm, "patchDimMax", sizeX ); patchThresh = mygetfield( prm, "patchThresh", 2.0 );
	patchEpsX = mygetfield( prm, "patchEpsX", 1e-5 );
	relFallThresh = mygetfield( prm, "relFallThresh", 0.5 );
	if (debugMode)
		assert( isrealscalar(omegaMin) );
		assert( isrealscalar(gNormTol) );
		assert( isrealscalar(iterLimit) );
		assert( isrealscalar(omegaFallTol) );
		assert( isrealscalar(stepSizeTol) );
		assert( isrealscalar(patchDimMax) );
		assert( isrealscalar(patchThresh) );
		assert( isrealscalar(patchEpsX) );
		assert( isrealscalar(relFallThresh) );
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
	endif
	%
	%
	trackHJBroyd = false;
	if ( trackHJBroyd )
		matJ_HJBroyd = matJ0
		matH_HJBroyd = matJ0'*matJ0
	endif
	%
	% Prep main loop.
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	omega = omega0
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
		%msg( __FILE__, __LINE__, sprintf( "~~~ ITERATION %d ~~~", iterCount ) );
		if ( iterCount > iterLimit )
			msgif( debugMode, __FILE__, __LINE__, "Reached iterLimit." );
			return;
		endif
		%
		%
		% Calculate Hessian and Patch.
		matJTJ = matJ'*matJ;
		%
		matP = zeros(sizeX,sizeX);
		if ( 1 <= patchDimMax )
			[ matPsi, matLambda ] = eig( matJTJ );
			[ vecAbsLambdaSortVal, vecAbsLambdaSortIndex ] = sort(abs(diag(matLambda)));
			vecEigToPatch = vecAbsLambdaSortIndex( vecAbsLambdaSortVal <= patchThresh*vecAbsLambdaSortVal(end) );
			if ( length(vecEigToPatch) > patchDimMax )
				vecEigToPatch = vecEigToPatch(1:patchDimMax);
			endif
			%
			for n=1:length(vecEigToPatch);
				m = vecEigToPatch(n);
				vecFP = funchFJ( vecX + patchEpsX*matPsi(:,m) );
				vecFM = funchFJ( vecX - patchEpsX*matPsi(:,m) );
				assert( abs( norm(matPsi(:,m)) - 1.0 ) < sqrt(eps) );
				h = (vecF'*( vecFP + vecFM - (2.0*vecF) )) / (patchEpsX*patchEpsX);
				if ( 0.0 >= h )
					msg( __FILE__, __LINE__, "Would-be h term is not positive. This scenario deserves futher attention!" );
				endif
				matP += abs(h)*(matPsi(:,m)*(matPsi(:,m)'));
			endfor
			%
			if (debugMode&&(length(vecEigToPatch)>0))
				funchOmegaFD = @(dummyX)( sumsq(funchFJ(dummyX))/2.0 );
				[ omegaFD, vecGFD, matHFD ] = evalFDGH( vecX, funchOmegaFD );
				patchCheckA = diag(matPsi'*matHFD*matPsi)(vecEigToPatch);
				patchCheckB = diag(matPsi'*(matJTJ+matP)*matPsi)(vecEigToPatch);
				if ( reldiff(patchCheckA,patchCheckB) > eps^0.2 )
					echo__patchCheckA = patchCheckA'
					echo__patchCheckB = patchCheckB'
					msg( __FILE__, __LINE__, "*** WARNING: patchCheckA and patchCheckB do not seem to agree. This should be impossible. ***" );
				endif
			endif
		endif
		%rcond(matJTJ)
		matH = matJTJ + matP;
		[ omegaFD, vecGFD, matHFD ] = evalFDGH( vecX, @(dummyX)( sumsq(funchFJ(dummyX))/2.0 ) );
		echo__matHFD = matHFD;
		%%%matH = matHFD;
		%rcond(matH)
		%msg( __FILE__, __LINE__, "RETURN!" ); return
		if ( 1==iterCount && trackHJBroyd )
			matH_HJBroyd += matP
		endif
		%
		%
		%
		% Pick step.
		% THIS IS VERY CRUDE.
		hScale = eigs(matH,1);
		assert( hScale > 0.0 );
		for mu = hScale*[ 0.0, 1e-8, 1e-4, 1e-2, 1e-1, 1e0, 2e0, 3e0, 1e1, 1e2, 1e3, 1e4 ]
		%for mu = 0
			[ matR, cholFlag ] = chol( matH + (mu*matIX) );
			if ( 0 ~= cholFlag )
				msg( __FILE__, __LINE__, "Cholesky factorization of Hessian failed." );
				continue;
			endif
			vecDeltaX_trial = -( matR \ ( matR' \ vecG ) );
			%
			vecX_trial = vecX + vecDeltaX_trial;
			omegaModel_trial = omega + vecG'*vecDeltaX_trial + 0.5*(vecDeltaX_trial'*matH*vecDeltaX_trial);
			omegaFallModel_trial = omega - omegaModel_trial;
			assert( 0.0 < omegaFallModel_trial );
			vecF_trial = funchFJ( vecX_trial );
			omega_trial = sumsq(vecF_trial)/2.0;
			omegaFall_trial = omega - omega_trial;
			if ( omegaFall_trial < relFallThresh*omegaFallModel_trial )
				% Fall is insufficient, possibly negative.
				% NOTE: As long as vecF isn't too crazy, we might be able to pick up some info about H.
				% This is for later.
				% IF OMEGAFALL_TRIAL INCREASES FROM PREV MU, USE THAT MU INSTEAD?
				continue;
			endif
			mu_final = mu
			break;
		endfor
		%
		%
		% Validate step.
		vecX_next = vecX_trial;
		[ vecF_next, matJ_next ] = funchFJ( vecX_next );
		omega_next = sumsq(vecF_next)/2.0;
		if ( omega_next >= omega )
			msg( __FILE__, __LINE__, "Omega increased." );
			echo__omega_next = omega_next
			return;
		endif
		vecG_next = matJ_next'*vecF_next;
		%
		if ( trackHJBroyd )
			vecDeltaX = vecX_next - vecX;
			vecDeltaF = vecF_next - vecF;
			matDeltaJ_HJBroyd = ( vecDeltaF - matJ_HJBroyd * vecDeltaX ) * (vecDeltaX') / sumsq(vecDeltaX);
			matJ_HJBroyd_next = matJ_HJBroyd + matDeltaJ_HJBroyd;
			assert( reldiff(matJ_HJBroyd_next*vecDeltaX,vecDeltaF) < sqrt(eps) );
			%
			vecDeltaG = matJ_HJBroyd_next'*vecF_next - matJ_HJBroyd'*vecF; % This eq is particularly suspect.
			foo = matH_HJBroyd*vecDeltaX;
			matDeltaH_HJBroyd = (vecDeltaG*(vecDeltaG'))/(vecDeltaG'*vecDeltaX) - (foo*(foo'))/(foo'*vecDeltaX);
			matH_HJBroyd_next = matH_HJBroyd + matDeltaH_HJBroyd;
			assert( reldiff(matH_HJBroyd_next*vecDeltaX,vecDeltaG) < sqrt(eps) );
			%
			matJ_HJBroyd = matJ_HJBroyd_next
			matH_HJBroyd = matH_HJBroyd_next
		endif
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
		omega = omega_next
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
%!	caseNum = 100;
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
%!	otherwise
%!		error( "Invalid caseNum." );
%!	endswitch
%!	%
%!	[ vecXF, datOut ] = findLocMin_alytJ_mk2( vecX0, funchFJ );
