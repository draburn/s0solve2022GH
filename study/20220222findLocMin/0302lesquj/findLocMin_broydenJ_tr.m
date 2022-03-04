% Function...

function [ vecXF, datOut ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchF, prm=[] )
	%
	%
	% Parse input.
	sizeX = size(vecX0,1);
	sizeF = size(vecF0,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		%
		doTest = true;
		if (doTest)
			vecDelta_test = 0.001*ones(sizeX,1);
			vecJF_testA = ( funchF( vecX0 + vecDelta_test ) - funchF( vecX0 - vecDelta_test ) ) / 2.0;
			vecJF_testB = matJ0*vecDelta_test;
			if ( sumsq(vecJF_testA-vecJF_testB) > sqrt(eps)*(sumsq(vecJF_testA)+sumsq(vecJF_testA)) )
				msg( __FILE__, __LINE__, "WARNING: Input Jacobian appears inconsistent with input function." );
			endif
			clear vecDelta_test;
			clear vecJF_testA;
			clear vecJF_testB;
		endif
	endif
	%
	%
	%
	omegaTol = 0.0;
	gNormTol = 0.0;
	iterLimit = 1000;
	omegaFallAbsTol = eps^1.5;
	omegaFallRelTol = eps^1.5;
	deltaJRelTol = eps;
	trustRegionSize0 = -1.0;
	%
	coeff_reduceTrustRegionOnFevalFail = 0.1;
	coeff_declareModelIsRadicallyWrong = 3.0;
	coeff_reduceTrustRegionOnRadicallyWrong = 0.1;
	coeff_declareModeIsHighlyAccurate = 0.1;
	coeff_increaseTrustRegionOnHighlyAccurate = 1.5;
	trialLimit = 10;
	useCDL = true;
	useOmegaModelMin = true; % Not sure why, omegaModelMin really hurts in test case.
	updateTRBasedOnCDL = false;
	useLesquj = false;
	%
	if (~isempty(prm))
		useCDL = mygetfield( prm, "useCDL", useCDL );
		useLesquj = mygetfield( prm, "useLesquj", useLesquj );
	endif
	if (debugMode)
		assert( isscalar(useCDL) );
		assert( isbool(useCDL) );
		assert( isscalar(useLesquj) );
		assert( isbool(useLesquj) );
	endif
	%
	%
	fevalCount = 0;
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	trustRegionSize = trustRegionSize0;
	iterCount = 0;
	if (useLesquj)
		lesquj_vecXVals = vecX0;
		lesquj_vecFVals = vecF0;
		lesquj_prm = [];
		lesquj_prm.jevalDat(1).vecX = vecX0;
		lesquj_prm.jevalDat(1).vecF = vecF0;
		lesquj_prm.jevalDat(1).matJ = matJ0;
	endif
	if ( nargout >= 2 )
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.vecXVals(:,iterCount+1) = vecX;
		datOut.vecFVals(:,iterCount+1) = vecF;
		datOut.omegaVals(iterCount+1) = sumsq(vecF0)/2.0;
		datOut.matJVals(:,:,iterCount+1) = matJ;
	endif
	while (1)
		omega = sumsq(vecF)/2.0;
		vecG = matJ'*vecF;
		matH = matJ'*matJ;
		%
		%
		% Check pre-iter success conditions.
		if ( omega <= omegaTol )
			msgif( debugMode, __FILE__, __LINE__, "FULL SUCCESS: Reached omegaTol." );
			break;
		elseif ( norm(vecG) <= gNormTol )
			msgif( debugMode, __FILE__, __LINE__, "PARTIAL SUCCESS: Reached gNormTol." );
			break;
		endif
		%
		%
		% Check pre-iter imposed limits.
		iterCount++;
		%msgif( debugMode, __FILE__, __LINE__, sprintf( "~~~ ITERATION %d ~~~", iterCount ) );
		if ( iterCount > iterLimit )
			msg( __FILE__, __LINE__, "IMPOSED STOP: Reached iterLimit." );
			break;
		endif
		%
		%
		%
		trialCount = 0;
		while ( 1 )
			trialCount++;
			if ( trialCount > trialLimit )
				msg( __FILE__, __LINE__, "IMPOSED STOP: Reached trialLimit." );
				break;
			endif
			%
			%
			if (useCDL)
				cdlPrm = [];
				if ( 0.0 < trustRegionSize )
					cdlPrm.deltaNormMax = trustRegionSize;
				endif
				if ( useOmegaModelMin )
					% DRaburn 2022.03.03:
					%cdlPrm.omegaModelMin = -sqrt(eps)*omega; % This would make sense, but leads to worse behavior in xBroydenComp2.
					cdlPrm.omegaModelMin = 0.0;
				else
					cdlPrm.omegaModelMin = [];
				endif
				[ vecDeltaX, cdlDatOut ] = calcDeltaLev( omega, vecG, matH, cdlPrm );
				if ( updateTRBasedOnCDL )
				if ( cdlDatOut.trustRegionShouldBeUpdated )
					msgif( debugMode, __FILE__, __LINE__, "Updating trust region size per calcDeltaLev." );
					trustRegionSize = norm(vecDeltaX);
					% This seems to make things worse.
					% Maybe don't do this?
					% Just let the J update handle things?
				endif
				endif
			else
				[ matR, cholFlag ] = chol( matH );
				if ( 0 ~= cholFlag )
					hScale = max(abs(diag(matH)));
					matR = chol( matH + sqrt(eps)*hScale*eye(sizeX,sizeX) );
				endif
				vecDeltaX = -( matR \ (matR'\vecG) );
				if ( 0.0 < trustRegionSize )
					vecDeltaX *= trustRegionSize/norm(vecDeltaX);
				endif
			endif
			%
			vecFModel = vecF + matJ*vecDeltaX;
			%
			vecX_trial = vecX + vecDeltaX;
			vecF_trial = funchF( vecX_trial );
			fevalCount++;
			if (isempty(vecF_trial))
				msg( __FILE__, __LINE__, "Feval failed." );
				trustRegionSize = coeff_reduceTrustRegionOnFevalFail*norm(vecDeltaX);
				continue;
			endif
			%
			if (useLesquj)
				lesquj_vecXVals = [ lesquj_vecXVals, vecX_trial ];
				lesquj_vecFVals = [ lesquj_vecFVals, vecF_trial ];
			endif
			%
			if ( norm( vecF_trial - vecFModel ) >= coeff_declareModelIsRadicallyWrong * norm(vecF) )
				msgif( debugMode, __FILE__, __LINE__, "Model was radically wrong." );
				trustRegionSize = coeff_reduceTrustRegionOnRadicallyWrong*norm(vecDeltaX);
				continue;
			endif
			%
			omega_trial = sumsq(vecF_trial)/2.0;
			%
			% Note that omega_trial might not necessarily be a decrease,
			% but, the evaluated value of vecF is reasonable enough that we should at least update matJ.
			% But, if the model was good on the first try, increase the trust region size.
			if ( 1 == trialCount )
			if ( 0.0 < trustRegionSize )
			if ( norm( vecF_trial - vecFModel ) <= coeff_declareModeIsHighlyAccurate * norm(vecF) )
				msgif( debugMode, __FILE__, __LINE__, "Increasing trust region size." );
				trustRegionSize *= coeff_increaseTrustRegionOnHighlyAccurate;
			endif
			endif
			endif
			%
			break;
		endwhile
		%
		%
		if (useLesquj)
			[ lesquj_vecX0, lesquj_vecF0, lesquj_matJ0 ] = calcLesquj_basic( lesquj_vecXVals, lesquj_vecFVals, lesquj_prm );
			matJ_trial = lesquj_matJ0;
			matDeltaJ = matJ_trial - matJ;
			%
			%% HACK... Useless.
			%vecX = lesquj_vecX0;
			%vecF = lesquj_vecF0;
			%omega_trial = sumsq(vecF)/2.0;
		else
			vecDeltaY = vecF_trial - vecF - matJ*vecDeltaX;
			matDeltaJ = (vecDeltaY*(vecDeltaX'))/(vecDeltaX'*vecDeltaX);
			matJ_trial = matJ + matDeltaJ;
			if (debugMode)
				assert( reldiff(matJ_trial*vecDeltaX,vecF_trial-vecF) < sqrt(eps) );
			endif
		endif
		%
		% In this version, always accept matJ.
		% In a future version, use TR; but, still accept and matDeltaJ that is "reasonable".
		matJ = matJ_trial;
		if ( omega_trial < omega )
			vecX = vecX_trial;
			vecF = vecF_trial;
			% Should we stop?
			if ( abs(omega_trial-omega) <= omegaFallAbsTol )
				msgif( debugMode, __FILE__, __LINE__, "IMPOSED STOP: Reached omegaFallAbsTol." );
				break;
			elseif ( abs(omega_trial-omega) <= omegaFallRelTol*omega )
				msgif( debugMode, __FILE__, __LINE__, "IMPOSED STOP: Reached omegaFallRelTol." );
				break;
			endif
		else
			% We're rejecting the new point.
			% Only keep going if the change in matJ is adequate.
			if ( sqrt(sum(sumsq(matDeltaJ))) < deltaJRelTol*sqrt(sum(sumsq(matJ))) )
				msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Trial objective is not better and Jacobian change is small." );
				break;
			end
		endif
		if ( nargout >= 2 )
			datOut.deltaNormVals(iterCount) = norm(vecDeltaX);
			datOut.trustRegionSize(iterCount+1) = trustRegionSize;
			datOut.fevalCountVals(iterCount+1) = fevalCount;
			datOut.vecXVals(:,iterCount+1) = vecX;
			datOut.vecFVals(:,iterCount+1) = vecF;
			datOut.omegaVals(iterCount+1) = sumsq(vecF)/2.0;
			datOut.matJVals(:,:,iterCount+1) = matJ;
		endif
	endwhile
	%
	%
	vecXF = vecX;
	if ( nargout >= 2 )
		datOut.fevalCount = fevalCount;
		datOut.vecFF = vecF;
		datOut.matJF = matJ;
	endif
endfunction


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
%!	matJ = zeros(sizeX,sizeX);
%!	for n=1:sizeX
%!	for m=1:sizeX
%!		matJ(n,m) = 2.0*matA(n,m)*(vecX(m)-vecXE(m));
%!	end
%!	end
%!endfunction


%!test
%!	funchFJ = @(dummyX) funcFJ_nonlin(dummyX);
%!	vecX0 = zeros(2,1);
%!	[ vecF0, matJ0 ] = funchFJ( vecX0 );
%!	prm = [];
%!	[ vecXF, datOut ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm );
%!	echo__vecXF = vecXF
%!	omega0 = sumsq(funchFJ(vecX0),1)/2.0
%!	omegaF = sumsq(funchFJ(vecXF),1)/2.0
%!	%
%!	figure();
%!	semilogy( ...
%!	  datOut.fevalCountVals, datOut.omegaVals, 'o-', ...
%!	  datOut.fevalCountVals(2:end), datOut.deltaNormVals, 'x-' );
%!	grid on;
%!	xlabel( "feval count" );
%!	legend( "omega", "||delta||", "location", "northeast" );
%!	title( "findLocMin broydenJ tr" );
