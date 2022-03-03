% Function...

function [ vecXF, datOut ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchF, prm=[] )
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
	%
	useCDL = true;
	useOmegaModelMin = true; % Not sure why, omegaModelMin really hurts in test case.
	%
	if (~isempty(prm))
		useCDL = mygetfield( prm, "useCDL", useCDL );
	endif
	if (debugMode)
		assert( isscalar(useCDL) );
		assert( isbool(useCDL) );
	endif
	%
	fevalCount = 0;
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	iterCount = 0;
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
		if (useCDL)
			cdlPrm = [];
			if (useOmegaModelMin)
				% DRaburn 2022.03.03:
				%cdlPrm.omegaModelMin = -sqrt(eps)*omega; % This would make sense, but leads to worse behavior in xBroydenComp2.
				cdlPrm.omegaModelMin = 0.0;
			else
				cdlPrm.omegaModelMin = [];
			endif
			vecDeltaX = calcDeltaLev( omega, vecG, matH, cdlPrm );
		else
			[ matR, cholFlag ] = chol( matH );
			if ( 0 ~= cholFlag )
				hScale = max(abs(diag(matH)));
				matR = chol( matH + sqrt(eps)*hScale*eye(sizeX,sizeX) );
			endif
			vecDeltaX = -( matR \ (matR'\vecG) );
		endif
		vecX_trial = vecX + vecDeltaX;
		vecF_trial = funchF( vecX_trial );
		fevalCount++;
		omega_trial = sumsq(vecF_trial)/2.0;
		%
		vecDeltaY = vecF_trial - vecF - matJ*vecDeltaX;
		matDeltaJ = (vecDeltaY*(vecDeltaX'))/(vecDeltaX'*vecDeltaX);
		matJ_trial = matJ + matDeltaJ;
		if ( nargout >= 2 )
			datOut.deltaNormVals(iterCount) = norm(vecDeltaX);
			datOut.fevalCountVals(iterCount+1) = fevalCount;
			datOut.omegaVals(iterCount+1) = omega;
		endif
		if (debugMode)
			assert( reldiff(matJ_trial*vecDeltaX,vecF_trial-vecF) < sqrt(eps) );
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
%!	[ vecXF, datOut ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm );
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
%!	title( "findLocMin broydenJ condi" );
