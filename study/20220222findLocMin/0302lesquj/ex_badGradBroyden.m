% Function...

function [ vecXF, datOut ] = ex_badGradBroyden( vecX0, vecF0, matJ0, funchF, prm=[] )
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
	iterLimit = 100;
	omegaFallAbsTol = eps;
	omegaFallRelTol = eps;
	deltaJRelTol = eps;
	%
	fevalCount = 0;
	if ( nargout >= 2 )
		datOut.fevalCountVals(1) = fevalCount;
		datOut.omegaVals(1) = sumsq(vecF0)/2.0;
	endif
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	iterCount = 0;
	while (1)
		omega = sumsq(vecF)/2.0;
		vecG = matJ'*vecF;
		matH = matJ'*matJ;
		%
		doGradCheck = true;
		if (doGradCheck)
			vecF_check = funchF( vecX - 1.0e-4*vecG );
			echo__omega = omega
			omega_check = sumsq(vecF_check)/2.0
			vecF_anticheck = funchF( vecX + 1.0e-4*vecG );
			omega_anticheck = sumsq(vecF_anticheck)/2.0
			assert( omega_anticheck > omega )
		endif
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
		msgif( debugMode, __FILE__, __LINE__, sprintf( "~~~ ITERATION %d ~~~", iterCount ) );
		if ( iterCount > iterLimit )
			msg( __FILE__, __LINE__, "IMPOSED STOP: Reached iterLimit." );
			break;
		endif
		%
		vecDeltaX = calcDeltaLev( omega, vecG, matH );
		vecX_trial = vecX + vecDeltaX;
		vecF_trial = funchF( vecX_trial );
		fevalCount++;
		omega_trial = sumsq(vecF_trial)/2.0;
		%
		doBTHack = true;
		if (doBTHack)
			dTreg = norm(vecDeltaX);
			btLimit = 10;
			btCount = 0;
			while ( omega_trial >= omega )
				btCount++;
				assert( btCount <= btLimit );
				dTreg /= 2.0;
				cdlPrm = [];
				cdlPrm.dTreg = dTreg;
				vecDeltaX = calcDeltaLev( omega, vecG, matH, cdlPrm );
				vecX_trial = vecX + vecDeltaX;
				vecF_trial = funchF( vecX_trial );
				fevalCount++;
				omega_trial = sumsq(vecF_trial)/2.0;
			endwhile
		endif
		%
		vecDeltaY = vecF_trial - vecF - matJ*vecDeltaX;
		matDeltaJ = (vecDeltaY*(vecDeltaX'))/(vecDeltaX'*vecDeltaX);
		matJ_trial = matJ + matDeltaJ;
		if (debugMode)
			assert( reldiff(matJ_trial*vecDeltaX,vecF_trial-vecF) < sqrt(eps) );
		endif
		if ( nargout >= 2 )
			datOut.fevalCountVals(iterCount+1) = fevalCount;
			datOut.omegaVals(iterCount+1) = omega;
			datOut.deltaNormVals(iterCount) = norm(vecDeltaX);
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
%!	msg( __FILE__, __LINE__, "This test illustrates how matJ'*vecF can be in the wrong direction." );
%!	format long;
%!	funchFJ = @(dummyX) funcFJ_nonlin(dummyX);
%!	vecX0 = zeros(2,1);
%!	[ vecF0, matJ0 ] = funchFJ( vecX0 );
%!	prm = [];
%!	[ vecXF, datOut ] = ex_badGradBroyden( vecX0, vecF0, matJ0, funchFJ, prm );
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
