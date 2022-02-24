% Function...
%  [ vecX, datOut ] = findLocMin_alytJ( vecX0, funchFJ, prm=[] )

function [ vecX, datOut ] = findLocMin_alytJ( vecX0, funchFJ, prm=[] )
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
	%
	stepSizeTol = sqrt(eps);
	omegaTol = 0.0;
	gNormTol = 0.0;
	omegaFallAbsTol = eps;
	omegaFallRelTol = sqrt(eps);
	iterLimit = 100;
	stepType = 110;
	if ( ~isempty(prm) )
		stepSizeTol = mygetfield( prm, "stepSizeTol", stepSizeTol );
		omegaTol = mygetfield( prm, "omegaTol", omegaTol );
		gNormTol = mygetfield( prm, "gNormTol", gNormTol );
		omegaFallAbsTol = mygetfield( prm, "omegaFallAbsTol", omegaFallAbsTol );
		omegaFallRelTol = mygetfield( prm, "omegaFallRelTol", omegaFallRelTol );
		iterLimit = mygetfield( prm, "iterLimit", iterLimit );
		stepType = mygetfield( prm, "stepType", stepType );
	endif
	if (debugMode)
		assert( isrealscalar(stepSizeTol) );
		assert( 0.0 <= stepSizeTol );
		assert( isrealscalar(omegaTol) );
		assert( 0.0 <= omegaTol );
		assert( isrealscalar(gNormTol) );
		assert( 0.0 <= gNormTol );
		assert( isrealscalar(omegaFallAbsTol) );
		assert( 0.0 <= omegaFallAbsTol );
		assert( isrealscalar(omegaFallRelTol) );
		assert( 0.0 <= omegaFallRelTol );
		assert( omegaFallRelTol < 1.0 );
		assert( isrealscalar(iterLimit) );
		assert( isrealscalar(stepType) );
	endif
	%
	%
	% Set dummy return values.
	vecX = vecX0;
	if ( nargout >= 2 )
		datOut = [];
	endif
	fevalCount = 0;
	jevalCount = 0;
	%
	%
	% Do pre-loop analysis.
	[ vecF0, matJ0 ] = funchFJ( vecX0 );
	fevalCount++;
	jevalCount++;
	sizeF = size(vecF0,1);
	if (debugMode)
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		vecDelta_test = 0.001*ones(sizeX,1);
		vecJF_testA = ( funchFJ( vecX0 + vecDelta_test ) - funchFJ( vecX0 - vecDelta_test ) ) / 2.0;
		vecJF_testB = matJ0*vecDelta_test;
		if ( sumsq(vecJF_testA-vecJF_testB) > sqrt(eps)*(sumsq(vecJF_testA)+sumsq(vecJF_testA)) )
			error( "Jacobian calculated by funchFJ appears to be incorrect." );
		endif
		clear vecDelta_test;
		clear vecJF_testA;
		clear vecJF_testB;
	endif
	omega0 = sumsq(vecF0,1)/2.0;
	vecG0 = matJ0'*vecF0;
	matJTJ0 = matJ0'*matJ0;
	matK0 = zeros(sizeX,sizeX);
	%
	%
	% Handle "zero corner" cases.
	if ( omega0 <= omegaTol )
		msg( __FILE__, __LINE__, "Initial omega is below target min." );
		return;
	elseif ( norm(vecG0) <= gNormTol )
		msg( __FILE__, __LINE__, "Initial gradient is below tolerance." );
		return;
	endif
	%
	%
	%
	% Prep main loop.
	vecX = vecX0; % Probably repeated, but, clarifying.
	vecF = vecF0;
	matJ = matJ0;
	omega = omega0;
	vecG = vecG0;
	matJTJ = matJTJ0;
	matK = matK0;
	dTreg = [];
	iterCount = 0;
	datOut.fevalCountVals(iterCount+1) = fevalCount;
	datOut.omegaVals(iterCount+1) = omega;
	while (1)
		% Check pre-iter success conditions.
		if ( omega <= omegaTol )
			msgif( debugMode, __FILE__, __LINE__, "Reached omegaTol." );
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
			fevalCount += 2;
			matH = matJTJ + abs(h)*(vecPhiHat*(vecPhiHat'));
			vecX_next = vecX - matH \ vecG;
		case 21
			% Single enforced patch then Newton with enforced regularization.
			epsX = 1e-4;
			[ matPsi, matLambda ] = eig(matJTJ);
			[ lambdaAbsMin, nOfAbsMin ] = min(abs(diag(matLambda)));
			vecPhiHat = matPsi(:,nOfAbsMin);
			h = vecF'*( funchFJ( vecX + epsX*vecPhiHat ) + funchFJ( vecX - epsX*vecPhiHat ) - 2.0*vecF )/(epsX*epsX);
			fevalCount += 2;
			matH = matJTJ + abs(h)*(vecPhiHat*(vecPhiHat'));
			hScale = max(abs(diag(matH)));
			vecX_next = vecX - ( matH + 1e-4*hScale*eye(sizeX) ) \ vecG;
		case 100
			% Baisc flmcj.
			flmcjPrm = [];
			[ vecX_next, flmcjDatOut ] = findLocMin_cnstJ( vecX, vecF, matJ, funchFJ, flmcjPrm );
			fevalCount += flmcjDatOut.fevalCount;
		case 110
			% flmcj with K and dTreg passing.
			flmcjPrm = [];
			flmcjPrm.matK = matK;
			flmcjPrm.deltaNormMax = dTreg; % Could be empgy;
			[ vecX_next, flmcjDatOut ] = findLocMin_cnstJ( vecX, vecF, matJ, funchFJ, flmcjPrm );
			dTreg = flmcjDatOut.deltaNormMax;
			matK = flmcjDatOut.matK;
			fevalCount += flmcjDatOut.fevalCount;
		otherwise
			error( "Invalid stepType." );
		endswitch		
		%
		%
		% Validate step.
		if ( ~isrealarray(vecX_next,[sizeX,1]) );
			msg( __FILE__, __LINE__, "vecX_next is not a valid vector." );
			return;
		elseif ( norm(vecX_next-vecX) <= stepSizeTol )
			msg( __FILE__, __LINE__, "Step is below tol." );
			return;
		elseif ( norm(vecX_next-vecX) > dTreg )
			msg( __FILE__, __LINE__, "Step moves outside trust region!" );
		endif
		%echo__vecX_next = vecX_next
		datOut.deltaNormVals(iterCount) = norm( vecX_next - vecX );
		[ vecF_next, matJ_next ] = funchFJ( vecX_next );
		fevalCount++;
		jevalCount++;
		omega_next = sumsq(vecF_next)/2.0;
		if ( omega_next >= omega )
			msg( __FILE__, __LINE__, "Omega did not decrease." );
			datOut.fevalCountVals(iterCount+1) = fevalCount;
			datOut.omegaVals(iterCount+1) = omega_next;
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
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.omegaVals(iterCount+1) = omega;
		%
		%
		% Check post-iter imposed limits.
		if ( abs(omega-omega_prev) <= omegaFallAbsTol )
			msg( __FILE__, __LINE__, "Reached omegaFallAbsTol." );
			return;
		elseif ( abs(omega-omega_prev) <= omegaFallRelTol*abs(omega) )
			msg( __FILE__, __LINE__, "Reached omegaFallRelTol." );
			return;
		endif
	endwhile
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
%!	prm = [];
%!	[ vecXF, datOut ] = findLocMin_alytJ( vecX0, funchFJ, prm );
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
