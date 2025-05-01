% Function ...
function [ vecX, datOut ] = findLocMin_alytJ( vecX0, funchFJ, prm=[] )
	% Use retcode?
	%
	sizeX = size(vecX0,1);
	%
	debugMode = true;
	omegaMin_apch = 1e-1;
	omegaMin_cnvg = 0.0;
	iterLim_apch = 100;
	iterLim_cnvg = 3;
	gNormSqTol_apch = 0.0; %!!!
	gNormSqTol_cnvg = 0.0;
	stepType_apch = 1;
	stepType_cnvg = 10;
	
	%stepSizeTol
	%omegaFallTol
	%stepsizeLim
	%modelAccuracyThresh
	
	if (~isempty(prm))
		debugMode = mygetfield( prm, "debugMode", debugMode );
		omegaMin_apch = mygetfield( prm, "omegaMin_apch", omegaMin_apch );
		omegaMin_cnvg = mygetfield( prm, "omegaMin_cnvg", omegaMin_cnvg );
		iterLim_apch = mygtfield( prm, "iterLim_apch", iterLim_apch );
		iterLim_cnvg = mygtfield( prm, "iterLim_cnvg", iterLim_cnvg );
		gNormSqTol_apch = mygetfield( prm, "gNormSqTol_apch", gNormSqTol_apch );
		gNormSqTol_cnvg = mygetfield( prm, "gNormSqTol_cnvg", gNormSqTol_cnvg );
		stepType_apch = mygetfield( prm, "stepType_apch", stepType_apch );
		stepType_cnvg = mygetfield( prm, "stepType_cnvg", stepType_cnvg );
	endif
	%
	if (debugMode)
		msg( __FILE__, __LINE__, "Using debugMode." );
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealscalar(omegaMin_apch) );
		assert( isrealscalar(omegaMin_cnvg) );
		assert( isrealscalar(iterLim_apch) );
		assert( isrealscalar(iterLim_cnvg) );
		assert( isrealscalar(gNormSqTol_apch) );
		assert( isrealscalar(gNormSqTol_cnvg) );
		assert( isrealscalar(stepType_apch) );
		assert( isrealscalar(stepType_cnvg) );
	endif
	%
	[ vecF0, matJ0 ] = funchFJ( vecX0 );
	sizeF = size(vecF0,1);
	if (debugMode)
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		if ( 0.0 ~= norm(vecF0) )
			vecDelta_test = 0.001*ones(sizeX,1);
			vecJF_testA = ( funchFJ( vecX0 + vecDelta_test ) - funchFJ( vecX0 - vecDelta_test ) ) / 2.0;
			vecJF_testB = matJ0*vecDelta_test;
			if ( sumsq(vecJF_testA-vecJF_testB) > sqrt(eps)*(sumsq(vecJF_testA)+sumsq(vecJF_testA)) )
				msg( __FILE__, __LINE__, "*** WARNING: Jacobian calculated by funchFJ appears to be incorrect. ***" );
			end
			clear vecDelta_test;
			clear vecJF_testA;
			clear vecJF_testB;
		end
	endif
	%
	omega0 = sumsq(vecF0,1)/2.0;
	vecG0 = matJ0'*vecF0;
	%
	%
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	omega = omega0;
	vecG = vecG0;
	gNormSq = sumsq(vecG);
	if ( nargout >= 2 )
		datOut.iterCount_apch = 0;
		datOut.iterCount_cnvg = 0;
		datOut.iterCount_tot = 0;
		datOut.vecXVals(:,datOut.iterCount_tot+1) = vecX;
		datOut.vecFVals(:,datOut.iterCount_tot+1) = vecF;
		datOut.omegaVals(datOut.iterCount_tot+1) = omega;
		datOut.vecGVals(:,datOut.iterCount_tot+1) = vecG;
		datOut.gNormSqVals(datOut.iterCount_tot+1) = gNormSq;
	endif
	msgif( debugMode, __FILE__, __LINE__, sprintf( "Start:  bigDeltaNormSq = %0.3e, omega = %0.3e,  gNormSq = %0.3e.", sumsq(vecX-vecX0), omega, gNormSq ) );
	%
	%
	%
	iterCount_apch = 0;
	while (1)
		% Check pre-iter success conditions.
		if ( omega <= omegaMin_apch )
			msgif( debugMode, __FILE__, __LINE__, "Reached omegaMin_apch." );
			approachWasSuccessful = true;
			break;
		endif
		%
		if ( gNormSq <= gNormSqTol_apch )
			msgif( debugMode, __FILE__, __LINE__, "Reached gNormSqTol_apch." );
			approachWasSuccessful = true;
			break;
		end
		assert( gNormSq ~= 0.0 );
		%
		%
		% Check pre-iter imposed limits.
		iterCount_apch++;
		if ( iterCount_apch > iterLim_apch )
			msgif( debugMode, __FILE__, __LINE__, "Reached iterLim_apch." );
			approachWasSuccessful = false;
			break;
		endif
		%
		%
		% Do iteration.
		switch (stepType_apch)
		case 0
			% "Gradient-descent line segment" step, to minimum along gradient-descent (with H = J' * J), no BT.
			% Note that, H = J' * J is inexact but yields...
			%  g' * H * g = (J*g)' * (J*g);
			%matH = matJ'*matJ;
			%p = calcLinishRootOfQuad( 0.5*(vecG'*matH*vecG), -(vecG'*vecG), omega )
			p = calcLinishRootOfQuad( sumsq(matJ*vecG)/2.0, -sumsq(vecG), omega );
			if ( p < 0.0 )
				error( "calcLinishRootOfQuad() returned a non-positive value; this should be impossible." );
			end
			vecDelta = -p*vecG;
		case 1
			% Scaled Gradient-descent line segment step.
			vecScl = sum(matJ.*matJ,1);
			matS = diag( vecScl + sqrt(eps)*max(vecScl) );
			vecGS = matS\vecG;
			p = calcLinishRootOfQuad( sumsq(matJ*vecGS)/2.0, -(vecG'*vecGS), omega );
			if ( p < 0.0 )
				error( "calcLinishRootOfQuad() returned a non-positive value; this should be impossible." );
			end
			vecDelta = -p*vecGS;
		case 10
			% Full Newton step, no BT.
			[ matR, cholFlag ] = chol( matJ'*matJ );
			if ( 0~=cholFlag )
				msg( __FILE__, __LINE__, "Cholesky factorization of Hessian failed." );
				return;
			end
			vecDelta = -( matR \ ( matR' \ vecG ) );
		otherwise
			error( "Invalid apchStep_type." );
			return;
		endswitch
		%
		vecX_next = vecX + vecDelta;
		[ vecF_next, matJ_next ] = funchFJ( vecX_next );
		omega_next = sumsq(vecF_next)/2.0;
		if ( omega_next >= omega )
			msg( __FILE__, __LINE__, "Omega increased." );
			approachWasSuccessful = false;
			break;
		end
		%
		vecX = vecX_next;
		vecF = vecF_next;
		matJ = matJ_next;
		omega = omega_next;
		vecG = matJ'*vecF;
		gNormSq = sumsq(vecG);
		if ( nargout >= 2 )
			datOut.iterCount_tot++;
			datOut.vecXVals(:,datOut.iterCount_tot+1) = vecX;
			datOut.vecFVals(:,datOut.iterCount_tot+1) = vecF;
			datOut.omegaVals(datOut.iterCount_tot+1) = omega;
			datOut.vecGVals(:,datOut.iterCount_tot+1) = vecG;
			datOut.gNormSqVals(datOut.iterCount_tot+1) = gNormSq;
			datOut.iterCount_apch = iterCount_apch;
		endif
	endwhile
	msgif( debugMode, __FILE__, __LINE__, sprintf( 
	  "Approach:  iterCount = %d,  bigDeltaNormSq = %0.3e,  omega = %0.3e,  gNormSq = %0.3e.",
	  iterCount_apch, sumsq(vecX-vecX0), omega, gNormSq ) );
	if ( ~approachWasSuccessful )
		msgif( debugMode, __FILE__, __LINE__, "Approach failed." );
		return;
	end
	%
	%
	%
	iterCount_cnvg = 0;
	while (1)
		% Check pre-iter success conditions.
		if ( omega <= omegaMin_cnvg )
			msgif( debugMode, __FILE__, __LINE__, "Reached omegaMin_cnvg." );
			cnvgWasSuccessful = true;
			break;
		endif
		%
		if ( gNormSq <= gNormSqTol_cnvg )
			msgif( debugMode, __FILE__, __LINE__, "Reached gNormSqTol_cnvg." );
			cnvgWasSuccessful = true;
			break;
		end
		assert( gNormSq ~= 0.0 );
		%
		%
		% Check pre-iter imposed limits.
		iterCount_cnvg++;
		if ( iterCount_cnvg > iterLim_cnvg )
			msgif( debugMode, __FILE__, __LINE__, "Reached iterLim_cnvg." );
			cnvgWasSuccessful = false;
			break;
		endif
		%
		%
		% Do iteration.
		switch (stepType_cnvg)
		case 0
			% "Gradient-descent line segment" step, to minimum along gradient-descent (with H = J' * J), no BT.
			% Note that, H = J' * J is inexact but yields...
			%  g' * H * g = (J*g)' * (J*g);
			%matH = matJ'*matJ;
			%p = calcLinishRootOfQuad( 0.5*(vecG'*matH*vecG), -(vecG'*vecG), omega )
			p = calcLinishRootOfQuad( sumsq(matJ*vecG)/2.0, -sumsq(vecG), omega );
			if ( p < 0.0 )
				error( "calcLinishRootOfQuad() returned a non-positive value; this should be impossible." );
			end
			vecDelta = -p*vecG;
		case 1
			% Scaled Gradient-descent line segment step.
			vecScl = sum(matJ.*matJ,1);
			matS = diag( vecScl + sqrt(eps)*max(vecScl) );
			vecGS = matS\vecG;
			p = calcLinishRootOfQuad( sumsq(matJ*vecGS)/2.0, -(vecG'*vecGS), omega );
			if ( p < 0.0 )
				error( "calcLinishRootOfQuad() returned a non-positive value; this should be impossible." );
			end
			vecDelta = -p*vecGS;
		case 10
			% Full Newton step, no BT.
			[ matR, cholFlag ] = chol( matJ'*matJ );
			if ( 0~=cholFlag )
				msg( __FILE__, __LINE__, "Cholesky factorization of Hessian failed." );
				return;
			end
			vecDelta = -( matR \ ( matR' \ vecG ) );
		otherwise
			error( "Invalid cnvgStep_type." );
			return;
		endswitch
		%
		vecX_next = vecX + vecDelta;
		[ vecF_next, matJ_next ] = funchFJ( vecX_next );
		omega_next = sumsq(vecF_next)/2.0;
		if ( omega_next >= omega )
			msg( __FILE__, __LINE__, "Omega increased." );
			cnvgWasSuccessful = false;
			break;
		end
		%
		vecX = vecX_next;
		vecF = vecF_next;
		matJ = matJ_next;
		omega = omega_next;
		vecG = matJ'*vecF;
		gNormSq = sumsq(vecG);
		if ( nargout >= 2 )
			datOut.iterCount_tot++;
			datOut.vecXVals(:,datOut.iterCount_tot+1) = vecX;
			datOut.vecFVals(:,datOut.iterCount_tot+1) = vecF;
			datOut.omegaVals(datOut.iterCount_tot+1) = omega;
			datOut.vecGVals(:,datOut.iterCount_tot+1) = vecG;
			datOut.gNormSqVals(datOut.iterCount_tot+1) = gNormSq;
			datOut.iterCount_cnvg = iterCount_cnvg;
		endif
	endwhile
	msgif( debugMode, __FILE__, __LINE__, sprintf( ...
	  "Convergence:  iterCount = %d,  bigDeltaNormSq = %0.3e,  omega = %0.3e,  gNormSq = %0.3e.",
	  iterCount_cnvg, sumsq(vecX-vecX0), omega, gNormSq ) );
	if ( ~cnvgWasSuccessful )
		msgif( debugMode, __FILE__, __LINE__, "Convergence failed." );
		return;
	end
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
%!	[ vecXF, datOut ] = findLocMin_alytJ( vecX0, funchFJ );
