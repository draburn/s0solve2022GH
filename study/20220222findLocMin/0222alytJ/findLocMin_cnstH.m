% Function...
%  [?] = findLocMin_cnstH(?)
% Returns vecDelta corresponding to the local min of the omega model,
%  possibly subject to a trust region and a "reasonableness" check such as omega >= 0.

function [ vecDelta, datOut ] = findLocMin_cnstH( omega0, vecG, matH, prm=[], datIn=[] )
	%
	%
	% Parse input.
	sizeX = size(vecG0,1);
	debugMode = mygetfield( prm, "debugMode", true );
	if (debugMode)
		assert( isrealscalar(omega0) );
		assert( isrealarray(vecG0,[sizeX,1]) );
		assert( isrealarray(matH,[sizeX,sizeX]) );
		assert( issymmetric(matH) );
	endif
	gNormSq = sumsq(vecG);
	gNorm = sqrt(gNormSq);
	hNorm = sqrt(sum(sumsq(matH)));
	omegaModelMin = mygetfield( prm, "omegaModelMin", [] );
	deltaNormMax = mygetfield( prm, "deltaNormMax", 100.0*gNorm/hNorm );
	%
	%
	% Set default return values.
	vecDelta = zeros(sizeX,1);
	if (nargout>=2)
		datOut = [];
	endif
	%
	%
	% Handle "corner zero" cases.
	if ( isrealscalar(omegaModelMin) )
	if ( omegaModelMin >= omega0 )
		msgif( debugMode, __FILE__, __LINE__, "Initial objective value is already below specified minimum." );
		return;
	endif
	endif
	if ( 0.0 == gNorm )
		msgif( debugMode, __FILE__, __LINE__, "Initial gradient is zero." );
		return;
	endif
	if ( 0.0 == hNorm )
		if ( isrealscalar(omegaModelMin) && isrealscalar(deltaNormMax) )
			vecDelta = min([ (omega0-omegaModelMin)/gNormSq, deltaNormMax/gNorm ]) * (-vecG);
		elseif ( isrealscalar(omeagMin) )
			vecDelta = ((omega0-omegaModelMin)/gNormSq) * (-vecG);
		elseif ( isrealscalar(deltaNormMax) )
			vecDelta = (deltaNormMax/gNorm) * (-vecG);
		else
			msgif( debugMode, __FILE__, __LINE__, "Hessian is zero and there is no constraint on step." );
		endif
		return;
	endif
	%
	%
	% Find the largest valid step we can (or, at least, something close).
	% Consider letting an inital guess for mu be passed in.
	mu = 0.0;
	matM = matH + mu*matI;
	[ matR, cholFlag ] = chol( matM );
	if ( 0 ~= cholFlag )
		msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with mu zero failed; trying with a small mu." );
		muReguCoeff = 1e-5; % Arbitrary.
		mu = muReguCoeff*hMaxAbs;
		matM = matH + mu*matI;
		[ matR, cholFlag ] = chol( matM );
		if ( 0 ~= cholFlag )
			msgif( debugMode, __FILE__, __LINE__, "Cholesky factorization with small mu failed; finding muCrit." );
			msgif( debugMode, __FILE__, __LINE__, "Calling eig(). This may be slow. Faster approaches may be possible, such as:" );
			msgif( debugMode, __FILE__, __LINE__, "  start with mu = upper bound for eigenvalue of H, and target omega = 0.0; or," );
			msgif( debugMode, __FILE__, __LINE__, "  increase mu exponentially until chol() works." );
			[ matPsi_eig, matLambda_eig ] = eig( matH );
			muCrit = -min(diag(matLambda_eig));
			mu = muCrit + muReguCoeff * ( muCrit + hNorm );
			matM = matH + mu*matI;
			matR = chol( matM );
			cholFlag = 0;
			%
			clear matPsi_eig;
			clear matLambda_eig;
			clear muCrit;
		endif
		clear muReguCoeff;
	endif
	clear matM;
	vecDelta = -( matR \ (matR'\vecG) );
	%
	%
	% Apply deltaNormMax (the trust region).
	if ( isrealscalar(deltaNormMax) )
	if ( norm(vecDelta) > deltaNormMax )
		deltaNormMin = 0.6*deltaNormMax; % Arbitrary-ish.
		findLocMin_cnstH__deltaNorm;
	endif
	endif
	
	msg( __FILE__, __LINE__, "HACK OUTPUT RETURN!" ); vecDelta = -0.0001*vecG0; return;

	error( "GARBAGE CODE BELOW!" );
	
	
	
	%
	%
	constrainDeltaNorm = false;
	if ( isrealscalar(deltaNormMax) )
	if ( norm(vecDelta) > deltaNormMax )
		constrainDeltaNorm = true;
	endif
	endif
	constrainOmegaModel = false;
	if ( isrealscalar(omegaModelMin) )
	if ( omgeaModel < omegaModelMin )
		constrainOmegaModel = true;
	endif
	endif
	%
	appplyConstraint = ( constrainDeltaNorm || constrainOmegaModel );
	if ( applyConstraint )
		% Find a high value of mu that satisfies the constraints.
		muHi = 10.0*hNorm + 10.0*gNorm/deltaNormMax; % Arbitrary-ish.
		
		
		
		muLo = mu;
		vecDeltaOfMuLo = vecDelta;
		omegaOfMuLo = omega;
		
	endif
	
	error( "GARBAGE CODE BELOW!" );
	
	
	% Apply deltaNormMax (the trust region).
	if ( isrealscalar(deltaNormMax) )
	if ( norm(vecDelta) > deltaNormMax )
		msgif( debugMode, __FILE__, __LINE__, sprintf( "Applying deltaNormMax: deltaNormMax = %f; mu = %f; norm(vecDelta) = %f.", deltaNormMax, mu, norm(vecDelta) ) );
		%
		% First, do a one-shot targetting center.
		deltaNormMin = 0.6*deltaNormMax; % Arbitrary-ish.
		deltaNormTrgt = ( deltaNormMax + deltaNormMin ) / 2.0;
		% Model: ||delta||^2 = ( a / (b+mu) )^2,
		%  match ||delta||^2 and d/dmu(||delta||^2) at previous mu.
		vecDeltaPrime = -( matR \ ( matR' \ vecDelta ) );
		dsq = sumsq(vecDelta,1);
		ddsqdmu = 2.0*(vecDelta'*vecDeltaPrime);
		assert( 0.0 > ddsqdmu );
		b = 2.0*dsq/(-ddsqdmu) - mu;
		a = 2.0*(dsq^1.5)/(-ddsqdmu);
		%
		mu = (a/deltaNormTrgt) - b;
		matM = matH + mu*matI;
		matR = chol( matM );
		vecDelta = -( matR \ ( matR' \ vecG ) );
		%
		if ( norm(vecDelta) > deltaNormMax || norm(vecDelta) < deltaNormMin )
			% Do loop as needed
			error( "Not implemented" );
		endwhile
	endif
	endif
	
	error( "END." );
	omegaModel = omegaModel = omega0 + vecG'*vecDelta + 0.5*(vecDelta'*matH*vecDelta);
	%
	%
	% Apply omegaMin and, if relevant, note updated trust region size.
	deltaNormMaxSuggested = [];
	if ( isrealscalar(omegaMin) )
	if ( omegaModel < omegaMin )
		msgif( debugMode, __FILE__, __LINE__, sprintf( "Applying omegaMin: omegaMin = %f; mu = %f; omegaModel = %f.", omegaMin, mu, omegaModel ) );
		assert( omegaMin < omega0 );
		omegaTrgt = omegaMin;
		% Model: omegaModel = omega0 - ( g^2 / (c+mu) ).
		%  match omega at previous mu.
		mu = mu + normGSq*(omegaTrgt-omegaModel)/((omega0-omegaTrgt)*(omega0-omegaModel));
		matM = matH + mu*matI;
		matR = chol( matM );
		vecDelta = -( matR \ ( matR' \ vecG ) );
		%
		msgif( debugMode, __FILE__, __LINE__, sprintf( "  Result: mu = %f; omegaModel = %f.", mu, omegaModel ) );
		if ( omegaModel < omegaMin )
			msg( __FILE__, __LINE__, sprintf( "WARNING: Failed to satisfy omegaMin (%f vs %f).", omegaMin, omegaModel ) );
			msg( __FILE__, __LINE__, "  Consider writing a loop to fix this." );
		elseif ( omegaModel > (omega0+omegaMin)/2.0 )
			msg( __FILE__, __LINE__, sprintf( "WARNING: Considerably overshot omegaMin (%f vs %f).", omegaMin, omegaModel ) );
			msg( __FILE__, __LINE__, "  Consider writing a loop to fix this." );
		endif
		%
		deltaNormMaxSuggested = norm(vecDelta);
		%
		clear omegaTrgt;
	endif
	endif
	%
return;
endfunction
