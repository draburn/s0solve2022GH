function [ vecX_final, datOut ] = linsolf_sja_looptr( funchMatAProd, vecB, vecX0, prm = [] )
	%error( "NOT FULLY IMPLEMENTED!" );
	time0 = time();
	fevalCount = 0;
	mydefs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__FLAGGED );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	fevalCount = 0;
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	if ( 0.0 ~= norm(vecX0) )
		warning( "vecX0 does nothing." );
	endif
	sizeF = size(vecB,1);
	assert( isrealarray(vecB,[sizeF,1]) );
	assert( sizeX == sizeF );
	%
	if ( 0.0 == norm(vecB) )
		vecX = zeros(sizeX,1);
		datOut.vecY = [];
		datOut.matV = zeros(sizeX,0);
		datOut.matW = zeros(sizeF,0);
		datOut.vecZ = zeros(sizeF,0);
		datOut.vecRho = zeros(sizeF,0);
		return;
	endif
	%
	tol = mygetfield( prm, "tol", 0.1 );
	assert( isrealscalar(tol) );
	assert( 0.0 <= tol );
	assert( tol <= 1.0 );
	% We'll replace "tol" with "tol_bound" below.
	%
	iterMax = mygetfield( prm, "iterMax", sizeX );
	assert( isposintscalar(iterMax) );
	%.
	cholTol = mygetfield( prm, "cholTol", 1e-6 );
	assert( isrealscalar(cholTol) );
	assert( 0.0 <= cholTol );
	assert( cholTol <= 1.0 );
	% DRaburn 2022-10-16: Shouldn't "cholTol" be "cholThresh"???
	%
	matP = mygetfield( prm, "matP", [] );
	useSJA = mygetfield( prm, "useSJA", true );
	assert( isbool(useSJA) );
	assert( isscalar(useSJA) );
	if ( useSJA )
		sja_prm = mygetfield( prm, "sja_prm", [] );
	endif
	sja_matJA = [];
	sja_matJAInv = [];
	%
	%
	% Trust Region stuff.
	tol_bound = mygetfield( prm, "tol_bound", tol );
	clear tol;
	tol_actual = mygetfield( prm, "tol_actual", (1.0+tol_bound)/2.0 );
	tol_unbound = mygetfield( prm, "tol_unbound", 100.0*eps );
	btCoeff = mygetfield( prm, "btCoeff", 0.2 );
	stepInverseTRLimit = mygetfield( prm, "stepInverseTRLimit", 0.0 );
	assert( isrealscalar( tol_bound ) );
	assert( isrealscalar( tol_unbound ) );
	assert( isrealscalar( btCoeff ) );
	assert( isrealscalar( stepInverseTRLimit ) );
	assert( 0.0 <= tol_unbound );
	assert( tol_unbound <= tol_bound );
	assert( tol_bound <= tol_actual );
	assert( tol_actual <= 1.0 );
	assert( 0.0 < btCoeff );
	assert( btCoeff < 1.0 );
	assert( 0.0 <= stepInverseTRLimit );
	boundStepPrm = mygetfield( prm, "boundStepPrm", [] );
	%
	stepRelTol = mygetfield( prm, "stepRelTol", 1.0e-8 );
	assert( isrealscalar(stepRelTol) );
	assert( 0.0 <= stepRelTol );
	assert( stepRelTol <= 1.0 );
	%
	%
	%
	% Everything past here is iterated upon.
	vecX_final = [];
	vecF_final = [];
	vecU = vecB;
	if ( ~isempty(matP) )
		vecU = matP*vecU;
	endif
	vecV = __orth( vecU, [], prm );
	assert( isrealarray(vecV,[sizeX,1]) );
	vecW = funchMatAProd( vecV ); fevalCount++;
	assert( isrealarray(vecW,[sizeF,1]) );
	matV = vecV;
	matW = vecW;
	successiveBTCount = 0;
	vecY_unbound = [];
	vecX_unbound = [];
	vecRho_unbound = vecB;
	rho_unbound = 1.0;
	vecY_bound = [];
	vecX_bound = [];
	vecRho_bound = vecB;
	rho_bound = 1.0;
	while (1)
		sizeV = size(matV,2);
		assert( reldiff( matV'*matV, eye(sizeV,sizeV) ) < sqrt(eps) );
		%
		%
		%
		matH = matW' * matW;
		[ matR, cholFlag ] = chol( matH );
		if ( 0 ~= cholFlag || min(diag(matR)) < cholTol*max(abs(diag(matR))) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Cholesky factorization of Hessian failed." );
			sizeV
			diag(matH)'
			break;
		endif
		vecG = matW' * vecB;
		%sizeV
		%stepInverseTRLimit
		%
		vecY_unbound = matR \ ( matR' \ vecG );
		vecX_unbound = matV * vecY_unbound;
		vecRho_unbound = vecB - matW*vecY_unbound; % "vecRho" and "rho" have inconsistent normalization!
		rho_unbound = norm( vecRho_unbound ) / norm(vecB);
		%
		vecY_bound = __getYBound( matR, vecG, vecY_unbound, stepInverseTRLimit, boundStepPrm );
		if ( isempty(vecY_bound) || ~issize(vecY_bound,[sizeV,1]) )
			msgif( verbLev >= VERBLEV__FLAG, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: __getYBound() failed." );
			sizeV
			break;
		endif
		if ( norm(vecY_bound)*stepInverseTRLimit > 1.0 + 100.0*eps )
			norm(vecY_bound)
			stepInverseTRLimit
			norm(vecY_bound)*stepInverseTRLimit
		endif
		assert( norm(vecY_bound)*stepInverseTRLimit <= 1.0 + 100.0*eps );
		vecX_bound = matV * vecY_bound;
		vecRho_bound = vecB - matW*vecY_bound; % "vecRho" and "rho" have inconsistent normalization!
		rho_bound = norm( vecRho_bound ) / norm(vecB);
		%
		if ( rho_bound <= tol_bound || rho_unbound <= tol_unbound || sizeV >= iterMax )
			vecF_actual = prm.funchFshift( vecX_bound );
			fevalCount++;
			rho_actual = norm(vecF_actual) / norm(vecB);
			if ( rho_actual <= tol_actual )
				vecX_final = vecX_bound;
				vecF_final = vecF_actual;
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: rho_actual <= tol_actual." );
				break;
			elseif ( norm(vecY_bound) <= stepRelTol*norm(vecY_unbound) )
				vecX_final = vecX_bound;
				vecF_final = vecF_actual;
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecY_bound) <= stepRelTol*norm(vecY_unbound)." );
				break;
			else
				%stepInverseTRLimit
				%btCoeff
				%norm(vecY_bound)
				stepInverseTRLimit = 1.0 / (btCoeff*norm(vecY_bound));
				%stepInverseTRLimit
				successiveBTCount++;
				assert( successiveBTCount < 100 );
				continue;
			endif
		endif
		successiveBTCount = 0;
		%
		%
		%
		if ( sizeV >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: sizeV >= iterMax." );
			break;
		endif
		if ( stopsignalpresent() )
			msgif( verbLev >= VERBLEV__FLAGGED, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			break;
		endif
		%
		vecU = vecW;
		if ( ~isempty(matP) )
			vecU = matP*vecU;
		endif
		vecV = __orth( vecU, matV, prm );
		if ( isempty(vecV) && ~isempty(sja_matJAInv) )
			% Let's try again...
			vecV = __orth( sja_matJAInv*vecB, matV, prm );
		endif
		if ( isempty(vecV) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: Data became linearly dependent." );
			break;
		endif
		assert( isrealarray(vecV,[sizeX,1]) );
		vecW = funchMatAProd( vecV ); fevalCount++;
		assert( isrealarray(vecW,[sizeF,1]) );
		matV = [ matV, vecV ];
		matW = [ matW, vecW ];
		%
		%if ( useSJA && ~isempty(sja_matJA) )
		%	msg( __FILE__, __LINE__, "Looks like sja_matJA didn't work. Kaboom!" );
		%	sja_matJA = [];
		%	sja_matJAInv = [];
		%endif
		%
		if ( useSJA && isempty(sja_matJA) )
			sizeK = size(matV,2);
			%%%sizeK_hidden = 1;
			%%%sizeK_margin = 1;
			%%%sizeK_hidden = ceil(sqrt(sizeK)/2.0);
			%%%sizeK_margin = ceil(sqrt(sizeK)/2.0);
			%%%sizeK_hidden = ceil(sqrt(sizeK/2.0));
			%%%sizeK_margin = ceil(sqrt(sizeK/2.0));
			%
			%sizeK_hidden = ceil(sqrt(sizeK));
			%sizeK_margin = 0;
			sizeK_hidden = ceil(sqrt(sizeK)/3.0);
			sizeK_margin = ceil(sqrt(sizeK)/3.0);
			%
			sizeK_pass = sizeK - sizeK_hidden;
			sizeK_nze = sizeK_pass - sizeK_margin;
			%%%sizeK_nze = max([ sizeK_pass - sizeK_margin, sqrt(N) ]);
			%%%if ( sizeK_nze >= 1 )
			if ( sizeK_nze >= 1 && 0 == mod(sizeK_nze,30) )
				sja_prm.maxNumNZEPerRow = sizeK_nze;
				sja_prm.abortOnBadRow = true;
				%%%sja_tol = 1.0e-3;
				sja_tol = 1.0e-2;
				%
				sjaMethod = mygetfield( prm, "sjaMethod", "corr" );
				switch (tolower(sjaMethod))
				case { "corr" }
					msg( __FILE__, __LINE__, "Using sja_corr()." );
					[ temp_matJA, sja_datOut ] = sja_corr( matV(:,1:sizeK_pass), matW(:,1:sizeK_pass), sja_prm );
				case { "omp" }
					msg( __FILE__, __LINE__, "Using \"sja_fast()\" (which is OMP)." );
					[ temp_matJA, sja_datOut ] = sja_fast( matV(:,1:sizeK_pass), matW(:,1:sizeK_pass), sja_prm );
				otherwise
					error( "Unsupported value of sjaMethod." );
				endswitch
				%%%[ temp_matJA, sja_datOut ] = sja_corr_oneshot( matV(:,1:sizeK_pass), matW(:,1:sizeK_pass), sja_prm );
				%[ sum(sumsq( temp_matJA*matV - matW )), sum(sumsq( matW )) ]
				%[ sum(sumsq( temp_matJA*matV(:,sizeK_pass+1:end) - matW(:,sizeK_pass+1:end) )), sum(sumsq( matW(:,sizeK_pass+1:end) )) ]
				if ( ~isempty(temp_matJA) ...
				 && sum(sumsq( temp_matJA*matV - matW )) < sja_tol^2 * sum(sumsq( matW )) ...
				 && sum(sumsq( temp_matJA*matV(:,sizeK_pass+1:end) - matW(:,sizeK_pass+1:end) )) < sja_tol^2 * sum(sumsq( matW(:,sizeK_pass+1:end) )) )
					msg( __FILE__, __LINE__, sprintf( "  Tentatively captured sparse Jacobian! ( %d / %d / %d ).", sizeK_nze, sizeK_pass, sizeK ) );
					%matV
					%sja_matJA
					%figure( 300 );
					%imagesc(sja_matJA);
					if ( rcond(temp_matJA) > sqrt(eps)  )
						temp_matJAInv = inv(temp_matJA);
						if ( sum(sumsq( matV - temp_matJAInv*matW )) < sja_tol^2 * sum(sumsq( matV )) ...
						  && sum(sumsq( matV(:,sizeK_pass+1:end) - temp_matJAInv*matW(:,sizeK_pass+1:end) )) ...
						   < sja_tol^2 * sum(sumsq( matV(:,sizeK_pass+1:end) )) )
							msg( __FILE__, __LINE__, "  Captured sparse Jacobian!" );
							sja_matJA = temp_matJA;
							matP = temp_matJAInv;
							sja_matJAInv = temp_matJAInv;
						else
							msg( __FILE__, __LINE__, "  Rejected sparse Jacobian on inverse check." );
						endif
					endif
				endif
				temp_matJA = [];
				temp_matJAInv = [];
			endif
		endif
	endwhile
	if ( isempty(vecX_final) )
		vecX_final = vecX_bound;
		vecF_final = [];
	endif
	if ( 2 <= nargout )
		datOut.fevalCount = fevalCount;
		%datOut.vecY = vecY;
		datOut.matV = matV;
		datOut.matW = matW;
		%datOut.vecZ = matW*vecY;
		datOut.vecRho = vecRho_bound;
		datOut.sja_matJA = sja_matJA;
		datOut.sja_matJAInv = sja_matJAInv;
		datOut.vecF_final = vecF_final;
	endif
return;
end


function [ vecV ] = __orth( vecU, matVPrev, prm=[] )
	uSq = sumsq(vecU);
	if ( uSq == 0.0 )
		vecV = [];
		return;
	endif
	%
	%
	vecV = vecU/sqrt(uSq);
	if (isempty(matVPrev))
		return;
	endif
	%
	numPasses = mygetfield( prm, "numPasses", 2 );
	for n=1:numPasses
		vecV -= matVPrev*(matVPrev'*vecV);
		vSq = sumsq(vecV);
		relTolVSq = mygetfield( prm, "relTolVSq", eps );
		if ( vSq <= relTolVSq*uSq )
			vecV = [];
			return;
		endif
		vecV /= sqrt(vSq);
	endfor
return;
endfunction


function vecY = __getY( matW, vecB, prm=[] )
	matH = matW'*matW;
	sizeV = size(matW,2);
	hNorm = sqrt(sum(sumsq(matH))/sizeV);
	assert( 0 ~= hNorm );
	cholTol = mygetfield( prm, "cholTol", 1e-6 );
	[ matR, cholFlag ] = chol(matH);
	if ( 0==cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		matHRegu = matH;
	else
		matHRegu = matH + hNorm*sqrt(eps)*eye(sizeV,sizeV);
		[ matR, cholFlag ] = chol(matHRegu);
		assert( 0 == cholFlag );
		assert( min(diag(matR)) > max(abs(diag(matR)))*cholTol );
	endif
	vecY = matR \ ( matR' \ (matW'*vecB) );
return;
endfunction


function [ vecY_bound ] = __getYBound( matR, vecG, vecY_newton, stepInverseTRLimit, boundStepPrm )
	if ( stepInverseTRLimit*norm(vecY_newton) <= 1.0 )
		% Even full Newton step is within trust region, so, take it.
		vecY_bound = vecY_newton;
		%msg( __FILE__, __LINE__, sprintf( ...
		%  "foo = %f = %f * %f", ...
		%  norm(vecY_bound)*stepInverseTRLimit, ...
		%  norm(vecY_bound), ...
		%  stepInverseTRLimit ) );
		return;
	endif
	p = (vecG'*vecG) / sumsq(matR*vecG);
	vecY_cauchy = -p*vecG;
	vecY_diff = vecY_newton - vecY_cauchy;
	%
	if ( stepInverseTRLimit*norm(vecY_cauchy) >= 1.0 )
		% Even full Cauchy step goes outside trust region, so, take a fraction of it.
		vecY_bound = vecY_cauchy/(stepInverseTRLimit*norm(vecY_cauchy));
		%msg( __FILE__, __LINE__, sprintf( ...
		%  "foo = %f = %f * %f", ...
		%  norm(vecY_bound)*stepInverseTRLimit, ...
		%  norm(vecY_bound), ...
		%  stepInverseTRLimit ) );
		return;
	endif
	%
	assert( 0.0 < stepInverseTRLimit );
	assert( isfinite(stepInverseTRLimit) );
	assert( isreal(stepInverseTRLimit) );
	assert( isscalar(stepInverseTRLimit) );
	targetStepSize = ( 1.0 - 100.0*eps ) / stepInverseTRLimit;
	a = sumsq(vecY_diff);
	b = 2.0*(vecY_cauchy'*vecY_diff);
	c = sumsq(vecY_cauchy) - (targetStepSize^2);
	discrim = (b^2) - (4.0*a*c);
	if ( discrim < 0.0 || a < 0.0 )
		if ( true )
			msg( __FILE__, __LINE__, "INTERNAL ERROR: Dog-leg was poorly behaved." );
			msg( __FILE__, __LINE__, "  Info dump..." );
			size(matR)
			sumsq(vecY_newton)
			sumsq(vecY_cauchy)
			sumsq(vecY_diff)
			stepInverseTRLimit
			a
			b
			c
			discrim
			vecY_bound = [];
			msg( __FILE__, __LINE__, "  End info dump." );
		endif
		return;
	endif
	p = (-b+sqrt(discrim))/(2.0*a);
	vecY_bound = vecY_cauchy + (p*vecY_diff);
	if ( norm(vecY_bound)*stepInverseTRLimit >= 1.0 + 10.0*eps )
		if ( true )
			msg( __FILE__, __LINE__, "INTERNAL ERROR: Dog-leg was poorly behaved." );
			msg( __FILE__, __LINE__, "  Info dump..." );
			size(matR)
			sumsq(vecY_newton)
			sumsq(vecY_cauchy)
			sumsq(vecY_diff)
			stepInverseTRLimit
			a
			b
			c
			discrim
			vecY_bound = [];
			msg( __FILE__, __LINE__, "  End info dump." );
		endif
	endif
	%msg( __FILE__, __LINE__, sprintf( ...
	%  "foo = %f = %f * %f", ...
	%  norm(vecY_bound)*stepInverseTRLimit, ...
	%  norm(vecY_bound), ...
	%  stepInverseTRLimit ) );
	assert( norm(vecY_bound)*stepInverseTRLimit <= 1.0 + 100.0*eps );
	assert( abs(norm(vecY_bound)-targetStepSize) < 0.01*targetStepSize );
return;
endfunction
