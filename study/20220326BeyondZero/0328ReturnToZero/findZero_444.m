% Dev
%  444 = 800 but do full Jacobian in first iter and update Jacobian within linsolf.
%  Note that _800__step may be identical to _700__step.
% THIS CODE IS BUGGY.

function [ vecXF, vecFF, datOut ] = findZero_444( vecX0, funchF, prm=[] )
	time0 = time();
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 ); fevalCount++;
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( 0~=norm(vecF0) );
	%
	%
	%
	matIX = eye(sizeX,sizeX);
	%
	%
	% Everything past here is iterated on.
	iterCount = 0;
	vecX_best = vecX0;
	vecF_best = vecF0;
	datOut.fNormVals(iterCount+1) = norm(vecF_best);
	datOut.fevalCountVals(iterCount+1) = fevalCount;
	datOut.iterCountVals(iterCount+1) = iterCount;
	%
	vecX = vecX0;
	vecF = vecF0;
	%
	%
	matJA = zeros( sizeF, sizeX );
	epsFD = mygetfield( prm, "epsFD", eps^0.4 );
	for n=1:sizeX
		vecXP = vecX + epsFD*matIX(:,n);
		vecFP = funchF( vecXP ); fevalCount++;
		matJA(:,n) = (vecFP-vecF)/(epsFD);
	endfor
	matW = matJA;
	matV = matIX;
	sizeV = sizeX;
	%
	dTreg = Inf; % Trust region size.
	msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "Initial dTreg = %f.", dTreg ) );
	initialFallRatio = [];
	findZero_800__step;
	initialFallRatio = max([ 0.1, norm(vecF_next)/norm(vecF) ]);
	%
	%
	%
	while (1)
		iterCount++;
		if ( norm(vecF_next) < norm(vecF_best) )
			vecX_best = vecX_next;
			vecF_best = vecF_next;
		endif
		datOut.fNormVals(iterCount+1) = norm(vecF_best);
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.iterCountVals(iterCount+1) = iterCount;
		%
		msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %4d, %10.3e;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
		  time()-time0, iterCount, fevalCount, ...
		  sizeV, norm(matW'*vecF), ...
		  norm(vecX_next-vecX0), norm(vecX_next-vecX0)-norm(vecX-vecX0), norm(vecX_next-vecX), ...
		  norm(vecF_next), norm(vecF)-norm(vecF_next), norm(vecF-vecF_next) ) );
		%
		%
		%
		fTol = mygetfield( prm, "fTol", eps );
		iterMax = mygetfield( prm, "iterMax", 50 );
		if ( norm(vecF_best) <= fTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: norm(vecF_best) <= fTol." );
			break;
		elseif ( norm(vecF_next) + fTol >= norm(vecF) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecF_next) + fTol >= norm(vecF)." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		vecX = vecX_next;
		vecF = vecF_next;
		%
		%
		%
		useCoasting = mygetfield( prm, "useCoasting", true );
		if (useCoasting)
		findZero_800__step;
		%
		% This criteria crude, but okay for now.
		if ( norm(vecF_next) < 0.5*norm(vecF) + 0.5*norm(vecFModel_pMax) )
			% Apply Broyden update.
			fooY = matV'*vecDelta;
			fooF = vecF_next - ( vecF + matW*fooY );
			fooW = fooF*(fooY')/(fooY'*fooY);
			matW += fooW;
			continue
		elseif ( norm(vecF_next) < norm(vecF) )
			vecX = vecX_next;
			vecF = vecF_next;
			% But, re-calculate Jacobian, below.
		endif
		endif
		%
		%
		%
		matJA += ( matW - (matJA*matV) ) * (matV');
		epsFD = mygetfield( prm, "epsFD", eps^0.4 );
		funchMatJProd = @(v)( ( funchF(vecX+epsFD*v) - vecF ) / epsFD );
		linsolf_prm = [];
		linsolf_prm.tol = mygetfield( prm, "linsolf_tol", 0.1*sqrt(norm(vecF)/norm(vecF0)) );
		%%
		if (1)
			linsolf_prm.matA = matJA;
			linsolf_prm = mygetfield( prm, "linsolf_prm", linsolf_prm );
			[ vecSSDeltaN, linsolf_datOut ] = linsolf_update( funchMatJProd, -vecF, linsolf_prm );
		else
			linsolf_prm.matP = pinv(matJA);
			linsolf_prm = mygetfield( prm, "linsolf_prm", linsolf_prm );
			[ vecSSDeltaN, linsolf_datOut ] = linsolf( funchMatJProd, -vecF, zeros(sizeX,1), linsolf_prm );
		endif
		%%
		fevalCount += linsolf_datOut.fevalCount;
		sizeV = size(linsolf_datOut.matV,2);
		matV = linsolf_datOut.matV;
		matW = linsolf_datOut.matW;
		%[ matLambda ] = eig( matW'*matW ); msg( __FILE__, __LINE__, sprintf( "phiMin = %e.", min(matLambda)/max(abs(matLambda)) ) );
		if ( 0 && verbLev >= VERBLEV__NOTIFY )
			msg( __FILE__, __LINE__, sprintf( "  [ rcond(A), frob(W-V)/frob(V), frob(W-A*V)/frob(V) ] = [ %f, %f, %f ].", ...
			  rcond(matJA), ...
			  sum(sumsq(matW-matV))/sum(sumsq(matV)), ...
			  sum(sumsq(matW-matJA*matV))/sum(sumsq(matV)) ) );
		endif
		%
		dTreg = Inf;
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, sprintf( "Reset dTreg = %f.", dTreg ) )
		initialFallRatio = [];
		findZero_800__step;
		initialFallRatio = norm(vecF_next)/norm(vecF);
		%
		% Apply Broyden update.
		fooY = matV'*vecDelta;
		fooF = vecF_next - ( vecF + matW*fooY );
		fooW = fooF*(fooY')/(fooY'*fooY);
		matW += fooW;
	endwhile
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
return;
endfunction


function p = __findPOfDeltaNorm( deltaNormMax, funchDeltaOfP )
	fzero_fun = @(p_dummy)( norm(funchDeltaOfP(p_dummy)) - deltaNormMax );
	if ( fzero_fun(1.0) <= 0.0 )
		p = 1.0;
		return;
	endif
	p = fzero( fzero_fun, [0.0, 1.0] );
return;
endfunction


function vecSSDelta = __funcSSDeltaOfP( p, matH, vecG )
	[ matR, cholFlag ] = chol( p*matH + (1.0-p)*eye(size(matH)) );
	assert( 0 == cholFlag );
	vecSSDelta = matR \ ( matR' \ (-p*vecG) );
return;
endfunction
