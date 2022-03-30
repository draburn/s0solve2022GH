% Dev
%  Pure Broyden.

function [ vecXF, vecFF, datOut ] = findZero_105( vecX0, funchF, prm=[] )
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
	matJ = zeros( sizeF, sizeX );
	epsFD = mygetfield( prm, "epsFD", eps^0.4 );
	for n=1:sizeX
		vecXP = vecX + epsFD*matIX(:,n);
		vecFP = funchF( vecXP ); fevalCount++;
		matJ(:,n) = (vecFP-vecF)/(epsFD);
	endfor
	%
	matH = matJ'*matJ;
	vecG = matJ'*vecF;
	hNorm = sqrt(sum(sumsq(matH))/sizeX);
	assert( 0 ~= hNorm );
	cholTol = mygetfield( prm, "cholTol", 1e-6 );
	[ matR, cholFlag ] = chol(matH);
	if ( 0==cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) )
		matHRegu = matH;
	else
		matHRegu = matH + hNorm*sqrt(eps)*matIX;
		matR = chol(matHRegu); % Ensure it's pos-def.
		assert( min(diag(matR)) > max(abs(diag(matR)))*cholTol );
	endif
	vecDelta = matR \ ( (matR') \ (-vecG) );
	vecX_next = vecX + vecDelta;
	vecF_next = funchF(vecX_next); fevalCount++;
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
		msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
		  time()-time0, iterCount, fevalCount, ...
		  norm(matJ'*vecF), ...
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
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		vecX = vecX_next;
		vecF = vecF_next;
		%
		%
		%
		% Try with JA.
		matH = matJ'*matJ;
		vecG = matJ'*vecF;
		hNorm = sqrt(sum(sumsq(matH))/sizeX);
		assert( 0 ~= hNorm );
		[ matR, cholFlag ] = chol(matH);
		if ( 0==cholFlag && min(diag(matR)) > cholTol*max(abs(diag(matR))) )
			matHRegu = matH;
		else
			msg( __FILE__, __LINE__, "Need regu!" );
			matHRegu = matH + hNorm*sqrt(eps)*matIX;
			matR = chol(matHRegu); % Ensure it's pos-def.
			assert( min(diag(matR)) > max(abs(diag(matR)))*cholTol );
		endif
		vecDelta = matR \ ( (matR') \ (-vecG) );
		vecX_next = vecX + vecDelta;
		vecF_next = funchF(vecX_next); fevalCount++;
		%
		if (1)
		% Apply Broyden update.
		fooX = vecX_next - vecX;
		fooF = vecF_next - ( vecF + matJ*vecDelta );
		fooJ = fooF*(fooX')/(fooX'*fooX);
		matJ += fooJ;
		endif
	endwhile
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
return;
endfunction
