% Dev

function [ vecXF, vecFF, datOut ] = findZero_shouldIterMin( vecX0, funchF, prm=[] )
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
	bestIsNot0 = false;
	vecX_best = vecX0;
	vecF_best = vecF0;
	%
	vecX = vecX0;
	vecF = vecF0;
	%
	matJ = zeros( sizeF, sizeX );
	epsFD = mygetfield( prm, "epsFD", eps^0.3 );
	for n=1:sizeX
		vecXP = vecX + epsFD*matIX(:,n);
		vecXM = vecX - epsFD*matIX(:,n);
		vecFP = funchF( vecXP ); fevalCount++;
		vecFM = funchF( vecXM ); fevalCount++;
		matJ(:,n) = (vecFP-vecFM)/epsFD;
	endfor
	%
	matH = matJ'*matJ;
	vecG = matJ'*vecF;
	hNorm = sqrt(sum(sumsq(matH))/sizeX);
	assert( 0 ~= hNorm );
	[ matR, cholFlag ] = chol(matH);
	if ( 0==cholFlag )
		matHRegu = matH;
	else
		matHRegu = matH + hNorm*sqrt(eps)*matIX;
		matR = chol(matHRegu); % Ensure it's pos-def.
	endif
	%
	funchDeltaOfP = @(p) ( p*matHRegu + (1.0-p)*matIX ) \ (-p*vecG);
	funchFNormOfP = @(p) norm(funchF(vecX+funchDeltaOfP(p)));
	fminbnd_options = optimset( "TolX", eps^2, "TolFun", eps^2 );
	[ fminbnd_x, fminbnd_fval, fminbnd_info, fminbnd_output ] = fminbnd( funchFNormOfP, 0.0, 1.0, fminbnd_options );
	fevalCount += fminbnd_output.funcCount;
	p = fminbnd_x;
	vecX_next = vecX + funchDeltaOfP(p);
	vecF_next = funchF(vecX_next); fevalCount++;
	%
	%
	%
	iterCount =1;
	while (1)
		if ( norm(vecF_next) < norm(vecF_best) )
			vecX_best = vecX_next;
			vecF_best = vecF_next;
		endif
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
		elseif ( norm(vecF_next) + fTol >= norm(vecF) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecF_next) + fTol >= norm(vecF)." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		iterCount++;
		vecX = vecX_next;
		vecF = vecF_next;
		%
		%
		%
		matJ = zeros( sizeF, sizeX );
		epsFD = mygetfield( prm, "epsFD", eps^0.3 );
		for n=1:sizeX
			vecXP = vecX + epsFD*matIX(:,n);
			vecXM = vecX - epsFD*matIX(:,n);
			vecFP = funchF( vecXP ); fevalCount++;
			vecFM = funchF( vecXM ); fevalCount++;
			matJ(:,n) = (vecFP-vecFM)/epsFD;
		endfor
		%
		matH = matJ'*matJ;
		vecG = matJ'*vecF;
		hNorm = sqrt(sum(sumsq(matH))/sizeX);
		assert( 0 ~= hNorm );
		[ matR, cholFlag ] = chol(matH);
		if ( 0==cholFlag )
			matHRegu = matH;
		else
			matHRegu = matH + hNorm*sqrt(eps)*matIX;
			matR = chol(matHRegu); % Ensure it's pos-def.
		endif
		%
		funchDeltaOfP = @(p) ( p*matHRegu + (1.0-p)*matIX ) \ (-p*vecG);
		funchFNormOfP = @(p) norm(funchF(vecX+funchDeltaOfP(p)));
		fminbnd_options = optimset( "TolX", eps^2, "TolFun", eps^2 );
		[ fminbnd_x, fminbnd_fval, fminbnd_info, fminbnd_output ] = fminbnd( funchFNormOfP, 0.0, 1.0, fminbnd_options );
		fevalCount += fminbnd_output.funcCount;
		p = fminbnd_x;
		%
		vecX_next = vecX + funchDeltaOfP(p);
		vecF_next = funchF(vecX_next); fevalCount++;
	endwhile
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
return;
endfunction
