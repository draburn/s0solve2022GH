% Dev

function [ vecXF, vecFF, datOut ] = findZero_basic( vecX0, funchF, prm=[] )
	time0 = time();
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 ); fevalCount++;
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	matIX = eye(sizeX,sizeX);
	%
	%
	vecX_best = vecX0;
	vecF_best = vecF0;
	%
	%
	vecX = [];
	vecF = [];
	matJ = [];
	vecX_next = vecX0;
	vecF_next = vecF0;
	iterCount = 0;
	datOut = [];
	while (1)
		if ( norm(vecF_next) < norm(vecF_best) )
			% Update best.
			vecX_best = vecX_next;
			vecF_best = vecF_next;
		endif
		%
		if ( ~isempty(vecX) && ~isempty(vecF) && ~isempty(matJ) )
			msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
			  time()-time0, iterCount, fevalCount, ...
			  norm(matJ'*vecF), ...
			  norm(vecX_next-vecX0), norm(vecX_next-vecX0)-norm(vecX-vecX0), norm(vecX_next-vecX), ...
			  norm(vecF_next), norm(vecF)-norm(vecF_next), norm(vecF-vecF_next) ) );
		else
			msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
			  time()-time0, iterCount, fevalCount, ...
			  0.0, ...
			  norm(vecX_next-vecX0), 0.0, 0.0, ...
			  norm(vecF_next), 0.0, 0.0 ) );
		endif
		
		if ( ~isempty(vecX) && ~isempty(vecF) && ~isempty(matJ) )
			assert( norm(vecF_next) < norm(vecF) );
			% Apply Broyden update.
			vecDeltaX = vecX_next - vecX;
			vecDeltaF = vecF_next - vecF - (matJ*vecDeltaX);
			matJ += vecDeltaF*(vecDeltaX')/(vecDeltaX'*vecDeltaX);
		endif
		% Move to next.
		iterCount++;
		vecX = vecX_next;
		vecF = vecF_next;
		datOut.iterCountVals(iterCount) = iterCount-1;
		datOut.fevalCountVals(iterCount) = fevalCount;
		datOut.fNormVals(iterCount) = norm(vecF_best);
		%
		%
		% Check stop.
		fTol = mygetfield( prm, "fTol", eps );
		iterMax = mygetfield( prm, "iterMax", 50 );
		if ( norm(vecF) <= fTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: norm(vecF) <= fTol." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		%
		%
		%
		if ( ~isempty(matJ) )
			% Try to find a good next guess using approximate Jacobian.
			vecG = matJ'*vecF;
			matH = matJ'*matJ;
			matR = chol( matH );
			vecDelta = - ( matR \ (matR'\vecG) );
			vecX_trial = vecX + vecDelta;
			vecF_trial = funchF( vecX_trial ); fevalCount++;
			if ( norm(vecF_trial) < norm(vecF) )
				msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Trial with approximate Jacobian was successful." );
				vecX_next = vecX_trial;
				vecF_next = vecF_trial;
				continue;
			endif
			clear vecX_trial;
			clear vecF_trial;
		endif
		%
		%
		%
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Trial with approximate Jacobian was unsuccessful; recalculating full local Jacobiain." );
		matJ = zeros( sizeF, sizeX );
		epsFD = mygetfield( prm, "epsFD", eps^0.4 );
		for n=1:sizeX
			vecXP = vecX + epsFD*matIX(:,n);
			vecFP = funchF( vecXP ); fevalCount++;
			matJ(:,n) = (vecFP-vecF)/epsFD;
		endfor
		%
		%
		%
		% Try to find a good next guess using approximate Jacobian.
		vecG = matJ'*vecF;
		matH = matJ'*matJ;
		matR = chol( matH );
		vecDelta = - ( matR \ (matR'\vecG) );
		p = 1.0;
		haveGoodTrial = false;
		while (~haveGoodTrial)
			vecX_trial = vecX + p*vecDelta;
			vecF_trial = funchF( vecX_trial ); fevalCount++;
			if ( norm(vecF_trial) < norm(vecF) )
				msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Trial with exact Jacobian was successful." );
				vecX_next = vecX_trial;
				vecF_next = vecF_trial;
				haveGoodTrial= true;
				break;
			else
				p /= 5.0;
				if ( p < 1e-8 )
					msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Step size got overly small." );
					break;
				endif
			endif
		endwhile
		%
		if ( haveGoodTrial )
			continue;
		else
			break;
		endif
	endwhile
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
return;
endfunction
