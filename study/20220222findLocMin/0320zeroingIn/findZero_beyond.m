% Goes beyond fsolve.

function [ vecXF, vecFF, datOut ] = findZero_beyond( vecX0, funchF, prm=[] )
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
	usePool = mygetfield( prm, "usePool", true );
	if ( usePool )
		poolDeltaX = [];
		poolDeltaF = [];
		poolJacobi = [];
	endif
	%
	vecX = [];
	vecF = [];
	matJ = [];
	trustRegionSize = [];
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
		if ( ~isempty(vecX) && ~isempty(vecF) && ~isempty(matJ) && ~isempty(trustRegionSize) )
			msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
			  time()-time0, iterCount, fevalCount, ...
			  norm(matJ'*vecF), trustRegionSize, ...
			  norm(vecX_next-vecX0), norm(vecX_next-vecX0)-norm(vecX-vecX0), norm(vecX_next-vecX), ...
			  norm(vecF_next), norm(vecF)-norm(vecF_next), norm(vecF-vecF_next) ) );
		else
			msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
			  time()-time0, iterCount, fevalCount, ...
			  0.0, 0.0, ...
			  norm(vecX_next-vecX0), 0.0, 0.0, ...
			  norm(vecF_next), 0.0, 0.0 ) );
		endif
		%
		if ( ~isempty(vecX) && ~isempty(vecF) && ~isempty(matJ) )
			assert( norm(vecF_next) < norm(vecF) );
			if ( usePool )
			if ( norm(vecX_next-vecX) > sqrt(eps) )
				poolTol = mygetfield( prm, "poolTol", 0.5 );
				poolDeltaX = [ vecX_next - vecX, poolDeltaX ];
				poolDeltaF = [ vecF_next - vecF, poolDeltaF ];
				[ matQ, rvecDrop ] = utorthdrop( poolDeltaX, poolTol );
				%
				% Apply drop vecs to matJSettled.
				for n=1:size( poolDeltaX, 2 )
				if ( rvecDrop(n) )
					fooX = poolDeltaX(:,n);
					fooF = poolDeltaF(:,n) - poolJacobi*fooX;
					poolJacobi += fooF*(fooX')/(fooX'*fooX);
				endif
				endfor
				poolDeltaX = poolDeltaX(:,~rvecDrop);
				poolDeltaF = poolDeltaF(:,~rvecDrop);
				%
				% Apply remaining vectors as a multi-rank orthonormal update.
				if ( size( poolDeltaX, 2 ) > 0 )
					matR = matQ'*poolDeltaX;
					matW = poolDeltaF / matR;
					matJ = poolJacobi*( matIX - matQ*(matQ') ) + matW * (matQ');
				else
					matJ = poolDat.matJSettled;
				endif
			endif
			else
				% Apply Broyden update.
				vecDeltaX = vecX_next - vecX;
				vecDeltaF = vecF_next - vecF - (matJ*vecDeltaX);
				matJ += vecDeltaF*(vecDeltaX')/(vecDeltaX'*vecDeltaX);
			endif
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
			assert( ~isempty(trustRegionSize) );
			% Try to find a good next guess using approximate Jacobian.
			vecG = matJ'*vecF;
			matH = matJ'*matJ;
			vecDeltaC = (-vecG) * (vecG'*vecG) / (vecG'*matH*vecG);
			vecDeltaN = matH \ (-vecG);
			%
			% First shot.
			vecDelta = calcDogLeg( vecDeltaC, vecDeltaN, trustRegionSize );
			vecX_trial = vecX + vecDelta;
			vecF_trial = funchF( vecX_trial ); fevalCount++;
			if ( norm(vecF_trial) < norm(vecF) )
				if ( norm( vecF + matJ*vecDelta - vecF_trial ) < 0.1 * norm( vecF ) )
					msgif( verbLev >= VERBLEV__PROGRESS+10, __FILE__, __LINE__, "Model was very accurate on first try." );
					trustRegionSize *= 2.0;
				endif
				msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Trial with approximate Jacobian was successful on first try." );
				vecX_next = vecX_trial;
				vecF_next = vecF_trial;
				continue;
			endif
			%
			if (0)
			% Second shot.
			trustRegionSize /= 5.0;
			vecDelta = calcDogLeg( vecDeltaC, vecDeltaN, trustRegionSize );
			vecX_trial = vecX + vecDelta;
			vecF_trial = funchF( vecX_trial ); fevalCount++;
			if ( norm(vecF_trial) < norm(vecF) )
				msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Trial with approximate Jacobian was successful." );
				vecX_next = vecX_trial;
				vecF_next = vecF_trial;
				continue;
			endif
			endif
		endif
		%
		%
		%
		msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Trial with approximate Jacobian was unsuccessful; recalculating full local Jacobiain." );
		matJ = zeros( sizeF, sizeX );
		fdOrder = mygetfield( prm, "fdOrder", 1 );
		switch (fdOrder)
		case 1
			epsFD = mygetfield( prm, "epsFD", eps^0.4 );
			for n=1:sizeX
				vecFP = funchF( vecX + epsFD*matIX(:,n) ); fevalCount++;
				matJ(:,n) = (vecFP-vecF)/epsFD;
			endfor
		case 2
			epsFD = mygetfield( prm, "epsFD", eps^0.25 );
			for n=1:sizeX
				vecFP = funchF( vecX + epsFD*matIX(:,n) ); fevalCount++;
				vecFM = funchF( vecX - epsFD*matIX(:,n) ); fevalCount++;
				matJ(:,n) = (vecFP-vecFM)/(2.0*epsFD);
			endfor
		otherwise
			error( "Invalid case." );
		endswitch
		if ( usePool )
			poolDeltaX = matIX;
			poolDeltaF = matJ;
			%poolDeltaX = [];
			%poolDeltaF = [];
			poolJacobi = matJ;
		endif
		%
		%
		%
		% Try to find a good next guess using actual local Jacobian.
		vecG = matJ'*vecF;
		matH = matJ'*matJ;
		vecDeltaC = (-vecG) * (vecG'*vecG) / (vecG'*matH*vecG);
		vecDeltaN = matH \ (-vecG);
		if ( isempty(trustRegionSize) )
			trustRegionSize = norm(vecDeltaN);
		endif
		%
		haveGoodTrial = false;
		while (~haveGoodTrial)
			vecDelta = calcDogLeg( vecDeltaC, vecDeltaN, trustRegionSize );
			vecX_trial = vecX + vecDelta;
			vecF_trial = funchF( vecX_trial ); fevalCount++;
			if ( norm(vecF_trial) < norm(vecF) )
				msgif( verbLev >= VERBLEV__COPIOUS, __FILE__, __LINE__, "Trial with exact Jacobian was successful." );
				vecX_next = vecX_trial;
				vecF_next = vecF_trial;
				haveGoodTrial= true;
				break;
			else
				trustRegionSize /= 5.0;
				if ( trustRegionSize < 1e-8 )
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


function vecDelta = calcDogLeg( vecDeltaC, vecDeltaN, stepSize )
	if ( stepSize >= norm(vecDeltaN) )
		vecDelta = vecDeltaN;
		return;
	elseif ( stepSize <= norm(vecDeltaC) )
		vecDelta = vecDeltaC*stepSize/norm(vecDeltaC);
		return;
	endif
	vecA = vecDeltaC;
	vecB = vecDeltaN-vecDeltaC;
	atb = vecA'*vecB;
	asq = vecA'*vecA;
	bsq = vecB'*vecB;
	discrim = atb^2 - (asq-stepSize^2)*bsq;
	assert( discrim >= 0.0 );
	s = ( sqrt(discrim) - atb ) / bsq;
	assert( s >= 0.0 - sqrt(eps) );
	assert( s <= 1.0 + sqrt(eps) );
	vecDelta = vecA + s*vecB;
	assert( reldiff(norm(vecDelta),stepSize) < sqrt(eps) );
endfunction
