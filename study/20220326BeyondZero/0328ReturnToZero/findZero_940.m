% Dev
%  9xx = approaching new JFNK (with phi-patch, etc),
%   for prototyping structure and refereshing memory.
%  940 = slisolf.
%
%  2022.04.28: Note that there is no coasting here.
% TODO:
%  Analysis: OSQU with Phi and Gamma???



function [ vecXF, vecFF, datOut ] = findZero_940( vecX0, funchF, prm=[] )
	time0 = time();
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = mygetfield( prm, "vecF0", [] );
	if (isempty(vecF0))
		vecF0 = funchF( vecX0 );
		fevalCount++;
	endif
	sizeF = size(vecF0,1);
	assert( 0 < sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	assert( 0~=norm(vecF0) );
	%
	%
	vecX_best = vecX0;
	vecF_best = vecF0;
	%
	%
	%
	%fModelType = mygetfield( prm, "fModelType", F_MODEL_TYPE__CONVENTIONAL );
	%fModelDat.modelType = fModelType;
	fModelDat.vecX = vecX0;
	fModelDat.vecF = vecF0;
	fModelDat.matA = eye(sizeF,sizeX); % Model for full-space Jacobian.
	fModelDat = mygetfield( prm, "fModelDat0", fModelDat );
	assert( reldiff(vecX0,fModelDat.vecX,eps) < eps );
	assert( reldiff(vecF0,fModelDat.vecF,eps) < eps );
	assert( isrealarray(fModelDat.matA,[sizeF,sizeX]) );
	%
	%
	%
	iterCount = 0;
	vecX = vecX0;
	vecF = vecF0;
	datOut.iterCountVals(iterCount+1) = iterCount;
	datOut.fevalCountVals(iterCount+1) = fevalCount;
	datOut.fNormVals(iterCount+1) = norm(vecF_best);
	datOut.vecXVals(:,iterCount+1) = vecX;
	datOut.vecFVals(:,iterCount+1) = vecF;
	%
	step_tol = sqrt(eps); % Use a tight solve on first iteration to get a large subspace.
	step_prm = mygetfield( prm, "step_prm", [] );
	step_prm.slinsolfver = mygetfield( prm, "slinsolfver", 100 );
	stepSearchDat = [];
	[ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, stepSearchDat, step_tol, step_prm );
	fevalCount += step_datOut.fevalCount;
	if (isempty(vecX_next))
		msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: __findStep() failed on starting values." );
		vecXF = [];
		vecFF = [];
		datOut.fevalCount = fevalCount;
		datOut.iterCount = iterCount;
		return;
	endif
	%
	%
	%
	while (1);
		iterCount++;
		if ( norm(vecF_next) < norm(vecF_best) )
			vecX_best = vecX_next;
			vecF_best = vecF_next;
		endif
		datOut.iterCountVals(iterCount+1) = iterCount;
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.fNormVals(iterCount+1) = norm(vecF_best);
		datOut.vecXVals(:,iterCount+1) = vecX_next;
		datOut.vecFVals(:,iterCount+1) = vecF_next;
		%
		sizeV = size( fModelDat_next.matV, 2 );
		msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d, %3d;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
		  time()-time0, iterCount, fevalCount, sizeV,...
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
		fModelDat = fModelDat_next;
		%
		%
		%
		% TODO: Add coasting.
		% If successful, "continue" back to start of loop.
		% If not, do following...
		%
		%
		%
		step_tol = max([ eps/norm(vecF), 0.1/norm(vecF0) ]);
		step_prm = mygetfield( prm, "step_prm", [] );
		step_prm.slinsolfver = mygetfield( prm, "slinsolfver", 100 );
		[ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, stepSearchDat, step_tol, step_prm );
		fevalCount += step_datOut.fevalCount;
		%
		if (isempty(vecX_next)||isempty(fModelDat_next))
			msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: __findStep() failed." );
			break;
		endif
	endwhile
	%
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
return;
endfunction



function [ vecX_next, vecF_next, fModelDat_next, stepSearchDat_next, step_datOut ] = __findStep( funchF, vecX, vecF, fModelDat, stepSearchDat, step_tol, step_prm )
	fevalCount = 0;
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	matA = fModelDat.matA;
	%
	%
	%
	% Interface with slinsolf().
	slinsolf_prm = step_prm;
	slinsolf_prm.dta_c0 = step_tol;
	slinsolf_datIn = [];
	slinsolf_datIn.preconDat.matA = matA;
	switch (mygetfield(step_prm,"slinsolfver",100))
	case 100
		[ vecX_next, vecF_next, slinsolf_datOut ] = slinsolf100( funchF, vecX, vecF, slinsolf_prm, slinsolf_datIn );
	case 200
		[ vecX_next, vecF_next, slinsolf_datOut ] = slinsolf200( funchF, vecX, vecF, slinsolf_prm, slinsolf_datIn );
	otherwise
		error( "Invalid slinsolfver." );
	endswitch
	fevalCount += slinsolf_datOut.fevalCount;
	%
	if (isempty(vecX_next))
		msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: slinsolf() failed." );
		fModelDat_next = [];
		stepSearchDat_next = [];
		step_datOut = [];
		step_datOut.slinsolf_datOut = slinsolf_datOut;
		step_datOut.fevalCount = fevalCount;
		return;
	endif
	if ( norm(vecX_next-vecX) < eps )
		msg( __FILE__, __LINE__, "ALGORITHM BREAKDOWN: slinsolf() failed." );
		fModelDat_next = [];
		stepSearchDat_next = [];
		step_datOut = [];
		step_datOut.slinsolf_datOut = slinsolf_datOut;
		step_datOut.fevalCount = fevalCount;
		return;
	endif
	%
	%
	%
	% Update matA.
	matV = slinsolf_datOut.localModelDat.matV;
	matW = slinsolf_datOut.localModelDat.matW;
	matPhi = slinsolf_datOut.localModelDat.matPhi;
	matGamma = slinsolf_datOut.localModelDat.matGamma;
	vecFModel_next = slinsolf_datOut.vecFModelF;
	sizeV = size(matV,2);
	sizePhi = size(matPhi,2);
	%
	matA += ( matW - (matA*matV) ) * (matV'); % 1: Basic multi-rank (Broydenesque) update.
	for n=1:sizePhi
		% 2: Advance model per quad terms.
		matA += matGamma(:,n) * ( matPhi(:,n)'*(vecX_next-vecX) ) * (matPhi(:,n)');
	endfor
	% 3: "On-step quadratic update", which is like Broyden x 2.
	fooX = vecX_next - vecX;
	fooF = vecF_next - vecFModel_next;
	%%%fooF = vecF_next - ( vecF + matA*fooX );
	matA += 2.0 * fooF * (fooX') / (fooX'*fooX);
	% TODO: This OSQU may be wrong when already including phi & gamma.
	%
	% Question: Is this equivalent to updating W then updating A???
	%
	%
	%
	fModelDat_next = [];
	fModelDat_next.matA = matA;
	fModelDat_next.matV = matV;
	fModelDat_next.maPhi = matPhi;
	%
	stepSearchDat_next = [];
	%
	step_datOut.fevalCount = fevalCount;
	step_datOut.sizeV = sizeV;
	step_datOut.slinsolf_datOut = slinsolf_datOut;
return;
endfunction
