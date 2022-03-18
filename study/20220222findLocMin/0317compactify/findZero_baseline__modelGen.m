function [ matJ, ary3Kappa, modelGen_datOut ] = findZero_baseline__modelGen( vecX, vecF, vecX_prev, vecF_prev, matJ_prev, funchF, modelGen_prm )
	modelGen_datOut = [];
	modelGen_datOut.fevalCount = 0;
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	%
	%
	% Update Jacobian.
	if ( isempty(matJ_prev) || ~mygetfield( modelGen_prm, "useInexactJ", false  ) )
		epsFD = mygetfield( modelGen_prm, "epsFD", eps^0.25 );
		for n = 1 : sizeX
			vecXP = vecX; vecXP(n) += epsFD;
			vecXM = vecX; vecXM(n) -= epsFD;
			vecFP = funchF(vecXP); modelGen_datOut.fevalCount++;
			vecFM = funchF(vecXM); modelGen_datOut.fevalCount++;
			matJ(:,n) = (vecFP-vecFM)/(2.0*epsFD);
		endfor
	else
		% Note: These cases may not be full supported in "_baseline".
		assert( issize(matJ_prev,[sizeF,sizeX]) );
		switch (tolower(mygetfield( modelGen_prm, "inexactJType", "broyden" )))
		case { "broyden" }
			fooX = vecX - vecX_prev;
			fooF = vecF - (vecF_prev+matJ_prev*fooX);
			matJ = matJ_prev + fooF*(fooX')/(fooX'*fooX);
		case { "pool" }
			error( "To-do: Implement pool update; this requires much data passing." );
		case { "none" }
			matJ = matJ_prev;
		otherwise
			error( "Unsupported case." );
		endswitch
		%
		if (mygetfield( modelGen_prm, "applyGradientUpdate", true ))
			% DRaburn 2022.03.18:
			%  I'm sure other approaches are possible.
			%  Preliminary investigation of trying to do this along both
			%   the gradient and Newton directions didn't work so well.
			%  Since inexactJ is not the focus of "_baseline", I'm moving on.
			%  I don't think this guarantees that the new gradient be downhill,
			%   but, this seems to do the trick.
			vecG = matJ'*vecF;
			assert( 0.0 ~= norm(vecG) );
			vecX_trial = vecX - (sqrt(eps))*(1.0+norm(vecX))*vecG/norm(vecG);
			vecF_trial = funchF( vecX_trial ); modelGen_datOut.fevalCount++;
			fooX = vecX_trial - vecX;
			fooF = vecF_trial - (vecF+matJ*fooX);
			matJ = matJ + fooF*(fooX')/(fooX'*fooX);
		endif
	endif
	%
	ary3Kappa = zeros( sizeX, sizeX, sizeF );
	% *** NOTE THAT ARY3KAPPA IS NOW XXF, NOT FXX! ***
	%
	% DRaburn 2022.03.18...
	% Regarding Kappa:
	%  - There's a question of what directions to consider, what basis to use,
	%   and what "off-diagonal" terms to consider.
	%  - This is an active topic of research and subject to heavy revision.
	%  - In this, "_baseline", sovler, Kappa is only ever re-calcualted from scratch.
	%  - There may be cases where above calculation of the Jacobian could be skipped, but, POITROME.
	%
	% One thing to do is find our "priority directions".
	% These include the steepest-descent direction, the step to the local minimum of the linear model (if it exists),
	%  and any directions for which matJ'*matJ has a non-positive eigenvaule (if they exist).
	% We can also consider scaling, such as Marquardt ( diag(sqrt(1.0./diag(matJ'*matJ))) ) scaling.
return;
endfunction
