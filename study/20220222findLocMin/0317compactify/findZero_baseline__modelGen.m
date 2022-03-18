function [ matJ, ary3Kappa, modelGen_datOut ] = findZero_baseline__modelGen( vecX, vecF, vecX_prev, vecF_prev, matJ_prev, funchF, modelGen_prm )
	modelGen_datOut = [];
	modelGen_datOut.fevalCount = 0;
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	%
	%
	% Update Jacobian.
	if ( isempty(matJ_prev) )
		% I expect this code is duplicated. Maybe de-duplicate later.
		epsFD = mygetfield( modelGen_prm, "epsFD", eps^0.25 );
		for n = 1 : sizeX
			vecXP = vecX; vecXP(n) += epsFD;
			vecXM = vecX; vecXM(n) -= epsFD;
			vecFP = funchF(vecXP); modelGen_datOut.fevalCount++;
			vecFM = funchF(vecXM); modelGen_datOut.fevalCount++;
			matJ(:,n) = (vecFP-vecFM)/(2.0*epsFD);
		endfor
	else
		assert( issize(matJ_prev,[sizeF,sizeX]) );
		switch (tolower(mygetfield( modelGen_prm, "jUpdate", "broyden"  )))
		case { "full" }
			% I expect this code is duplicated. Maybe de-duplicate later.
			epsFD = mygetfield( modelGen_prm, "epsFD", eps^0.25 );
			for n = 1 : sizeX
				vecXP = vecX; vecXP(n) += epsFD;
				vecXM = vecX; vecXM(n) -= epsFD;
				vecFP = funchF(vecXP); modelGen_datOut.fevalCount++;
				vecFM = funchF(vecXM); modelGen_datOut.fevalCount++;
				matJ(:,n) = (vecFP-vecFM)/(2.0*epsFD);
			endfor
		case { "broyden" }
			fooX = vecX - vecX_prev;
			fooF = vecF - (vecF_prev+matJ_prev*fooX);
			matJ = matJ_prev + fooF*(fooX')/(fooX'*fooX);
			if (mygetfield( modelGen_prm, "doUphillCheck", true ))
				vecG = matJ'*vecF;
				if ( norm(vecG) > 0.0 )
					vecX_trial = vecX - (eps^0.25)*vecG/norm(vecG);
					vecF_trial = funchF( vecX_trial ); modelGen_datOut.fevalCount++;
					if ( norm(vecF_trial) > norm(vecF) )
						msg( __FILE__, __LINE__, "Oh, snap! The 'gradient' was going uphill!" );
						% Do *another* Broyden update, to make "gradient" at least be downhill.
						fooX = vecX_trial - vecX;
						fooF = vecF_trial - (vecF+matJ_prev*fooX);
						matJ = matJ + fooF*(fooX')/(fooX'*fooX);
					endif
				else
					msg( __FILE__, __LINE__, "Gradient is apparently zero. What should I do?" );
				endif
			endif
		case { "pool" }
			error( "To-do: Implement pool update; this requires much data passing." );
		case { "none" }
			matJ = matJ_prev;
		otherwise
			error( "Invalid value of jUpdate." );
		endswitch
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
