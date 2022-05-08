% DRaburn 2022.05.02
%  zlinsolf100
%  First concrete attempt at MVP.

function [ vecX_best, vecF_best, datOut ] = zlinsolf100( funchF, vecX_initial, vecF_initial=[], datIn=[], prmIn=[] )
	%
	error( "END OF VALID CODE." );
	%
	initStuff;
	% such as...
	zlinsolf100_setConstants;
	prm = prmIn;
	prm.verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	prm.dbugLev = mygetfield( prm, "dbugLev", DBUGLEV__COPIOUS );
	%
	fevalCount = 0;
	vecX_best = vecX_initial;
	vecF_best = vecF_initial;
	datOut = [];
	%
	vecX = vecX_initial;
	vecF = vecF_initial;
	iterCount = 0;
	%
	while (1)
		iterCount++;
		%
		% Simple stoping criteria.
		if ( sumsq(vecF_best) <= fTol^2 )
			msgif( prm.verbLev >= VERBLEV__MAIN, "SUCCESS: sumsq(vecF_best) <= fTol^2." );
			break;
		endif
		%
		iterLimit = mygetfield( prm, "iterLimit", ceil(20 + 10*sqrt(sizeX) + 0.01*sizeX) );
		if ( iterCount > iterLimit )
			msgif( prm.verbLev >= VERBLEV__MAIN, "IMPOSED STOP: iterCount >= iterLimit." );
			break;
		endif
		%
		if ( stopsignalpresent()) )
			msgif( prm.verbLev >= VERBLEV__MAIN, "IMPOSED STOP: Received stop signal." );
			break;
		endif
		%
		%
		% Niche cases.
		if ( 0 == sizeV )
			__expandSpace( vecF );
			fevalCount += ?;
			continue;
		endif
		%
		c = 0.5;
		bestIsAcceptable = (  sumsq(vecF_best) < sumsq(vecF) - c*( sumsq(vecF) - sumsq(vecF_idealUnbound) );
		if ( bestIsAcceptable )
			__acceptPoint( vecX_best, vecF_best );
			fevalCount += ?;
			continue;
		endif
		%
		if ( hessianIsSingular )
			msgif( prm.verbLev >= VERBLEV__NOTICE, "TODO: Implement handling for hessianIsSingular." );
			break;
		endif
		%
		%
		% We're in striking distance?
		idealUnbound_isWithinTR = (  norm(matB*vecY_idealUnbound) < 1.0  );
		if ( idealUnbound_isWithinTR )
			c = 0.5;
			plausiblyWithinTol = (  sumsq(vecF_idealUnbound) > c*(fTol^2)  );
			c = 0.5;
			likelyToDecrease = (  sumsq(vecF_idealUnbound) + fSumSq_deltaW_idealUnbound < c*(fTol^2)  );
			c = 0.5;
			pbIsNotMuchBetter = (  c*(sumsq(vecF_idealUnbound) + fSumSq_deltaW_idealUnbound) < sumsq(vecF_pessimisticBound) + fSumSq_deltaW_pessismisticBound  );
			if ( plausiblyWithinTol && likelyToDecrease && pbIsNotMuchBetter )
				__tryPoint( vecY_idealUnbound );
					% If bad: if needed, refresh W; if (refreshed) fModel is bad, increase B.
					% If excellent: if norm(matB*vecY_pessimisticBound) == 1 (within fp), decrease B.
				fevalCount += ?;
				if ( norm(vecF_trial) < norm(vecF_best) )
					vecX_best = vecX_trial; % = vecX + matV * vecY.
					vecF_best = vecF_trial;
				endif
				continue
			endif;
			% else, see below.
		endif
		%
		%
		% The routine cases.
		c = 0.5;
		needLargerSubspace = (  sumsq(vecF_idealUnbound) > c * (fTol^2)  );
		% This needLargerSubspace is based on the idea that we'll eventually need an appropriately large subspace,
		%  so we'd do better by expanding to that dimensionality early, to use as a preconditioner.
		if ( needLargerSubspace )
			__expandSubspace( vecF_idealUnbound );
			fevalCount += ?;
			continue;
		endif
		%
		c = 1.0e-2;
		fSumSqFallThresh = c*sumsq(vecF);
		noGoodStepInTR = (  sumsq(vecF_idealBound) > sumsq(vecF) - ffSumSqFallThresh  );
		if ( idealBoundIsUnacceptable )
			msgif( prm.verbLev >= VERBLEV__NOTICE, "IMPOSED STOP: noGoodStepInTR." );
			break;
		endif
		%
		c = 0.5;
		fSumSqFallThresh = c * ( sumsq(vecF) - sumsq(vecF_idealBound) );
		fSumSq_upperBoundEst_pessimisticBound = sumsq(vecF_pessimisticBound) + fSumSq_deltaW_pessismisticBound;
		needBetterAim = (  fSumSq_upperBoundEst_pessimisticBound > sumsq(vecF) - fSumSqFallThresh  );
		if ( needBetterAim )
			__refreshSubspace( vecY_pessimisticBound );
			fevalCount += ?;
			continue;
		endif
		%
		%
		% Everything looks good. Try the step.
		__tryPoint( vecY_pessimisticBound );
			% If bad: if needed, refresh W; if (refreshed) fModel is bad, increase B.
			% If excellent: if norm(matB*vecY_pessimisticBound) == 1 (within fp), decrease B.
		fevalCount += ?;
		if ( norm(vecF_trial) < norm(vecF_best) )
			vecX_best = vecX_trial; % = vecX + matV * vecY.
			vecF_best = vecF_trial;
		endif
		continue;
	endwhile
	%
	%
	%
	datOut.fevalCount = fevalCount;
return
endfunction
