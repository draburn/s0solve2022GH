% DRaburn 2022.05.02
%  zlinsolf100
%  First concrete attempt at MVP.

function [ vecX_best, vecF_best, datOut ] = zlinsolf100( funchF, vecX_initial, vecF_initial=[], datIn=[], prm=[] )
	%
	error( "END OF VALID CODE." );
	%
	initStuff;
	% such as...
	zlinsolf100_setConstants;
	prm.verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	prm.dbugLev = mygetfield( prm, "dbugLev", DBUGLEV__COPIOUS );
	fevalCount = 0;
	%
	vecX_best = vecX_initial;
	vecF_best = vecF_initial;
	datOut = [];
	%
	%
	%
	vecX = vecX_initial;
	vecF = vecF_initial;
	iterCount = 0;
	%
	while (1)
		iterCount++;
		if ( norm(vecF_best) <= fTol )
			msgif( prm.verbLev >= VERBLEV__MAIN, "SUCCESS: norm(vecF_best) <= fTol." );
			break;
		elseif ( iterCount > mygetfield( prm, "iterLimit", 100 ) )
			msgif( prm.verbLev >= VERBLEV__MAIN, "IMPOSED STOP: iterCount >= iterLimit." );
			break;
		endif
		%
		%
		%
		if ( bestIsGoodEnough )
			__moveToBest();
			continue;
		endif
		%
		%
		%
		if ( looksLikeWereNearABadLocalMin ) %% Needs refinement.
			__tryToEscapeBadLocalMin();
			continue;
		endif
		%
		%
		%
		if ( omegaMeanStrikeIsTooLargeAndHaveUnexploredSpace )
			__expandSpace();
			if ( vecWasLinearlyIndepdentOfLocalBasis )
				continue;
			endif
		endif
		%
		%
		%
		if ( strikePtLooksGood )
			__evalStep( vecYStrike );
			% The pt should only look good if better than current frontrunner.
			% Internally...
			%  if step is new best...
			%    set it as such
			%    and, if step was near edge of TR and model was very good...
			%      reduce B (increase TR size).
			%  if step is not new best...
			%    if uncertainty in that direction is not small...
			%      refresh local subspace in that direction.
			%    Then, if (refreshed) model is bad compared to actual feval...
			%      increase B (decrease TR size)
			continue;
		endif
		%
		if ( strikePtLooksMaybe )
			__refreshLocalSubspace( vecYStrike );
			if ( vecWasLinearlyIndepdentOfLocalBasis )
				% We reduce the uncertainty in that direction,
				% perhaps to zero.
				continue;
			endif
		endif
		%
		if ( strikePtAlmostGoodEnough )
			__refreshLocalSubspace( __precon( __fModel( vecYStrike ) ) );
			if ( vecWasLinearlyIndepdentOfLocalBasis )
				continue;
			endif
		endif
		%
		%
		%
		if ( approachPtLooksGood )
			__evalStep( vecYApproach );
			continue;
		endif
		%
		if ( approachPtLooksMaybe )
			__refreshLocalSubspace( vecYApproach );
			if ( vecWasLinearlyIndepdentOfLocalBasis )
				continue;
			endif
		endif
		%
		if ( approachPtAlmostGoodEnough )
			__refreshLocalSubspace( __precon( __fModel( vecYApproach ) ) );
			if ( vecWasLinearlyIndepdentOfLocalBasis )
				continue;
			endif
		endif
		%
		%
		%
		if ( exploreExpandSpaceSeemsPromising )
			__expandSpace();
			if ( vecWasLinearlyIndepdentOfLocalBasis )
				continue;
			endif
		endif
		%
		if ( refreshSeemsPromising )
			__refreshLocalSubspace( __krylovesque() );
			if ( vecWasLinearlyIndepdentOfLocalBasis )
				continue;
			endif
		endif
		%
		%
		%
		msgif( prm.verbLev >= VERBLEV__MAIN, "ALGORITHM BREAKDOWN: No worthwhile action." );
		break;
	endwhile
	%
	%
	%
	datOut.fevalCount = fevalCount;
return
endfunction
