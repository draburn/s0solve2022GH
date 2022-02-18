msg( __FILE__, __LINE__, "WIP!" );
vecX_next = vecX - 0.0001*vecG;
return;

%
matK = zeros(sizeX,sizeX);
trialCount = 0;
vecX_next = []; % Keep track of best so far.
vecF_next = []; % Keep track of best so far.
omega_next = []; % Keep track of best so far.
dGenMin = dTreg;
while (1)
	% Check pre-iter imposed limits.
	trialCount++;
	trialLimit = 100;
	if ( trialCount > trialLimit )
		msgif( debugMode, __FILE__, __LINE__, "Reached trialLimit." );
		return;
	endif
	%
	%
	%
	% Generate step
	genLimit = 100;
	matD = matI(sizeX,sizeX);
	%
	genCount = 0;
	mu = 0.0;
	matH = matJTJ + matK;
	while (1)
		genCount++;
		if ( genCount > genLimit )
			msgif( debugMode, __FILE__, __LINE__, "Reached genLimit." );
			return;
		endif
		%
		[ matR, cholFlag ] = chol( matH + mu*matD );
		if ( 0~=cholFlag )
			error( "TODO: Increase mu to guarantee pos-def." );
			continue;
		endif
		vecDelta_trial = -( matR \ ( matR' \ vecG ) );
		%
		if ( norm(vecDelta_trial) > dGenMin )
			error( "TODO: Increase mu to be under dGenMin." );
			continue;
		endif
		dGenMin = norm(vecDelta_trial);
		%
		omegaModel_trial = omega + vecG'*vecDelta_trial + 0.5*(vecDelta_trial'*matH*vecDelta_trial);
		if ( omegaModel_trial < -eps*abs(omega) )
			error( "TODO: Increase mu so that omegaModel is not negative." );
			continue;
		endif
		%
		%
		% Check if step is too small.
		% Should we also check omegaFallTol and stepSizeTol???
		wolfeC = 0.9; % MAKE PARAM
		if ( abs(vecDelta_trial'*(vecG+matH*vecDelta_trial)) >= wolfeC * abs(vecDelta_trial'*vecG) )
			msgif( debugMode, __FILE__, __LINE__, "Step is too small per Wolfe curvature condition. Giving up." );
			fooLHS = vecDelta_trial'*(vecG+matH*vecDelta_trial)
			fooRHS = vecDelta_trial'*vecG
		endif
		%
		break;
	endwhile
	vecFModel_trial = vecF + matJ*vecDelta_trial; % No quad info here.
	omegaModel_trial = omega + vecG'*vecDelta_trial + 0.5*(vecDelta_trial'*matH*vecDelta_trial); % Repeated. POITROME.
	vecGModel_trial = vecG + matH*vecDelta_trial; % Some quad info here.
	%
	%
	%
	vecX_trial = vecX + vecDelta_trial;
	vecF_trial = funchFJ( vecX_trial );
	if ( ~isrealarray(vecF_trial,[sizeF,1]) )
		msgif( debugMode, __FILE__, __LINE__, "Function evaluation failed." );
		if (~isempty(vecX_next))
			msgif( __FILE__, __LINE__, "Function evaluation failed even though it succeeded for an earlier trial step!" );
			msgif( __FILE__, __LINE__, "This suggests the invalid space is not simply at distant values of x." );
		endif
		dGenFailureBTFactor = 0.1; % MAKE PARAM.
		dGenMin *= dGenFailureBTFactor;
		msgif( debugMode, __FILE__, __LINE__, "Updating dTreg." );
		dTreg = dGenMin;
		continue;
	endif
	omega_trial = sumsq(vecF_trial,1)/2.0;
	%
	%
	%
	% Find how far vecF_trial is from vecF, vecFModel.
	% Since ||vecFModel-vecF|| <= 2.0*||vecF|| is guaranteed,
	%  we can simply compare against vecF.
	c_hugeChangeToF = 10.0; % MAKE PARAM.
	assert( c_hugeChangeToF >= 2.0 );
	if ( norm(vecF_trial-vecF) >= c_hugeChangeToF*norm(vecF) )
		msgif( debugMode, __FILE__, __LINE__, "Function was radically different at trial point." );
		if (~isempty(vecX_next))
			msgif( __FILE__, __LINE__, "Function was radically different even though it was reasonable for an earlier trial step!" );
			msgif( __FILE__, __LINE__, "This suggests the presence of a sepratrix." );
		endif
		dGenRadicallyWrongBTFactor = 0.1; % MAKE PARAM.
		dGenMin *= dGenRadicallyWrongBTFactor;
		dTreg = dGenMin;
		continue;
	endif
	%
	%
	%
	% Check outright acceptance criterion.
	if ( omega_trial <= (1.0-c_acceptance)*omega + c_acceptance*omegaModel_trial )
		msgif( debug, __FILE__, __LINE__, "Step is very good!" );
		if (~isempty(omega_next))
		if ( omega_trial >= omega_next )
			error( "Somehow this very good step is no better than an earlier one. This should be impossible." );
		endif
		endif
		vecX_next = vecX_trial;
		vecF_next = vecF_trial;
		omega_next = omega_trial;
		return;
	endif
	%
	%
	% Having gotten here, the step is reasonable enough that we can use the results to update matK.
	% We also conditionally might accept the step.
endwhile
