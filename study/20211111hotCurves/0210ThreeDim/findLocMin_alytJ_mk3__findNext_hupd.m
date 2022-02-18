msg( __FILE__, __LINE__, "WIP!" );
vecX_next = vecX - 0.0001*vecG;
return;

trialLimit = 100;
wolfeC = 0.9;
dGenFailureBTFactor = 0.1;
c_hugeChangeToF = 10.0;
dGenRadicallyWrongBTFactor = 0.1;
%
matK = zeros(sizeX,sizeX);
trialCount = 0;
vecX_next = []; % Keep track of best so far.
dGenMin = dTreg;
while (1)
	% Check pre-iter imposed limits.
	trialCount++;
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
		dGenMin *= dGenFailureBTFactor;
		msgif( debugMode, __FILE__, __LINE__, "Updating dTreg." );
		dTreg = dGenMin;
		continue;
	end
	omega_trial = sumsq(vecF_trial,1)/2.0;
	%
	%
	%
	% Check outright acceptance criterion.
	if ( omega_trial <= (1.0-c_acceptance)*omega + c_acceptance*omegaModel_trial )
		msgif( debug, __FILE__, __LINE__, "Step is very good!" );
		vecX_next = vecX_trial;
		return;
	end
	%
	%
	%
	% Find how far vecF_trial is from (1-b)*vecF + b*vecFModel, minimized over b.
	foo1 = ( vecF_trial - vecF )' * ( vecFModel_trial - vecF );
	foo2 = sumsq( vecFModel_trial - vecF );
	if ( foo1 <= 0.0 )
		b = 0.0;
	elseif ( foo1 >= foo2 )
		b = 1.0;
	else
		b = foo1 / foo2;
	endif
	vecRho = vecF_trial - ( (1.0-b)*vecF + b*vecFModel_trial );
	if ( norm(vecRho) > c_hugeChangeToF * norm(vecF) )
		msgif( debugMode, __FILE__, __LINE__, "Function was radially different at trial point." );
		dGenMin *= dGenRadicallyWrongBTFactor;
		dTreg = dGenMin;
		continue;
	end
	%
	%
	% Having gotten here, the step is reasonable enough that we can use the results to update matK.
	% We also conditionally might accept the step.
endwhile
