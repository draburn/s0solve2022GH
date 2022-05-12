% DRaburn 2022.05.07
%  zlinsolf100
%  First concrete attempt at MVP.
%
% DRaburn 2022.05.09...
%   Oh, I've been forgetting to update this change list. Whoops.
%
% DRaburn 2022.05.08...
%   - Overhaulled trust region / boundary stuff.
%
% DRaburn 2022.05.07...
%   Todo:
%    - Implement partial-quadratic update to A when taking a step.
%    - Refactor trial step handling.
%   Soon:
%    - Implement halding of BLM / overly small TR as momentum / restart jump, consider quadratic terms?
%    - Allow for a constant matrix preconditioner.
%    - Allow dog-leg intstead of Levenberg curve.
%    - Refactor from scratch, test, compare, etc.
%   Much later:
%    - Optimize Octave code.
%    - Code in C/C++.
%    - Try with PIES.
%    - Test against available solvers.
%   Post-MVP:
%    ~ Reconsider anything and everything.
%    - Sepratrix domain enumeration for "basin" BLM handling?
%    - Non-smooth functions: automatically adjust epsFD for differentiation, what else?
%    - Automatic scaling in X, (makes subspace basis non-orthogonal)?

function [ vecX_best, vecF_best, datOut ] = zlinsolf100( funchF, vecX_initial, vecF_initial=[], prm=[] )
	% Init...
	datOut = [];
	fevalCount = 0;
	if (isempty(vecF_initial))
		vecF_initial = funchF( vecX_initial );
		fevalCount++;
	endif
	prm = __initPrm( vecX_initial, vecF_initial, prm );
	%
	msgif( prm.msgCopious, __FILE__, __LINE__, "Initializing subspace." );
	[ fModelDat, datOut_initModel ] = __initModel( funchF, vecX_initial, vecF_initial, prm );
	fevalCount += datOut_initModel.fevalCount;
	%
	% Current local vecX and vecF are stored only in fModelDat to avoid redundancy.
	iterCount = 0;
	vecX_cand = []; % Candidate for next vecX.
	vecF_cand = [];
	vecX_best = vecX_initial;
	vecF_best = vecF_initial;
	%
	stepCount = 0;
	datOut.iterCountOfSteps(stepCount+1) = 1;
	%
	while (1)
		%
		datOut.iterCountVals(iterCount+1) = iterCount;
		datOut.fevalCountVals(iterCount+1) = fevalCount;
		datOut.stepCountVals(iterCount+1) = stepCount;
		datOut.fNormVals(iterCount+1) = norm(vecF_best);
		datOut.vecXVals(:,iterCount+1) = fModelDat.vecX;
		datOut.vecFVals(:,iterCount+1) = fModelDat.vecF;
		iterCount++;
		%
		%
		fModelDat = __analyzeModel( fModelDat, prm );
		if ( fModelDat.omegaModelVarIU < -sqrt(eps)*fModelDat.omega ...
		  || fModelDat.omegaModelVarIB < -sqrt(eps)*fModelDat.omega ...
		  || fModelDat.omegaModelVarPB < -sqrt(eps)*fModelDat.omega )
			msgif( prm.msgNotice, __FILE__, __LINE__, "WARNING omegaModelVar (of at least one flavor) was significantly negative." );
			[ fModelDat.omega, fModelDat.omegaModelAvgPB, fModelDat.omegaModelVarIU, fModelDat.omegaModelVarIB, fModelDat.omegaModelVarPB ]
			fModelDat.omegaModelVarIU = max([ fModelDat.omegaModelVarIU, 0.0 ]);
			fModelDat.omegaModelVarIB = max([ fModelDat.omegaModelVarIB, 0.0 ]);
			fModelDat.omegaModelVarPB = max([ fModelDat.omegaModelVarPB, 0.0 ]);
		endif
		%
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( ...
		  "  iter: %d/%d;  step: %d;  feval: %d;  omega: %0.2e/%0.2e;  sizeV: %d/%d/%dx%d;  sizeB: %d.", ...
		  iterCount, prm.iterMax, stepCount, fevalCount, sumsq(fModelDat.vecF)/2.0, prm.omegaTol, ...
		  size(fModelDat.matVLocal,2), size(fModelDat.matV,2), size(vecX_initial,1), size(vecF_initial,1), size(fModelDat.matB,2) ) );
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( ...
		  "  omega IU: %8.2e ~ %8.2e (%8.2e);  IB: %8.2e ~ %8.2e (%8.2e);  PB: %8.2e ~ %8.2e (%8.2e).", ...
		  fModelDat.omegaModelAvgIU, fModelDat.omegaModelVarIU, fModelDat.omega-fModelDat.omegaModelAvgIU, ...
		  fModelDat.omegaModelAvgIB, fModelDat.omegaModelVarIB, fModelDat.omega-fModelDat.omegaModelAvgIB, ...
		  fModelDat.omegaModelAvgPB, fModelDat.omegaModelVarPB, fModelDat.omega-fModelDat.omegaModelAvgPB ) );
		msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( ...
		  "  step  IU: %8.2e @ %8.2e b %8.2e;  IB: %8.2e @ %8.2e b %8.2e;  PB: %8.2e @ %8.2e b %8.2e.", ...
		  norm(fModelDat.vecYIU), fModelDat.lIU, fModelDat.bIU, ...
		  norm(fModelDat.vecYIB), fModelDat.lIB, fModelDat.bIB, ...
		  norm(fModelDat.vecYPB), fModelDat.lPB, fModelDat.bPB ) );
		%
		if ( fModelDat.bPB > 1.0 + sqrt(eps) || fModelDat.bIB > 1.0 + sqrt(eps) )
			msg( __FILE__, __LINE__, "ERROR: fModelDat.bPB or bIB > 1.0 + sqrt(eps)" );
			break;
		endif
		if ( fModelDat.omegaModelAvgIU > fModelDat.omega*(1.0+sqrt(eps)) ...
		  || fModelDat.omegaModelAvgIB > fModelDat.omega*(1.0+sqrt(eps)) ...
		  || fModelDat.omegaModelAvgPB > fModelDat.omega*(1.0+sqrt(eps)) )
			msg( __FILE__, __LINE__, "ERROR: omegaModelAvg (of at least one flavor) > omega*(1.0+sqrt(eps))." );
			break;
		endif
		if (mygetfield(fModelDat,"haltAfterNextReport",false))
			error( "HALT!" );
		endif
		%
		% Simple stoping criteria.
		if ( norm(vecF_best) <= prm.fTol )
			msgif( prm.msgMain, __FILE__, __LINE__, "SUCCESS: sumsq(vecF_best) <= prm.fTol^2." );
			break;
		endif
		%
		if ( iterCount > prm.iterMax )
			msgif( prm.msgMain, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= prm.iterMax." );
			if (prm.msgCopious)
				msgif( prm.msgCopious, __FILE__, __LINE__, "vvvvv Data dump..." );
				__dumpModel( fModelDat, prm );
				msgif( prm.msgCopious, __FILE__, __LINE__, "^^^^^ End data dump." );
			endif
			break;
		endif
		%
		if (prm.debugMode)
			if (~isempty(fModelDat.matVLocal))
				assert( reldiff( fModelDat.matVLocal, fModelDat.matV*(fModelDat.matV'*fModelDat.matVLocal), eps ) < sqrt(eps) );
				foo = size(fModelDat.matVLocal,2);
				assert( reldiff( eye(foo,foo), fModelDat.matVLocal'*fModelDat.matVLocal, eps ) < sqrt(eps) );
				clear foo;
				%
				rd = reldiff( fModelDat.matWLocal, fModelDat.matW*(fModelDat.matV'*fModelDat.matVLocal) );
				if ( rd > sqrt(eps) )
					error( "matWLocal does not agree with matW." );
				endif
				%
				if (0)
				fooZ1 = fModelDat.matV'*fModelDat.matVLocal*(fModelDat.matVLocal'*fModelDat.vecYIU);
				if ( norm(fooZ1) > 0.0 )
					fooZ1 *= sqrt(eps)/norm(fooZ1);
					fooZ2 = 2.0*fooZ1;
					fooM1 = fModelDat.vecF + fModelDat.matW*fooZ1;
					fooM2 = fModelDat.vecF + fModelDat.matW*fooZ2;
					fooF1 = funchF( fModelDat.vecX + fModelDat.matV*fooZ1 );
					fooF2 = funchF( fModelDat.vecX + fModelDat.matV*fooZ2 );
					rho1 = norm(fooM1-fooF1);
					rho2 = norm(fooM2-fooF2);
					assert( rho2 < 1.5*rho1 );
				endif
				endif
			endif
		endif
		%
		if ( stopsignalpresent() )
			msgif( prm.msgMain, __FILE__, __LINE__, "IMPOSED STOP: Received stop signal." );
			break;
		endif
		%
		%%%if ( fModelDat.omegaModelAvgIU > prm.omegaTol )
		if ( fModelDat.omegaModelAvgIU > fModelDat.omega/100.0 )
		%%%if ( fModelDat.omegaModelAvgIU > max([ prm.omegaTol, 0.1*fModelDat.omega*(fModelDat.omega/(sumsq(vecF_initial)/2.0)) ]) )
		%%%if (  ( 0==stepCount && fModelDat.omegaModelAvgIU > prm.omegaTol ) ...
		%%%  ||  ( 0<stepCount  && fModelDat.omegaModelAvgIU > max([ prm.omegaTol, 0.1*fModelDat.omega*(fModelDat.omega/(sumsq(vecF_initial)/2.0)) ]) )  )
		if ( size(fModelDat.matV,2) < size(vecX_initial,1) )
			% Conceptually, we could also consider "refreshing" something already in our space;
			%  alternatively, we could be more conservative about expanding the subspace.
			% This is an area for future analysis.
			msgif( prm.msgCopious, __FILE__, __LINE__, "Expanding subspace." );
			[ fModelDat, datOut_expandModel ] = __expandModel( fModelDat.vecFModelIU, funchF, fModelDat, prm );
			fevalCount += datOut_expandModel.fevalCount;
			continue;
		endif
		endif
		%
		if ( fModelDat.bIU <= 1.0 && fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU <= prm.omegaTol )
			msgif( prm.msgCopious, __FILE__, __LINE__, "Trying ideal unbound step." );
			vecX_trial = fModelDat.vecXIU;
			vecF_trial = funchF( vecX_trial );
			fevalCount++;
			omega_trial = sumsq(vecF_trial)/2.0;
			%msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  omega_trial = %g.", omega_trial ) );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  >> trial... omega: %8.2e to %8.2e, expected %8.2e ~ %8.2e;", ...
			  fModelDat.omega, omega_trial, ...
			  fModelDat.omegaModelAvgPB, fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB ) );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "              fall: %8.2e, expected %8.2e ~ %8.2e.", ...
			  fModelDat.omega-omega_trial, ...
			  fModelDat.omega - fModelDat.omegaModelAvgPB, fModelDat.omega - (fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB) ) );
			if ( norm(vecF_trial) < norm(vecF_best) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Step is new best." );
				vecX_best = vecX_trial;
				vecF_best = vecF_trial;
			endif
			if ( ~isempty(vecF_cand) )
			if ( norm(vecF_trial) >= norm(vecF_cand) )
				msgif( prm.msgNotice, __FILE__, __LINE__, "Current trial is worse than earlier candidate; forcing acceptance of earlier candidate." );
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_cand)/2.0, fModelDat.omega-sumsq(vecF_cand)/2.0, fevalCount ) );
				fModelDat = __moveTo( vecX_cand, vecF_cand, fModelDat, prm );
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			endif
			%
			avefaThresh = mygetfield( prm, "avefaThresh", 0.5 ); % Actual vs expect fall acceptace threshold
			assert( 0.0 < avefaThresh );
			assert( avefaThresh < 1.0 );
			%msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  expected fall = %g.", fModelDat.omega - (fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU) ) );
			%msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  actual fall = %g.", fModelDat.omega - omega_trial ) );
			if ( omega_trial <= fModelDat.omega - avefaThresh*( fModelDat.omega - (fModelDat.omegaModelAvgIU + fModelDat.omegaModelVarIU) ) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Accepting step." );
				excellentThresh = mygetfield( prm, "excellentThresh", 0.1 );
				if ( norm(vecF_trial-fModelDat.vecFModelIU) <= excellentThresh*norm(fModelDat.vecFModelIU) )
					msgif( prm.msgCopious, __FILE__, __LINE__, "  Model was very accurate; removing wall(s)." );
					bRemoveFactor = mygetfield( prm, "bRemoveFactor", 1.5 );
					fModelDat = __removeB( bRemoveFactor*fModelDat.vecYIU, fModelDat, prm );
				endif
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_trial)/2.0, fModelDat.omega-sumsq(vecF_trial)/2.0, fevalCount ) );
				fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm );
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			msgif( prm.msgCopious, __FILE__, __LINE__, "  Rejecting step." );
			%
			if ( norm(vecF_trial) < norm(fModelDat.vecF) )
				% Trial is better than current, at least, so it's a candidate.
				vecX_cand = vecX_trial;
				vecF_cand = vecF_trial;
			endif
			%
			% Consider adding a new barrier/wall.
			% But, first, make sure uncertainty in W was not the issue.
			vecY = fModelDat.vecYIU;
			vecU = fModelDat.matV*vecY;
			vecV = __calcOrthonorm( vecU, fModelDat.matVLocal, prm );
			refreshedTrialStep = false;
			if ( 0.0 ~= norm(vecV) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Refreshing trial step before adding wall." );
				[ fModelDat, datOut_refresh ] = __refresh( vecY, funchF, fModelDat, prm );
				fevalCount += datOut_refresh.fevalCount;
				refreshedTrialStep = true;
			endif
			vecFModel = fModelDat.vecF + fModelDat.matW*vecY;
			badThresh = mygetfield( prm, "badThresh", 0.5 );
			if ( norm( vecF_trial - vecFModel ) > badThresh * norm(vecFModel) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Model was very bad; adding wall." );
				bAddFactor = mygetfield( prm, "bAddFactor", 0.5 );
				fModelDat = __addB( bAddFactor*vecY, fModelDat, prm );
			elseif (~refreshedTrialStep)
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Did not refresh trial step but model was bed enough; adding wall." );
				bAddFactor = mygetfield( prm, "bAddFactor", 0.5 );
				fModelDat = __addB( bAddFactor*vecY, fModelDat, prm );
			endif
			%
			clear vecX_trial;
			clear vecF_trial;
			clear omega_trial;
			continue;
		endif
		%
		minRelFallThresh = mygetfield( prm, "minRelFallThresh", 1.0E-4 );
		if ( fModelDat.omegaModelAvgIB > fModelDat.omega*(1.0-minRelFallThresh) )
			msgif( prm.msgMain, __FILE__, __LINE__, "We seem to have no way to reduce omega much." );
			msgif( prm.msgMain, __FILE__, __LINE__, "  This is expected to happen near a bad local minimum." );
			msgif( prm.msgMain, __FILE__, __LINE__, "  Todo: add handling for this case." );
			break;
		endif
		%
		%
		%
		practicalRelFallThresh = mygetfield( prm, "practicalRelFallThresh", 0.5 );
		if ( fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB <= fModelDat.omega - practicalRelFallThresh*(fModelDat.omega-fModelDat.omegaModelAvgIB) )
		if ( fModelDat.omegaModelAvgPB <= 2.0 * fModelDat.omegaModelAvgIB )
			msgif( prm.msgCopious, __FILE__, __LINE__, "Trying practical bound step." );
			%
			%
			%%%
			makePlotsAndHalt = false;
			if (1)
			if ( prm.msgCopious && prm.debugMode )
			makePlotsAndHalt = ( size(fModelDat.matVLocal,2)>=1 && size(fModelDat.matB,2)>=20 );
			if ( makePlotsAndHalt )
				vecY = fModelDat.vecYPB;
				fModelDat.haltAfterNextReport = true;
				msg( __FILE__, __LINE__, "Making plots before refreshing and/or adding wall!" );
				numFigs = 0;
				%
				%
				%
				numPts = 101;
				%sPts = linspace(0.99999,0.9999,numPts);
				sPts = linspace(0.9999,0.0,numPts);
				sPts = ( 1.0 - sPts.^2 ).^2;
				%%%sPts = linspace( 0.99, 1.00, numPts );
				sizeV = size(fModelDat.matV,2);
				sizeF = size(fModelDat.vecF,1);
				vecMG = -(fModelDat.matW'*fModelDat.vecF);
				matH0 = fModelDat.matW'*fModelDat.matW;
				matH1 = matH0 + fModelDat.matA;
				matSCurve = eye(sizeV,sizeV);
				%
				for n=1:numPts
					vecYIPts(:,n) = ( sPts(n)*matH0 + (1.0-sPts(n))*matSCurve ) \ (sPts(n)*vecMG);
					vecYPPts(:,n) = ( sPts(n)*matH1 + (1.0-sPts(n))*matSCurve ) \ (sPts(n)*vecMG);
					vecFModelIPts(:,n) = fModelDat.vecF + fModelDat.matW*vecYIPts(:,n);
					vecFModelPPts(:,n) = fModelDat.vecF + fModelDat.matW*vecYPPts(:,n);
					omegaModelAvgIPts(n) = sumsq(vecFModelIPts(:,n))/2.0;
					omegaModelAvgPPts(n) = sumsq(vecFModelPPts(:,n))/2.0;
					omegaModelPVarIPts(n) = omegaModelAvgIPts(:,n) + 0.5*vecYIPts(:,n)'*fModelDat.matA*vecYIPts(:,n);
					omegaModelPVarPPts(n) = omegaModelAvgPPts(:,n) + 0.5*vecYPPts(:,n)'*fModelDat.matA*vecYPPts(:,n);
					vecFActualIPts(:,n) = funchF( fModelDat.vecX + fModelDat.matV*vecYIPts(:,n) );
					vecFActualPPts(:,n) = funchF( fModelDat.vecX + fModelDat.matV*vecYPPts(:,n) );
					omegaActualIPts(n) = sumsq(vecFActualIPts(:,n))/2.0;
					omegaActualPPts(n) = sumsq(vecFActualPPts(:,n))/2.0;
					rhoIPts(n) = sumsq( vecFActualIPts(:,n) - vecFModelIPts(:,n) )/2.0;
					rhoPPts(n) = sumsq( vecFActualPPts(:,n) - vecFModelPPts(:,n) )/2.0;
					bIPts(n) = max(abs(vecYIPts(:,n)'*fModelDat.matB));
					bPPts(n) = max(abs(vecYPPts(:,n)'*fModelDat.matB));
					%lIPts(n) = (eps+norm(fModelDat.matVLocal'*fModelDat.matV*vecYIPts(:,n)))/(eps+norm(vecYIPts(:,n)));
					%lPPts(n) = (eps+norm(fModelDat.matVLocal'*fModelDat.matV*vecYPPts(:,n)))/(eps+norm(vecYPPts(:,n)));
					lIPts(n) = norm(fModelDat.matVLocal'*fModelDat.matV*vecYIPts(:,n))/norm(vecYIPts(:,n));
					lPPts(n) = norm(fModelDat.matVLocal'*fModelDat.matV*vecYPPts(:,n))/norm(vecYPPts(:,n));
				endfor
				%
				%
				indexOfPB = numPts;
				while (bPPts(indexOfPB)>1.0)
					indexOfPB--;
				endwhile
				%
				numFigs++; figure(numFigs);
				axLo = max([ 0.0, min(omegaModelAvgIPts) - 0.5*(fModelDat.omega - min(omegaModelAvgIPts)) ]);
				axHi = min([ 2.0*fModelDat.omega, fModelDat.omega + 0.5*(fModelDat.omega - min(omegaModelAvgIPts)) ]);
				cap_omegaModelAvgIPts = cap( omegaModelAvgIPts, 0.0, axHi );
				cap_omegaModelAvgPPts = cap( omegaModelAvgPPts, 0.0, axHi );
				cap_omegaModelPVarIPts = cap( omegaModelPVarIPts, 0.0, axHi );
				cap_omegaModelPVarPPts = cap( omegaModelPVarPPts, 0.0, axHi );
				cap_omegaActualIPts = cap( omegaActualIPts, 0.0, axHi );
				cap_omegaActualPPts = cap( omegaActualPPts, 0.0, axHi );
				plot( ...
				  sPts, cap_omegaModelAvgIPts, 'o-', 'markersize', 30, ...
				  sPts, cap_omegaModelAvgPPts, 'x-', 'markersize', 26, ...
				  sPts, cap_omegaModelPVarIPts, 's-', 'markersize', 22, ...
				  sPts, cap_omegaModelPVarPPts, '+-', 'markersize', 18, ...
				  sPts, cap_omegaActualIPts, 'p-', 'markersize', 14, ...
				  sPts, cap_omegaActualPPts, '*-', 'markersize', 10, ...
				  sPts(indexOfPB)*[1.0,1.0], [axLo,axHi], 'r^-' );
				axis([ 0.0, 1.0, axLo, axHi ]);
				grid on;
				xlabel( "s" );
				ylabel( "omega" );
				legend( ...
				  "omega model avg I", ...
				  "omega model avg P", ...
				  "omega model pvar I", ...
				  "omega model pvar P", ...
				  "omega actual I", ...
				  "omega actual P", ...
				  "boundary for P", ...
				  "location", "northeast" );
				%
				numFigs++; figure(numFigs);
				plot( ...
				  sPts, bIPts, 'o-', ...
				  sPts, bPPts, 'x-', ...
				  sPts, 1.0+(0.0*sPts), 'k-', ...
				  sPts, max(abs(vecY'*fModelDat.matB)), 'rp-', ...
				  sPts(indexOfPB)*[1.0,1.0], [0.0,2.0], 'r^-' );
				axis([ 0.0, 1.0, 0.0, 3.0 ]);
				grid on;
				xlabel( "s" );
				ylabel( "b" );
				legend( ...
				  "bI", ...
				  "bP", ....
				  "1.0", ...
				  "trial step", ...
				  "boundary for P", ...
				  "location", "northwest" );
				%
				numFigs++; figure(numFigs);
				plot( ...
				  sPts, lIPts, 'o-', ...
				  sPts, lPPts, 'x-' );
				grid on;
				xlabel( "s" );
				ylabel( "||VLocal'*V*y||/||y||" );
			endif
			endif
			endif
			%%%
			%
			%
			vecX_trial = fModelDat.vecXPB;
			vecF_trial = funchF( vecX_trial );
			fevalCount++;
			omega_trial = sumsq(vecF_trial)/2.0;
			%msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  omega_trial = %g.", omega_trial ) );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  >> trial... omega: %8.2e to %8.2e, expected %8.2e ~ %8.2e;", ...
			  fModelDat.omega, omega_trial, ...
			  fModelDat.omegaModelAvgPB, fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB ) );
			msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "              fall: %8.2e, expected %8.2e ~ %8.2e.", ...
			  fModelDat.omega-omega_trial, ...
			  fModelDat.omega - fModelDat.omegaModelAvgPB, fModelDat.omega - (fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB) ) );
			if ( norm(vecF_trial) < norm(vecF_best) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Step is new best." );
				vecX_best = vecX_trial;
				vecF_best = vecF_trial;
			endif
			if ( ~isempty(vecF_cand) )
			if ( norm(vecF_trial) >= norm(vecF_cand) )
				msgif( prm.msgNotice, __FILE__, __LINE__, "Current trial is worse than earlier candidate; forcing acceptance of earlier candidate." );
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_cand)/2.0, fModelDat.omega-sumsq(vecF_cand)/2.0, fevalCount ) );
				fModelDat = __moveTo( vecX_cand, vecF_cand, fModelDat, prm );
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			endif
			%
			avefaThresh = mygetfield( prm, "avefaThresh", 0.5 ); % Actual vs expect fall acceptace threshold
			assert( 0.0 < avefaThresh );
			assert( avefaThresh < 1.0 );
			%msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  expected fall = %g.", fModelDat.omega - (fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB) ) );
			%msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( "  actual fall = %g.", fModelDat.omega - omega_trial ) );
			if ( omega_trial <= fModelDat.omega - avefaThresh*( fModelDat.omega - (fModelDat.omegaModelAvgPB + fModelDat.omegaModelVarPB) ) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Accepting step." );
				excellentThresh = mygetfield( prm, "excellentThresh", 0.1 );
				if ( norm(vecF_trial-fModelDat.vecFModelPB) <= excellentThresh*norm(fModelDat.vecFModelPB) )
					msgif( prm.msgCopious, __FILE__, __LINE__, "  Model was very accurate; removing wall(s)." );
					bRemoveFactor = mygetfield( prm, "bRemoveFactor", 1.5 );
					fModelDat = __removeB( bRemoveFactor*fModelDat.vecYPB, fModelDat, prm );
				endif
				msgif( prm.msgProgress, __FILE__, __LINE__, sprintf( "Moving from %10.3e to %10.3e (fall = %10.3e, fevalCount = %d).", ...
				  fModelDat.omega, sumsq(vecF_trial)/2.0, fModelDat.omega-sumsq(vecF_trial)/2.0, fevalCount ) );
				%
				%
				%%%
				fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm );
				%msg( __FILE__, __LINE__, "HACK! Reinitializing fModelDat!" );
				%[ fModelDat, datOut_initModel ] = __initModel( funchF, vecX_trial, vecF_trial, prm );
				%fevalCount += datOut_initModel.fevalCount;
				%%%
				%
				%
				stepCount++;
				datOut.iterCountOfSteps(stepCount+1) = iterCount+1;
				clear vecX_trial;
				clear vecF_trial;
				clear omega_trial;
				vecX_cand = [];
				vecF_cand = [];
				continue;
			endif
			msgif( prm.msgCopious, __FILE__, __LINE__, "  Rejecting step." );
			%
			if ( norm(vecF_trial) < norm(fModelDat.vecF) )
				% Trial is better than current, at least, so it's a candidate.
				vecX_cand = vecX_trial;
				vecF_cand = vecF_trial;
			endif
			%
			%
			% Consider adding a new barrier/wall.
			% But, first, make sure uncertainty in W was not the issue.
			vecY = fModelDat.vecYPB;
			vecU = fModelDat.matV*vecY;
			vecV = __calcOrthonorm( vecU, fModelDat.matVLocal, prm );
			refreshedTrialStep = false;
			if ( 0.0 ~= norm(vecV) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Refreshing trial step before adding wall." );
				%assert( reldiff( norm(fModelDat.matVLocal'*vecV), 0.0, eps ) < sqrt(eps) );
				%assert( reldiff( norm(fModelDat.matVLocal'*vecV), 0.0, eps ) < sqrt(eps) );
				%assert( 0.0 ~= norm(vecV) );
				%assert( fModelDat.lPB < 1.0-sqrt(eps) ); % Conceptually equiv to 0=norm(vecV), but less reliable?
				[ fModelDat, datOut_refresh ] = __refresh( vecY, funchF, fModelDat, prm );
				fevalCount += datOut_refresh.fevalCount;
				refreshedTrialStep = true;
			endif
			vecFModel = fModelDat.vecF + fModelDat.matW*vecY;
			badThresh = mygetfield( prm, "badThresh", 0.5 );
			if ( norm( vecF_trial - vecFModel ) > badThresh * norm(vecFModel) )
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Model was very bad; adding wall." );
				bAddFactor = mygetfield( prm, "bAddFactor", 0.5 );
				fModelDat = __addB( bAddFactor*vecY, fModelDat, prm );
			elseif (~refreshedTrialStep)
				msgif( prm.msgCopious, __FILE__, __LINE__, "  Did not refresh trial step but model was bed enough; adding wall." );
				bAddFactor = mygetfield( prm, "bAddFactor", 0.5 );
				fModelDat = __addB( bAddFactor*vecY, fModelDat, prm );
			endif
			%
			%
			%%%
			if (1)
			if ( makePlotsAndHalt )
				msg( __FILE__, __LINE__, "Making plots after refreshing and/or adding wall!" );
				numFigs = 10;
				%
				%
				%
				numPts = 101;
				sPts = linspace(0.999,0.0,numPts);
				sPts = ( 1.0 - sPts.^2 ).^2;
				sizeV = size(fModelDat.matV,2);
				sizeF = size(fModelDat.vecF,1);
				vecMG = -(fModelDat.matW'*fModelDat.vecF);
				matH0 = fModelDat.matW'*fModelDat.matW;
				matH1 = matH0 + fModelDat.matA;
				matSCurve = eye(sizeV,sizeV);
				%
				vecYIPts = zeros(sizeV,numPts);
				vecYPPts = zeros(sizeV,numPts);
				vecFModelIPts = zeros(sizeF,numPts);
				vecFModelPPts = zeros(sizeF,numPts);
				vecFActualIPts = zeros(sizeF,numPts);
				vecFActualPPts = zeros(sizeF,numPts);
				omegaModelAvgIPts = zeros(1,numPts);
				omegaModelAvgPPts = zeros(1,numPts);
				omegaModelPVarIPts = zeros(1,numPts);
				omegaModelPVarPPts = zeros(1,numPts);
				omegaActualIPts = zeros(1,numPts);
				omegaActualPPts = zeros(1,numPts);
				rhoIPts = zeros(1,numPts);
				rhoPPts = zeros(1,numPts);
				bIPts = zeros(1,numPts);
				bPPts = zeros(1,numPts);
				for n=1:numPts
					vecYIPts(:,n) = ( sPts(n)*matH0 + (1.0-sPts(n))*matSCurve ) \ (sPts(n)*vecMG);
					vecYPPts(:,n) = ( sPts(n)*matH1 + (1.0-sPts(n))*matSCurve ) \ (sPts(n)*vecMG);
					vecFModelIPts(:,n) = fModelDat.vecF + fModelDat.matW*vecYIPts(:,n);
					vecFModelPPts(:,n) = fModelDat.vecF + fModelDat.matW*vecYPPts(:,n);
					omegaModelAvgIPts(n) = sumsq(vecFModelIPts(:,n))/2.0;
					omegaModelAvgPPts(n) = sumsq(vecFModelPPts(:,n))/2.0;
					omegaModelPVarIPts(n) = omegaModelAvgIPts(:,n) + 0.5*vecYIPts(:,n)'*fModelDat.matA*vecYIPts(:,n);
					omegaModelPVarPPts(n) = omegaModelAvgPPts(:,n) + 0.5*vecYPPts(:,n)'*fModelDat.matA*vecYPPts(:,n);
					vecFActualIPts(:,n) = funchF( fModelDat.vecX + fModelDat.matV*vecYIPts(:,n) );
					vecFActualPPts(:,n) = funchF( fModelDat.vecX + fModelDat.matV*vecYPPts(:,n) );
					omegaActualIPts(n) = sumsq(vecFActualIPts(:,n))/2.0;
					omegaActualPPts(n) = sumsq(vecFActualPPts(:,n))/2.0;
					rhoIPts(n) = sumsq( vecFActualIPts(:,n) - vecFModelIPts(:,n) )/2.0;
					rhoPPts(n) = sumsq( vecFActualPPts(:,n) - vecFModelPPts(:,n) )/2.0;
					bIPts(n) = max(abs(vecYIPts(:,n)'*fModelDat.matB));
					bPPts(n) = max(abs(vecYPPts(:,n)'*fModelDat.matB));
					%lIPts(n) = (eps+norm(fModelDat.matVLocal'*fModelDat.matV*vecYIPts(:,n)))/(eps+norm(vecYIPts(:,n)));
					%lPPts(n) = (eps+norm(fModelDat.matVLocal'*fModelDat.matV*vecYPPts(:,n)))/(eps+norm(vecYPPts(:,n)));
					lIPts(n) = norm(fModelDat.matVLocal'*fModelDat.matV*vecYIPts(:,n))/norm(vecYIPts(:,n));
					lPPts(n) = norm(fModelDat.matVLocal'*fModelDat.matV*vecYPPts(:,n))/norm(vecYPPts(:,n));
				endfor
				%
				%
				indexOfPB = numPts;
				while (bPPts(indexOfPB)>1.0)
					indexOfPB--;
				endwhile
				%
				numFigs++; figure(numFigs);
				axLo = max([ 0.0, min(omegaModelAvgIPts) - 0.5*(fModelDat.omega - min(omegaModelAvgIPts)) ]);
				axHi = min([ 2.0*fModelDat.omega, fModelDat.omega + 0.5*(fModelDat.omega - min(omegaModelAvgIPts)) ]);
				cap_omegaModelAvgIPts = cap( omegaModelAvgIPts, 0.0, axHi );
				cap_omegaModelAvgPPts = cap( omegaModelAvgPPts, 0.0, axHi );
				cap_omegaModelPVarIPts = cap( omegaModelPVarIPts, 0.0, axHi );
				cap_omegaModelPVarPPts = cap( omegaModelPVarPPts, 0.0, axHi );
				cap_omegaActualIPts = cap( omegaActualIPts, 0.0, axHi );
				cap_omegaActualPPts = cap( omegaActualPPts, 0.0, axHi );
				plot( ...
				  sPts, cap_omegaModelAvgIPts, 'o-', 'markersize', 30, ...
				  sPts, cap_omegaModelAvgPPts, 'x-', 'markersize', 26, ...
				  sPts, cap_omegaModelPVarIPts, 's-', 'markersize', 22, ...
				  sPts, cap_omegaModelPVarPPts, '+-', 'markersize', 18, ...
				  sPts, cap_omegaActualIPts, 'p-', 'markersize', 14, ...
				  sPts, cap_omegaActualPPts, '*-', 'markersize', 10, ...
				  sPts(indexOfPB)*[1.0,1.0], [axLo,axHi], 'r^-' );
				axis([ 0.0, 1.0, axLo, axHi ]);
				grid on;
				xlabel( "s" );
				ylabel( "omega" );
				legend( ...
				  "omega model avg I", ...
				  "omega model avg P", ...
				  "omega model pvar I", ...
				  "omega model pvar P", ...
				  "omega actual I", ...
				  "omega actual P", ...
				  "boundary for P", ...
				  "location", "northeast" );
				%
				numFigs++; figure(numFigs);
				plot( ...
				  sPts, bIPts, 'o-', ...
				  sPts, bPPts, 'x-', ...
				  sPts, 1.0+(0.0*sPts), 'k-', ...
				  sPts, max(abs(vecY'*fModelDat.matB)), 'rp-', ...
				  sPts(indexOfPB)*[1.0,1.0], [0.0,2.0], 'r^-' );
				axis([ 0.0, 1.0, 0.0, 3.0 ]);
				grid on;
				xlabel( "s" );
				ylabel( "b" );
				legend( ...
				  "bI", ...
				  "bP", ....
				  "1.0", ...
				  "trial step", ...
				  "boundary for P", ...
				  "location", "northwest" );
				%
				numFigs++; figure(numFigs);
				plot( ...
				  sPts, lIPts, 'o-', ...
				  sPts, lPPts, 'x-' );
				grid on;
				xlabel( "s" );
				ylabel( "||VLocal'*V*y||/||y||" );
				%
				%
				%
				if (0)
				numPts = 101;
				pPts = linspace(0.0,1.0,numPts).^2;
				fModelNormPts = zeros(1,numPts);
				fModelResNormPts = zeros(1,numPts);
				for n=1:numPts
					omegaModelAvgPts(n) = sumsq( fModelDat.vecF + fModelDat.matW*vecY*pPts(n) )/2.0;
					omegaModelPVarPts(n) = omegaModelAvgPts(n) + (pPts(n)^2)*(vecY'*fModelDat.matA*vecY)/2.0;
					%omegaModelPVar0Pts(n) = omegaModelAvgPts(n) + (pPts(n)^2)*(vecY'*fModelDat.matA0*vecY)/2.0;
					rhoVsTrialPts(n) = sumsq( fModelDat.vecF + fModelDat.matW*vecY*pPts(n) - vecF_trial )/2.0;
					vecXAtPts(:,n) = fModelDat.vecX + fModelDat.matV*vecY*pPts(n);
					vecFAtPts(:,n) = funchF( vecXAtPts(:,n) );
					omegaAtPts(n) = sumsq( vecFAtPts(:,n) )/2.0;
					rhoAtPts(n) = sumsq( fModelDat.vecF + fModelDat.matW*vecY*pPts(n) - vecFAtPts(:,n) )/2.0;
				endfor
				%rd = reldiff( vecF_trial, vecFAtPts(:,numPts) )
				%assert( rd < sqrt(eps) );
				%
				%
				numFigs++; figure(numFigs);
				plot( ...
				  pPts, omegaModelAvgPts, 'o-', 'markersize', 20, ...
				  pPts, omegaModelPVarPts, 's-', 'markersize', 12, ...
				  pPts, omegaAtPts, 'p-', 'markersize', 4 );
				grid on;
				xlabel("p");
				legend( ...
				  "omega model avg", ...
				  "omega model pvar", ...
				  "omega actual", ...
				  "location", "northeast" );
				%
				%
				numFigs++; figure(numFigs);
				plot( ...
				  pPts, omegaAtPts, 'p-' );
				grid on;
				xlabel("p");
				legend( ...
				  "omega actual", ...
				  "location", "northeast" );
				%
				%
				numFigs++; figure(numFigs);
				plot( ...
				  pPts, rhoVsTrialPts, 'x-', 'markersize', 20, ...
				  pPts, rhoAtPts, '+-', 'markersize', 12 );
				grid on;
				xlabel("p");
				legend( ...
				  "rho vs trial", ...
				  "rho actual", ...
				  "location", "northeast" );
				%
				endif
			endif
			endif
			%%%
			%
			%
			%
			clear vecX_trial;
			clear vecF_trial;
			clear omega_trial;
			continue;
		endif
		endif
		%
		msgif( prm.msgCopious, __FILE__, __LINE__, "Refreshing subspace." );
		%%%
		%%% STUDY ME!
		%%%[ fModelDat, datOut_refresh ] = __refresh( fModelDat.vecYPB, funchF, fModelDat, prm );
		[ fModelDat, datOut_refresh ] = __refresh( fModelDat.vecYIB, funchF, fModelDat, prm );
		%%%
		%%%
		fevalCount += datOut_refresh.fevalCount;
		continue;
	endwhile
	%
	%
	datOut.iterCountVals(iterCount+1) = iterCount;
	datOut.fevalCountVals(iterCount+1) = fevalCount;
	datOut.fNormVals(iterCount+1) = norm(vecF_best);
	datOut.vecXVals(:,iterCount+1) = fModelDat.vecX;
	datOut.vecFVals(:,iterCount+1) = fModelDat.vecF;
	datOut.iterCountOfSteps(stepCount+1) = iterCount;
	%
	datOut.fevalCount = fevalCount;
	datOut.iterCount = iterCount;
	datOut.stepCount = stepCount;
return;
endfunction


function prm = __initPrm( vecX, vecF, prm )
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	prm.msgCopious = mygetfield( prm, "msgCopious", verbLev >= VERBLEV__COPIOUS );
	prm.msgProgress = mygetfield( prm, "msgProgress", verbLev >= VERBLEV__PROGRESS );
	prm.msgMain = mygetfield( prm, "msgMain", verbLev >= VERBLEV__MAIN );
	prm.msgNotice = mygetfield( prm, "msgNotice", verbLev >= VERBLEV__NOTICE );
	prm.debugMode = mygetfield( prm, "debugMode", true );
	%
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	%fTol = max([ sqrt(eps)*norm(vecF), eps ]);
	fTol = 100.0*eps;
	prm.fTol = mygetfield( prm, "fTol", fTol );
	prm.iterMax = mygetfield( prm, "iterMax", ceil(20 + 10*sqrt(sizeX) + 0.01*sizeX) );
	%
	prm.omegaTol = fTol^2/2.0;
return;
endfunction


function [ fModelDat, datOut ] = __initModel( funchF, vecX, vecF, prm )
	datOut = [];
	fevalCount = 0;
	%
	fModelDat = mygetfield( prm, "fModelDat_initial", [] );
	if (~isempty(fModelDat))
		assert( ~isempty(fModelDat,"vecX",[]) );
		assert( ~isempty(fModelDat,"vecF",[]) );
		assert( reldiff(fModelDat.vecX,vecX,eps^2) <= eps );
		if ( isempty(vecF) )
			vecF = fModelDat.vecF;
		else
			assert( reldiff(fModelDat.vecF,vecF,eps^2) <= eps );
		endif
		%
		% Do more checks here?
		%
		return;
	endif
	%
	if (isempty(vecF))
		vecF = funchF(vecX);
		fevalCount++;
	endif
	if (prm.debugMode)
		sizeX = size(vecX,1);
		sizeF = size(vecF,1);
		assert( sizeX >= 1 );
		assert( sizeF >= 1 );
		assert( isrealarray(vecX,[sizeX,1]) );
		assert( isrealarray(vecF,[sizeF,1]) );
	endif
	%
	%
	fNorm = norm(vecF);
	if ( 0.0 == fNorm )
		error( "Initial vecF is zero." );
	endif
	vecV = vecF/fNorm;
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecW) )
		error( "Initial vecW is zero." );
	endif
	%
	fModelDat.vecX = vecX;
	fModelDat.vecF = vecF;
	fModelDat.matVLocal = [ vecV ];
	fModelDat.matV = [ vecV ];
	fModelDat.matW = [ vecW ];
	fModelDat.matWLocal = [ vecW ];
	fModelDat.matA = [ 0.0 ];
	fModelDat.matA0 = [ 0.0 ];
	fModelDat.matB = [];
	%
	datOut.fevalCount = fevalCount;
return;
endfunction



function vecW = __calcJV( vecV, funchF, vecX, vecF, prm )
	v = norm(vecV);
	assert( 0.0 < v );
	epsFD = mygetfield( prm, "epsFD", eps^0.3 );
	assert( 0.0 < epsFD );
	vecFP = funchF( vecX + epsFD*vecV );
	vecW = ( vecFP - vecF ) / (epsFD*norm(vecV));
return;
endfunction


function fModelDat = __analyzeModel( fModelDat, prm )
	% Unpack.
	%matVLocal = fModelDat.matVLocal;
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	% Basics.
	matWTW = matW'*matW;
	vecMG = -(matW'*vecF);
	matSCurve = eye(sizeV,sizeV);
	%
	vecYIU = __calcStep( matWTW, vecMG, prm );
	vecYIB = __calcBoundStep( matWTW, vecMG, matB, matSCurve, prm );
	vecYPB = __calcBoundStep( matWTW + matA, vecMG, matB, matSCurve, prm );
	%
	vecFModelIU = vecF + (matW*vecYIU);
	vecFModelIB = vecF + (matW*vecYIB);
	vecFModelPB = vecF + (matW*vecYPB);
	%
	fModelDat.vecYIU = vecYIU;
	fModelDat.vecYIB = vecYIB;
	fModelDat.vecYPB = vecYPB;
	%
	fModelDat.vecXIU = vecX + (matV*vecYIU);
	fModelDat.vecXIB = vecX + (matV*vecYIB);
	fModelDat.vecXPB = vecX + (matV*vecYPB);
	%
	fModelDat.vecFModelIU = vecFModelIU;
	fModelDat.vecFModelIB = vecFModelIB;
	fModelDat.vecFModelPB = vecFModelPB;
	%
	fModelDat.omegaModelAvgIU = sumsq(vecFModelIU)/2.0;
	fModelDat.omegaModelAvgIB = sumsq(vecFModelIB)/2.0;
	fModelDat.omegaModelAvgPB = sumsq(vecFModelPB)/2.0;
	%
	fModelDat.omegaModelVarIU = vecYIU'*matA*vecYIU;
	fModelDat.omegaModelVarIB = vecYIB'*matA*vecYIB;
	fModelDat.omegaModelVarPB = vecYPB'*matA*vecYPB;
	%
	if ( isempty(matB) )
		fModelDat.bIU = 0.0;
		fModelDat.bIB = 0.0;
		fModelDat.bPB = 0.0;
	else
		fModelDat.bIU = max(abs(vecYIU'*matB));
		fModelDat.bIB = max(abs(vecYIB'*matB));
		fModelDat.bPB = max(abs(vecYPB'*matB));
	endif
	if (isempty(matVLocal))
		fModelDat.lIU = 0.0;
		fModelDat.lIB = 0.0;
		fModelDat.lPB = 0.0;
	else
		fModelDat.lIU = norm(matVLocal'*matV*vecYIU)/(eps+norm(vecYIU));
		fModelDat.lIB = norm(matVLocal'*matV*vecYIB)/(eps+norm(vecYIB));
		fModelDat.lPB = norm(matVLocal'*matV*vecYPB)/(eps+norm(vecYPB));
	endif
	%
	fModelDat.omega = sumsq(vecF)/2.0;
return;
endfunction


function vecY = __calcStep( matH, vecMG, prm )
	[ matR, cholFlag ] = chol( matH );
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	if ( 0 == cholFlag && min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
		vecY = matR \ (matR'\vecMG);
	else
		msgif( prm.msgNotice, __FILE__, __LINE__, "Extrapolating step to singular point. (Perhaps should bail instead?)" );
		hScale = max(max(abs(matH)));
		assert( 0.0 < hScale );
		sizeV = size(matH,1);
		matIV = eye(sizeV,sizeV);
		epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
		matH1 = matH + ( 1.0 * epsRelRegu * hScale ) * matIV;
		matH2 = matH + ( 2.0 * epsRelRegu * hScale ) * matIV;
		matR1 = chol( matH1 );
		matR2 = chol( matH2 );
		vecY1 = matR1 \ (matR1'\vecMG);
		vecY2 = matR2 \ (matR2'\vecMG);
		vecY = (2.0*vecY1) - vecY2;
	endif
return;
endfunction


function vecY = __calcBoundStep( matH, vecMG, matB, matSCurve, prm );
	if ( isempty(matB) || 0.0==max(max(abs(matB))) )
		vecY = __calcStep( matH, vecMG, prm );
		return;
	endif
	%
	%cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	%
	%[ matR, cholFlag ] = chol( matH );
	%if ( 0 == cholFlag && min(diag(matR)) > cholSafeTol*max(abs(diag(matR))) )
	% Looks like chol isn't accurate enough, and/or "\" triggers a check based on rcond()?
	% So, we'l use rcond too.
	rc = rcond( matH );
	if ( rcond(matH) > sqrt(eps) )
		s1 = 1.0;
	else
		epsRelRegu = mygetfield( prm, "epsRelRegu", sqrt(eps) );
		hScale = max(max(abs(matH)));
		sScale = max(max(abs(matSCurve)));
		assert( 0.0 < hScale );
		assert( 0.0 < sScale );
		s1 = 1.0 - (epsRelRegu*hScale/(epsRelRegu*hScale+sScale));
	endif
	assert( s1 >= 0.0 );
	%
	funchYOfS = @(s)( ( s*matH + (1.0-s)*matSCurve ) \ (s*vecMG) );
	funchBOfY = @(y)( max(abs(y'*matB)) );
	vecY1 = funchYOfS(s1);
	b1 = funchBOfY(vecY1);
	if ( b1 <= 1.0 )
		vecY = vecY1;
		return;
	endif
	%
	funchBM1OfS = @(s)( funchBOfY(funchYOfS(s)) - 1.0 );
	%
	s = fzero( funchBM1OfS, [0.0, s1] );
	vecY = funchYOfS(s);
	assert( reldiff(max(abs(vecY'*matB)),1.0) < sqrt(eps) );
return;
endfunction


function [ fModelDat, datOut ] = __expandModel( vecU, funchF, fModelDat, prm )
	fevalCount = 0;
	% Unpack.
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matA0 = fModelDat.matA0;
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,2);
	%
	%
	uNorm = norm(vecU);
	assert( 0.0 < uNorm );
	vecV = __calcOrthonorm( vecU, matV, prm );
	if ( 0.0 == norm(vecV) )
		error( "Failed to generate new subspace basis vector." );
	endif
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecV) )
		error( "Jacobian along new subspace basis vector is zero." );
	endif
	%
	%
	if (isempty(matVLocal))
		fModelDat.matVLocal = [ vecV ];
		fModelDat.matWLocal = [ vecW ];
	else
		fModelDat.matVLocal = [ matVLocal, vecV ];
		fModelDat.matWLocal = [ matWLocal, vecW ];
	endif
	%
	fModelDat.matV = [ matV, vecV ];
	fModelDat.matW = [ matW, vecW ];
	fModelDat.matA = zeros(sizeV+1,sizeV+1);
	fModelDat.matA(1:sizeV,1:sizeV) = matA;
	fModelDat.matB = [ matB; zeros(1,sizeB) ];
	fModelDat.matA0 = zeros(sizeV+1,sizeV+1);
	fModelDat.matA0(1:sizeV,1:sizeV) = matA0;
	%
	%
	datOut.fevalCount = fevalCount;
return;
endfunction



function vecV = __calcOrthonorm( vecU, matV, prm )
	numPasses = 2;
	u0 = norm(vecU);
	if (0.0==u0)
		vecV = 0.0*vecU;
		return;
	endif
	if (isempty(matV))
		vecV = vecU/u0;
		return;
	endif
	orthoTol = mygetfield( prm, "orthoTol", 1.0e-10 );
	vecV = vecU;
	for n=1:numPasses
		vecV -= matV*(matV'*vecV);
		v = norm(vecV);
		if ( v <= orthoTol*u0 )
			vecV(:) = 0.0;
			return;
		else
			vecV /= v;
		endif
	endfor
return;
endfunction


function fModelDat = __moveTo( vecX_trial, vecF_trial, fModelDat, prm )
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matA0 = fModelDat.matA0;
	%matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	msgif( prm.msgCopious, __FILE__, __LINE__, sprintf( " ( ||F||: %10.3e -> %10.3e. )", norm(vecF), norm(vecF_trial) ) );
	%
	%
	vecDeltaX = vecX_trial - vecX;
	vecY = matV'*vecDeltaX;
	yNorm = norm(vecY);
	assert( 0.0 < yNorm );
	%
	%
	%
	msgif( prm.msgCopious, __FILE__, __LINE__, "TODO: consider partial quadratic update (OSQU)." );
	msgif( prm.msgCopious, __FILE__, __LINE__, "  Here, only linear (Broyden) update is applied." );
	vecFModel_trial = vecF + matW*vecY;
	vecRho = vecF_trial - vecFModel_trial;
	matW_plus = matW + vecRho * (vecY')/(yNorm^2);
	%
	%
	%
	%%%
	%%% STUDY ME!
	vecYHat = vecY/yNorm;
	if (0)
	stepUpdateAccuracyCoeff = mygetfield( prm, "stepUpdateAccuracyCoeff", 0.0 );
	assert( 0.0 <= stepUpdateAccuracyCoeff );
	assert( stepUpdateAccuracyCoeff <= 1.0 );
	matE = eye(sizeV,sizeV) - (stepUpdateAccuracyCoeff*vecYHat)*(vecYHat');
	%
	matD = diag(max([ abs(diag(matW'*matW)), abs(diag(matW_plus'*matW_plus)) ]'));
	foo1 = sumsq( matW_plus*vecY - matW*vecY ) - vecY'*matA*vecY;
	foo2 = vecY'*matD*vecY;
	if ( foo1 <= 0.0 )
		s = 0.0;
	elseif ( foo2 <= foo1 )
		s = 1.0;
	else
		s = foo1 / foo2;
		assert( 0.0 <= s );
		assert( s <= 1.0 );
	endif
	coeffD = mygetfield( prm, "coeffD", 10.0 );
	s*=coeffD;
	matA0 = matE*( matA + s*matD )*matE;
	endif
	%%%matA0 = matA + s*matD;
	%
	matWTW = matW'*matW;
	matA0 = 100.0*( diag(diag(matWTW)) + eye(sizeV,sizeV)*max(max(matWTW))*0.1 );
	%%%matA0 += matA + 0.00001*diag(diag(matWTW));
	%%%
	%%%
	%
	if (0)
		msg( __FILE__, __LINE__, "Data dump..." );
		vecD = diag(matWTW);
		matD = diag(vecD);
		dAvg = sum(vecD)/sizeV;
		dSqAvg = sum(vecD.^2)/sizeV;
		dVar = sqrt( dSqAvg - dAvg^2 );
		dMin = min(vecD);
		dMax = max(vecD);
		[ sumsq(vecRho), dAvg, dVar, dMin, dMax ]
		[ sumsq(vecRho)/sumsq(vecY), sumsq(vecRho)*dAvg/sumsq(matW*vecY), sumsq(vecRho)*dAvg/(vecY'*matD*vecY) ]
		error( "HALT!" );
	endif
	%
	%
	%
	fModelDat.matVLocal = [];
	fModelDat.matWLocal = [];
	fModelDat.vecX = vecX_trial;
	fModelDat.vecF = vecF_trial;
	fModelDat.matW = matW_plus;
	fModelDat.matA = matA0;
	fModelDat.matA0 = matA0;
return;
endfunction


function [ fModelDat, datOut ] = __refresh( vecY, funchF, fModelDat, prm )
	fevalCount = 0;
	% Unpack.
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matWLocal = fModelDat.matWLocal;
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA0 = fModelDat.matA0; % Hessian variation matrix, < (delta W)' * (delta W) >.
	%matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeVLocal = size(matVLocal,2);
	%sizeB = size(matB,2);
	%
	%
	yNorm = norm(vecY);
	assert( 0.0 < yNorm );
	vecYHat = vecY/yNorm;
	vecU = matV*vecYHat;
	vecV = __calcOrthonorm( vecU, matVLocal, prm );
	%%
	doSaveMeHack20220511 = true;
	if (doSaveMeHack20220511)
		clear vecYHat;
		if ( 0.0 == norm(vecV) )
			vecU = matV*fModelDat.vecYPB;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		if ( 0.0 == norm(vecV) )
			vecU = matV*fModelDat.vecYIB;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		if ( 0.0 == norm(vecV) )
			vecU = matV*fModelDat.vecYIU;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		% Expand model?
		if ( 0.0 == norm(vecV) )
			vecU = matW'*fModelDat.vecF;
			vecV = __calcOrthonorm( vecU, matVLocal, prm );
		endif
		% Try a random vector???
		if ( 0.0 == norm(vecV) )
			error( "Failed to generate local subspace basis vector." );
		endif
	endif
	%%
	if ( 0.0 == norm(vecV) )
		error( "Failed to generate local subspace basis vector." );
	endif
	if (prm.debugMode)
		prmTemp.orthoTol = sqrt(eps);
		foo = __calcOrthonorm( vecV, matV, prmTemp );
		if ( norm(foo) != 0.0 )
			msg( __FILE__, __LINE__, "Data dump..." );
			[ norm(vecU), norm(vecV), norm(foo) ]
			[ norm(vecU), norm(vecV), norm(foo) ] - 1.0
			[ vecV'*vecU, foo'*vecU, foo'*vecV ]
			[ norm( matV'*vecU ), norm( matV'*vecV ), norm( matV'*foo ) ]
			[ norm( vecV - matV*(matV'*vecV) ), norm( vecU - matV*(matV'*vecU) ), norm( foo - matV*(matV'*foo) ) ]
			[ norm( vecV - matV*(matV'*vecV) ), norm( vecU - matV*(matV'*vecU) ), norm( foo - matV*(matV'*foo) ) ] - 1.0
			error( "Subspace basis vector was thrown out of own space?!?!" );
		endif
	endif
	%
	vecW = __calcJV( vecV, funchF, vecX, vecF, prm );
	fevalCount++;
	if ( 0.0 == norm(vecW) )
		error( "Jacobian along local subspace basis vector is zero." );
	endif
	%
	matVLocal = [ matVLocal, vecV ];
	matWLocal = [ matWLocal, vecW ];
	sizeVLocal++;
	if ( prm.debugMode )
		assert( reldiff( matVLocal, matV*(matV'*matVLocal), eps ) < sqrt(eps) );
		assert( reldiff( eye(sizeVLocal,sizeVLocal), matVLocal'*matVLocal, eps ) < sqrt(eps) );
	endif
	matVLocal = matV*(matV'*matVLocal);
	for n=1:sizeVLocal
		matVLocal(:,n) /= norm(matVLocal(:,n));
	endfor
	%
	matE = eye(sizeV,sizeV) - matV'*matVLocal*(matVLocal')*matV;
	%
	fModelDat.matVLocal = matVLocal;
	fModelDat.matWLocal = matWLocal;
	%
	%%% The yHats are NOT perpendicular!
	%%%fModelDat.matW = matW + (vecW - matW*vecYHat)*(vecYHat');
	%%% This z thing also works..
	%%%vecZ = matV'*vecV;
	%%%fModelDat.matW = matW + (vecW - matW*vecZ)*(vecZ');
	% But, to be safe, let's do this.
	fModelDat.matW = matW + ( matWLocal - matW*(matV')*matVLocal)*(matVLocal'*matV); % Yeah?
	%
	fModelDat.matA = matE * matA0 * matE;
	%
	%
	datOut.fevalCount = fevalCount;
return;
endfunction



function fModelDat = __addB( vecY, fModelDat, prm )
	% Unpack.
	%vecX = fModelDat.vecX;
	%vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%matV = fModelDat.matV; % Subspace basis matrix.
	%matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	%matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	%sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	ysumsq = sumsq(vecY);
	assert( 0.0 < ysumsq );
	if (isempty(matB))
		fModelDat.matB = vecY/ysumsq;
	else
		fModelDat.matB = [ matB, vecY/ysumsq ];
	endif
return;
endfunction



function fModelDat = __removeB( vecY, fModelDat, prm )
	% Unpack.
	%vecX = fModelDat.vecX;
	%vecF = fModelDat.vecF;
	%omega = fModelDat.omega;
	%matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	%matV = fModelDat.matV; % Subspace basis matrix.
	%matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	%matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	%sizeX = size(vecX,1);
	%sizeF = size(vecF,1);
	%sizeV = size(matV,2);
	%sizeB = size(matB,2);
	%
	if ( isempty(matB) )
		return;
	endif
	msk = logical( abs(vecY'*matB) >= 1.0 );
	fModelDat.matB = matB( :, msk );
return;
endfunction

function __dumpModel( fModelDat, prm )
	msg( __FILE__, __LINE__, "vvvvvvvvvvvvvvvvvvvv Begin __dumpModel()..." );
	vecX = fModelDat.vecX;
	vecF = fModelDat.vecF;
	omega = fModelDat.omega;
	matVLocal = fModelDat.matVLocal; % Locally evaluated subspace basis matrix.
	matV = fModelDat.matV; % Subspace basis matrix.
	matW = fModelDat.matW; % Projected subspace basis matrix, J*V.
	matA = fModelDat.matA; % Hessian variation matrix, < (delta W)' * (delta W) >.
	matB = fModelDat.matB; % Boundary / trust region matrix; steps must satify max(abs(y'*B)) <= 1.
	sizeX = size(vecX,1);
	sizeF = size(vecF,1);
	sizeV = size(matV,2);
	sizeB = size(matB,2);
	%
	echo__omega = omega
	sizeVLocal = size(matVLocal,2)
	sizeV = size(matV,2)
	%
	vecYIU = fModelDat.vecYIU;
	vecYIB = fModelDat.vecYIB;
	vecYPB = fModelDat.vecYPB;
	%
	vecFModelIU = vecF + (matW*vecYIU);
	vecFModelIB = vecF + (matW*vecYIB);
	vecFModelPB = vecF + (matW*vecYPB);
	%
	omegaModelAvgIU = sumsq(vecFModelIU)/2.0
	omegaModelAvgIB = sumsq(vecFModelIB)/2.0
	omegaModelAvgPB = sumsq(vecFModelPB)/2.0
	%
	omegaModelVarIU = vecYIU'*matA*vecYIU
	omegaModelVarIB = vecYIB'*matA*vecYIB
	omegaModelVarPB = vecYPB'*matA*vecYPB
	%
	if (isempty(matB))
		bIU = 0.0
		bIB = 0.0
		bPB = 0.0
	else
		bIU = max(abs(vecYIU'*matB))
		bIB = max(abs(vecYIB'*matB))
		bPB = max(abs(vecYPB'*matB))
	endif
	%
	msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^ End __dumpModel()." );
return;
endfunction
