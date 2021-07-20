function [ s, p, datOut, retCode ] = extFit_iterSolve( ...
  s0, p0, xVals, fVals, nFit, wVals=[], prm=[], datIn=[] )
	commondefs;
	thisFile = "extFit_iterSolve";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	datOut = [];
	%
	numPts = size(xVals,2);
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	doChecks = mygetfield( prm, "doChecks", true );
	if ( doChecks )
		assert( isrealscalar(s0) );
		assert( isrealscalar(p0) );
		assert( 0.0 < p0 );
		%
		assert( numPts >= 3 ); % Really, we should have 5+ points.
		assert( isrealarray(xVals,[1,numPts]) );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
		assert( isrealarray(fVals,[1,numPts]) );
		%
		assert( isrealscalar(nFit) );
		assert( fleq(nFit,round(nFit)) );
		assert( 1 <= nFit );
		assert( nFit <= numPts );
		%
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	prm_getOmega = mygetfield( prm, "prm_getOmega", [] );
	omega0 = extFit_getOmega( s0, p0, xVals, fVals, nFit, wVals, prm_getOmega );
	%
	iterLimit = mygetfield( prm, "iterLimit", 10 );
	omegaTol = mygetfield( prm, "omegaTol", eps*sum(fVals.^2) )
	deltaSTol = mygetfield( prm, "deltaSTol", sqrt(eps)*(max(xVals)-min(xVals)) );
	deltaPTol = mygetfield( prm, "deltaPTol", sqrt(sqrt(eps)) );
	epsS = mygetfield( prm, "epsS", (eps^0.75)*(max(xVals)-min(xVals)) );
	epsP = mygetfield( prm, "epsP", eps^0.75 );
	if ( doChecks )
		assert( isrealscalar(iterLimit) );
		assert( iterLimit >= 1 );
		assert( isrealscalar(omegaTol) );
		assert( omegaTol > 0.0 );
		assert( isrealscalar(deltaSTol) );
		assert( deltaSTol > 0.0 );
		assert( isrealscalar(deltaPTol) );
		assert( deltaPTol > 0.0 );
		assert( isrealscalar(epsS) );
		assert( epsS > 0.0 );
		assert( isrealscalar(epsP) );
		assert( epsP > 0.0 );
	end
	%
	%
	%
	s = s0;
	p = p0;
	omega = omega0;
	iterCount = 0;
	while (1)
		if ( omega <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged (%g <= %g).", omega, omegaTol ) );
			retCode = RETCODE__SUCCESS;
			return;
		elseif ( iterCount >= iterLimit )
			msg_notify( verbLev, thisFile, __LINE__, sprintf( ...
			  "Reached iterLimit (%d).", iterLimit ) );
			retCode = RETCODE__IMPOSED_STOP;
			return;
		end
		%
		rhoPPVals = extFit_getRho( s+epsS, p+epsP, xVals, fVals, nFit, prm_getOmega );
		rhoP0Vals = extFit_getRho( s+epsS, p,      xVals, fVals, nFit, prm_getOmega );
		rhoPMVals = extFit_getRho( s+epsS, p-epsP, xVals, fVals, nFit, prm_getOmega );
		rho0PVals = extFit_getRho( s,      p+epsP, xVals, fVals, nFit, prm_getOmega );
		rho00Vals = extFit_getRho( s,      p,      xVals, fVals, nFit, prm_getOmega );
		rho0MVals = extFit_getRho( s,      p-epsP, xVals, fVals, nFit, prm_getOmega );
		rhoMPVals = extFit_getRho( s-epsS, p+epsP, xVals, fVals, nFit, prm_getOmega );
		rhoM0Vals = extFit_getRho( s-epsS, p,      xVals, fVals, nFit, prm_getOmega );
		rhoMMVals = extFit_getRho( s-epsS, p-epsP, xVals, fVals, nFit, prm_getOmega );
		%
		rhoDSVals = ( rhoP0Vals - rhoM0Vals ) / ( 2.0 * epsS );
		rhoDPVals = ( rho0PVals - rho0MVals ) / ( 2.0 * epsP );
		%rhoDSDSVals = ( rhoP0Vals - 2.0*rho00Vals + rhoM0Vals ) / ( epsS^2 );
		%rhoDPDPVals = ( rho0PVals - 2.0*rho00Vals + rho0MVals ) / ( epsP^2 );
		%rhoDSDPVals = ( rhoPPVals + rhoMMVals - rhoPMVals - rhoMPVals ) / ( 4.0 * epsS * epsP );
		%
		sigma0S = sum( wVals .* rho00Vals .* rhoDSVals );
		sigma0P = sum( wVals .* rho00Vals .* rhoDPVals );
		sigmaSS = sum( wVals .* rhoDSVals .* rhoDSVals );
		sigmaPP = sum( wVals .* rhoDPVals .* rhoDPVals );
		sigmaSP = sum( wVals .* rhoDSVals .* rhoDPVals );
		%sigma0SS = sum( wVals .* rho00Vals .* rhoDSDSVals );
		%sigma0PP = sum( wVals .* rho00Vals .* rhoDPDPVals );
		%sigma0SP = sum( wVals .* rho00Vals .* rhoDSDPVals );
		%
		vecG = [ sigma0S; sigma0P ];
		matH = [ sigmaSS, sigmaSP; sigmaSP, sigmaPP ];
		%matH2 = matH1 + [ sigma0SS, sigma0SP; sigma0SP, sigma0PP ];
		matD = diag(abs(diag(matH))); % There are other possibilities.
		%
		for mu = [ 0.00, 0.01, 0.1, 1.0, 10.0, 1.0E2, 1.0E3, 1.0E4, 1.0E5 ]
			vecDelta = -( matH + mu*matD ) \ vecG;
			deltaS = vecDelta(1);
			deltaP = vecDelta(2);
			s_trial = s + deltaS;
			p_trial = p + deltaP;
			omega_trial = extFit_getOmega( s_trial, p_trial, xVals, fVals, nFit, prm_getOmega );
			if ( omega_trial < omega )
				break;
			end
			if (  (abs(deltaS) < deltaSTol)  &&  (abs(deltaP) < deltaSTol)  )
				msg_main( verbLev, thisFile, __LINE__, sprintf( ...
				  "Hit local extermum (|%g| <= %g, |%g| <= %g).", ...
				  deltaS, deltaSTol, deltaP, deltaPTol ) );
				retCode = RETCODE__ALGORITHM_BREAKDOWN;
				return;
			end
		end
		iterCount++;
		%
		msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "%3d;  %11.3e, %11.3e, %11.3e;  %11.3e, %11.3e", ...
		  iterCount, ...
		  s, ...
		  p, ...
		  omega, ...
		  s_trial - s, ...
		  p_trial - p ) );
		%
		s = s_trial;
		p = p_trial;
		omega = omega_trial;
	end
% This line should be unreachable.
end
