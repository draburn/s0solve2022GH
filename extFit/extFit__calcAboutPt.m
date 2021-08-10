function [ rhoVals, bigF0, bigF1, omega, vecG, matH, retCode ] = extFit__calcAboutPt( ...
  s, p, xVals, fVals, wVals, prm=[] )
	%
	commondefs;
	thisFile = "extFit__calcAboutPt";
	doChecks = mygetfield( prm, "doChecks", false );
	%
	if ( doChecks )
		assert( isrealscalar(s) );
		assert( isrealscalar(p) );
	end
	%
	epsS = mygetfield( prm, "epsS", eps^0.25*(max(xVals)-min(xVals)) );
	epsP = mygetfield( prm, "epsP", eps^0.25 );
	numPts = size(xVals,2);
	if ( doChecks )
		assert( isrealscalar(epsS) );
		assert( epsS > 0.0 );
		assert( isrealscalar(epsP) );
		assert( epsP > 0.0 );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	if ( mygetfield(prm,"useTurbo",true) )
		dVals = sqrt(wVals);
		[ rhoVals, bigF0, bigF1, omega, flag_00 ] = extFit__calcAtPt_turbo( s, p, xVals, fVals, dVals );
		[ rhoVals_p0, flag_p0 ] = extFit__calcRhoVals( s+epsS, p, xVals, fVals, dVals );
		[ rhoVals_m0, flag_m0 ] = extFit__calcRhoVals( s-epsS, p, xVals, fVals, dVals );
		[ rhoVals_0p, flag_0p ] = extFit__calcRhoVals( s, p+epsP, xVals, fVals, dVals );
		[ rhoVals_0m, flag_0m ] = extFit__calcRhoVals( s, p-epsP, xVals, fVals, dVals );
		if ( flag_00 || flag_p0 || flag_m0 || flag_0m || flag_0p )
			retCode = RETCODE__BAD_INPUT;
			return;
		end
	else % NOT TURBO
	%
	prm_calcAtPt = mygetfield( prm, "prm_calcAtPt", [] );
	[ rhoVals, bigF0, bigF1, omega, retCode ] = extFit__calcAtPt(
	  s, p, xVals, fVals, wVals, prm_calcAtPt );
	if ( RETCODE__SUCCESS ~= retCode )
		return;
	end
	%
	[ rhoVals_p0, f0, f1, f2, retCode ] = extFit__calcAtPt( s+epsS, p, xVals, fVals, wVals, prm_calcAtPt );
	if ( RETCODE__SUCCESS ~= retCode )
		return;
	end
	[ rhoVals_m0, f0, f1, f2, retCode ] = extFit__calcAtPt( s-epsS, p, xVals, fVals, wVals, prm_calcAtPt );
	if ( RETCODE__SUCCESS ~= retCode )
		return;
	end
	[ rhoVals_0p, f0, f1, f2, retCode ] = extFit__calcAtPt( s, p+epsP, xVals, fVals, wVals, prm_calcAtPt );
	if ( RETCODE__SUCCESS ~= retCode )
		return;
	end
	[ rhoVals_0m, f0, f1, f2, retCode ] = extFit__calcAtPt( s, p-epsP, xVals, fVals, wVals, prm_calcAtPt );
	if ( RETCODE__SUCCESS ~= retCode )
		return;
	end
	end % NOT TURBO
	rhoDSVals = ( rhoVals_p0 - rhoVals_m0 ) / ( 2.0*epsS );
	rhoDPVals = ( rhoVals_0p - rhoVals_0m ) / ( 2.0*epsP );
	%
	sigma0S = sum( wVals .* rhoVals .* rhoDSVals );
	sigma0P = sum( wVals .* rhoVals .* rhoDPVals );
	sigmaSS = sum( wVals .* rhoDSVals .* rhoDSVals );
	sigmaPP = sum( wVals .* rhoDPVals .* rhoDPVals );
	sigmaSP = sum( wVals .* rhoDSVals .* rhoDPVals );
	%
	vecG = [ sigma0S; sigma0P ];
	matH = [ sigmaSS, sigmaSP; sigmaSP, sigmaPP ];
	%
	retCode = RETCODE__SUCCESS;
return;
end
