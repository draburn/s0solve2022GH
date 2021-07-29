function [ rhoVals, bigF0, bigF1, omega, vecG, matH ] = extFit__calcAboutPt( ...
  s, p, xVals, fVals, nExactFit, wVals=[], prm=[] )
	%
	thisFile = "extFit__calcAboutPt";
	doChecks = mygetfield( prm, "doChecks", true );
	%
	if ( isempty(wVals) )
		wVals = ones(size(xVals));
	end
	if ( doChecks )
		assert( isrealscalar(s) );
		assert( isrealscalar(p) );
		assert( 0.0 < p );
		numPts = size(xVals,2);
		assert( numPts >= 3 );
		assert( isrealarray(xVals,[1,numPts]) );
		xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
		assert( xValsAreStrictlyIncreasing );
		assert( isrealarray(fVals,[1,numPts]) );
		assert( isposintscalar(nExactFit) );
		assert( 1 <= nExactFit );
		assert( nExactFit <= numPts );
		assert( isrealarray(wVals,[1,numPts]) );
		noWValIsNegative = (0==sum(wVals<0.0));
		assert( noWValIsNegative );
		atLeastOneWValIsPositive = (1<=sum(wVals>0.0));
		assert( atLeastOneWValIsPositive );
	end
	%
	epsS = mygetfield( prm, "epsS", eps^0.25*(max(xVals)-min(xVals)) );
	epsP = mygetfield( prm, "epsP", eps^0.25 );
	if ( doChecks )
		assert( isrealscalar(epsS) );
		assert( epsS > 0.0 );
		assert( isrealscalar(epsP) );
		assert( epsP > 0.0 );
	end
	%
	prm_calcAtPt = mygetfield( prm, "prm_calcAtPt", [] );
	[ rhoVals, bigF0, bigF1, omega ] = extFit__calcAtPt( ...
	  s, p, xVals, fVals, nExactFit, wVals, prm_calcAtPt );
	rhoVals_p0 = extFit__calcAtPt( s+epsS, p, xVals, fVals, nExactFit, wVals, prm_calcAtPt );
	rhoVals_m0 = extFit__calcAtPt( s-epsS, p, xVals, fVals, nExactFit, wVals, prm_calcAtPt );
	rhoVals_0p = extFit__calcAtPt( s, p+epsP, xVals, fVals, nExactFit, wVals, prm_calcAtPt );
	rhoVals_0m = extFit__calcAtPt( s, p-epsP, xVals, fVals, nExactFit, wVals, prm_calcAtPt );
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
return;
end