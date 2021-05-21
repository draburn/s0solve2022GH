thisFile = "findGoodCand__fofprime__sub";
msg( thisFile, __LINE__, "Hey!" );
%
% Make sure we have enough points monotonic...
mCheckLo = n + 1 - numPtsToCheck;
mCheckHi = n;
for m = mCheckLo+1 : mCheckHi
if ( gVals(m-1) <= gVals(m) )
	thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
	return;
end
end
assert( isrealarray(chVals(mCheckLo:mCheckHi),[1,numPtsToCheck]) );


% Build linear model.
mFitLo = n + 1 - numPtsForFit;
mFitHi = n;
cxFitLo = cxVals(mFitLo);
cxFitHi = cxVals(mFitHi);
%
% A bit excessive to recalc? 
cyVals = (cxVals - cxFitLo) / (cxFitHi - cxFitLo);
msg_copious( verbLev, thisFile, __LINE__, "Looking at..." );

if ( verbLev >= VERBLEV__COPIOUS )
	echo__myFitXVals = cxVals(mFitLo:mFitHi)
	echo__myFitYVals = cyVals(mFitLo:mFitHi)
	echo__myFitHVals = chVals(mFitLo:mFitHi)
	echo__myCheckXVals = cxVals(mCheckLo:mCheckHi)
	echo__myCheckYVals = cyVals(mCheckLo:mCheckHi)
	echo__myCheckHVals = chVals(mCheckLo:mCheckHi)
end
%
cyVec = cyVals(mFitLo:mFitHi)';
chVec = chVals(mFitLo:mFitHi)';
assert( isrealarray(cyVec,[numPtsForFit,1]) );
assert( isrealarray(chVec,[numPtsForFit,1]) );
cyMat = [ ones(numPtsForFit,1), cyVec ];
coeffVec = cyMat \ chVec;
assert( isrealarray(coeffVec,[2,1]) );
c0 = coeffVec(1);
c1 = coeffVec(2);

% Also excessive, but, convenient.
chModelVals = c0 + c1*cyVals;
if ( verbLev >= VERBLEV__COPIOUS )
	echo__cyVec = cyVec
	echo__chVec = chVec
	echo__coeffVec = coeffVec
	echo__myCheckHModelVals = chModelVals(mCheckLo:mCheckHi)
end

% Validate fit.
sumSqRes = 0.0;
sumSqVal = 0.0;
for m = mCheckLo : mCheckHi
	hModel = c0 + c1 * cyVals(m);
	sumSqVal += chVals(m)^2;
	sumSqRes += ( chVals(m) - chModelVals(m) )^2;
end
if ( verbLev >= VERBLEV__COPIOUS )
	echo__sumSqVal = sumSqVal
	echo__sumSqRes = sumSqRes
end
if ( sumSqRes >= sumSqVal/10.0 )
	thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
	return;
end
if ( abs(c1) < eps*abs(c0) )
	msg_flagged( verbLev, thisFile, __LINE__, "abs(c1) < eps*abs(c0)!" );
	thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
	return;
end
yTemp = -c0 / c1;
xTemp = cxFitLo + yTemp*(cxFitHi-cxFitLo);

% Only do this once, for now.
distCheck = abs( xVals-xTemp );
if ( min(abs(xVals-xTemp)) > sqrt(eps)*(cxFitHi-cxFitLo) )
	xCand = xTemp;
	meritCand = 1.0;
	msg_copious( verbLev, thisFile, __LINE__, "Accepting." );
	msg_copious( verbLev, thisFile, __LINE__, "  xNew via 'f over f prime'." );
	thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
	return;
end
msg_copious( verbLev, thisFile, __LINE__, "Rejecting." );
