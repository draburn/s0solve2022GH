thisFile = "findGoodCand__fofprime__sub";
%
% Make sure we have enough points monotonic...
mCheckLo = n - numPtsToCheck;
mCheckHi = n - 1;
for m = mCheckLo : mCheckHi
if ( gVals(m) <= gVals(m+1) )
	msg_copious( verbLev, thisFile, __LINE__, "Insufficient points." );
	thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
	return;
end
end
assert( isrealarray(chVals(mCheckLo:mCheckHi),[1,numPtsToCheck]) );


% Build one-sided linear model.
% Later, could implement bounded quadratic, etc.
mFitLo = n - numPtsForFit;
mFitHi = n - 1;
cxFitLo = cxVals(mFitLo);
cxFitHi = cxVals(mFitHi);
%
% A bit excessive to recalc? 
cyVals = (cxVals - cxFitLo) / (cxFitHi - cxFitLo);
%msg_copious( verbLev, thisFile, __LINE__, "Looking at..." );

if 0&&( verbLev >= VERBLEV__COPIOUS )
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
if 0&&( verbLev >= VERBLEV__COPIOUS )
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
if 0&&( verbLev >= VERBLEV__COPIOUS )
	echo__sumSqVal = sumSqVal
	echo__sumSqRes = sumSqRes
end
if ( sumSqRes >= sumSqVal/10.0 )
	msg_copious( verbLev, thisFile, __LINE__, "Bad fit." );
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
msg_copious( verbLev, thisFile, __LINE__, ...
  sprintf("xTemp = %f.",xTemp) );


% Well-behaved?
if ( (xVals(n+1)-xTemp) * (xTemp-xVals(n)) > 0.0 )
	xCand = xTemp;
	meritCand = 1.0;
	msg_copious( verbLev, thisFile, __LINE__, "Accepting." );
	thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
	return;
end

% Not well behaved. So, simple hack...
if ( xVals(n+1) < xVals(n) )
	xCand = (xVals(n+1) + xVals(n))/2.0;
	meritCand = 1.0;
	msg_copious( verbLev, thisFile, __LINE__, "Forcing bisection." );
	thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
	return;
end
msg_copious( verbLev, thisFile, __LINE__, "OOB." );
thisFile = "RETURNING FROM findGoodCand__fofprime__sub";
return;
