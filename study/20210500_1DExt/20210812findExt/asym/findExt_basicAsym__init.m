% Do standard stuff.
commondefs;
thisFile = "findExt_basicAsym__init";
verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
%
%
% Set default return values.
xExt = [];
fExt = [];
retCode = RETCODE__NOT_SET;
datOut = [];
%
%
% Unpack & validate primary data.
numPts = size(xVals,2);
if ( valLev >= VALLEV__LOW )
	assert( numPts >= 3 );
	assert( isrealarray(xVals,[1,numPts]) );
	assert( isrealarray(fVals,[1,numPts]) );
end
xMax = max(xVals);
xMin = min(xVals);
if ( valLev >= VALLEV__LOW )
	notMax = (  abs(xVals-xMax)  > eps075*(abs(xVals)+abs(xMax)) );
	notMin = (  abs(xVals-xMin)  > eps075*(abs(xVals)+abs(xMin)) );
	numValsNotMaxOrMin = sum(notMax .* notMin);
	haveAtLeast3UniqueXVals = ( numValsNotMaxOrMin >= 1 );
	assert( haveAtLeast3UniqueXVals );
	clear numValsNotMaxOrMin;
	clear notMin;
	clear notMax;
end
%
%
% Unpack & validate initial guess.
if (isempty(xExt_initial))
	xExt_initial = (xMax+xMin)/2.0;
end
if ( valLev >= VALLEV__LOW )
	assert( isrealscalar(xExt_initial) );
	assert( isrealscalar(pL_initial) );
	assert( isrealscalar(pR_initial) );
end
%
%
% Unpack & validate bounds.
xExtMin = mygetfield( prm, "xExtMin", xMin - 10.0*(xMax-xMin) );
xExtMax = mygetfield( prm, "xExtMax", xMax + 10.0*(xMax-xMin) );
if ( valLev >= VALLEV__LOW )
	assert( isrealorinfscalar(xExtMin) );
	assert( isrealorinfscalar(xExtMax) );
end
assert( xExtMin <= xExt_initial );
assert( xExtMax >= xExt_initial );
%
pLMin = mygetfield( prm, "pLMin", 0.0 );
pLMax = mygetfield( prm, "pLMax", 20.0 );
if ( valLev >= VALLEV__LOW )
	assert( isrealorinfscalar(pLMin) );
	assert( isrealorinfscalar(pLMax) );
end
assert( pLMin <= pL_initial );
assert( pLMax >= pL_initial );
%
pRMin = mygetfield( prm, "pRMin", 0.0 );
pRMax = mygetfield( prm, "pRMax", 20.0 );
if ( valLev >= VALLEV__LOW )
	assert( isrealorinfscalar(pRMin) );
	assert( isrealorinfscalar(pRMax) );
end
assert( pRMin <= pR_initial );
assert( pRMax >= pR_initial );
%
%
% Unpack & validate weights.
wVals = mygetfield( prm, "wVals", ones(size(xVals))/numPts );
if ( valLev >= VALLEV__LOW )
	assert( isrealarray(wVals,[1,numPts]) );
	noWValIsNegative = (0==sum(wVals<0.0));
	assert( noWValIsNegative );
	assert( sum(wVals)>0.0 );
end
dVals = sqrt(wVals);
omegaScale = 0.5*sum(wVals.*(fVals.^2));
%
%
% Set quality parameters, halting criteria for "success".
qualLev = mygetfield( prm, "qualLev", QUALLEV__MEDIUM );
if ( qualLev >= QUALLEV__EXTREME )
	tempTol = 0.0;
elseif ( qualLev >= QUALLEV__VERY_HIGH )
	tempTol = 1.0E-15;
elseif ( qualLev >= QUALLEV__HIGH )
	tempTol = 1.0E-10;
elseif ( qualLev >= QUALLEV__MEDIUM )
	tempTol = 1.0E-6;
elseif ( qualLev >= QUALLEV__LOW )
	tempTol = 1.0E-3;
elseif ( qualLev >= QUALLEV__VERY_LOW )
	tempTol = 1.0E-1;
else
	tempTol = 1.0;
end
xExtTol = mygetfield( prm, "xExtTol", tempTol*(xMax-xMin) );
pLTol = mygetfield( prm, "pLTol", tempTol );
pRTol = mygetfield( prm, "pRTol", tempTol );
clear tempTol;
if ( valLev >= VALLEV__LOW )
	assert( isrealscalar(xExtTol) );
	assert( isrealscalar(pLTol) );
	assert( isrealscalar(pRTol) );
	assert( 0.0 <= xExtTol );
	assert( 0.0 <= pLTol );
	assert( 0.0 <= pRTol );
end
%
%
% Set investment parameters, halting criteria for "imposed stop".
investLev = mygetfield( prm, "investLev", INVESTLEV__MEDIUM );
if ( investLev >= INVESTLEV__EXTREME )
	maxNumIter = Inf;
elseif ( investLev >= INVESTLEV__VERY_HIGH )
	maxNumIter = 10000;
elseif ( investLev >= INVESTLEV__HIGH )
	maxNumIter = 500;
elseif ( investLev >= INVESTLEV__MEDIUM )
	maxNumIter = 20;
elseif ( investLevalLev >= INVESTLEV__LOW )
	maxNumIter = 3;
elseif ( investLev >= INVESTLEV__VERY_LOW )
	maxNumIter = 1;
else
	maxNumIter = 0;
end
maxNumIter = mygetfield( prm, "maxNumIter", maxNumIter );
if ( valLev >= VALLEV__LOW )
	assert( isrealscalar(maxNumIter) );
	assert( 0 <= maxNumIter );
end
%
thisfile = [ "RETURNING FROM " thisFile ];
return;
