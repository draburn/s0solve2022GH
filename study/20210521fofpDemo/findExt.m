thisFile = "findExt";

% DRaburn 2021.05.24.2335
% Analyze monotonicity of points.
assert( isrealarray(xVals,[1,numPts]) );
assert( isrealarray(fVals,[1,numPts]) );
assert( 2 <= numPts );
xValsAreStrictlyIncreasing = (0==sum( 0.0 >= diff(xVals) ));
assert(xValsAreStrictlyIncreasing);
numExtrema = 0;
extremaIndices = [];
for n=2:numPts-1
if ( (fVals(n+1)-fVals(n))*(fVals(n)-fVals(n-1)) <= 0.0 )
	numExtrema++;
	extremaIndices = [ extremaIndices, n ];
end
end
assert( 1 >= numExtrema );
if ( 0 == numExtrema )
	hasPtWiseExtrema = false;
else
	hasPtWiseExtrema = true;
	ptWiseExtremaIndex = extremaIndices
	assert( isposintscalar(ptWiseExtremaIndex) );
end

% Calculate our "super points".
% These are, in some sense (f')/(f'').
numSuperPts = numPts-2;
for m=1:numSuperPts
	vecX = xVals(m:m+2)';
	matX = [ ones(3,1), vecX, vecX.^2 ];
	vecF = fVals(m:m+2)';
	vecC = matX \ vecF;
	c0 = vecC(1);
	c1 = vecC(2);
	c2 = vecC(3);
	%
	% What should super_xVal be???
	super_xVals(m) = 0.5*( vecX(1) + vecX(3) );
	super_xVals_probablyBad(m) = vecX(2);
	%
	super_hVals(m) = 0.5*c1/c2 + super_xVals(m);
end
%
% From our "super points", find a best-fit bigX and p.
super_matH = [ ones(numSuperPts,1), super_hVals' ];
super_vecRes = [ super_hVals' + super_xVals' ];
super_vecC = super_matH \ super_vecRes;
%
%
bigX_initial = super_vecC(1);
bigP_initial = super_vecC(2);
msg( thisFile, __LINE__, sprintf( "bigX_initial = %g.", bigX_initial ) );


leftIndex_initial = 0;
while (1)
	if ( leftIndex_initial == numPts )
		break;
	elseif ( xVals(leftIndex_initial+1) >= bigX_initial )
		break;
	else
		leftIndex_initial++;
	end
end
msg( thisFile, __LINE__, sprintf( "leftIndex_initial = %d.", leftIndex_initial ) );
if (hasPtWiseExtrema)
	if ( leftIndex_initial == ptWiseExtremaIndex-1  )
		leftIndex_alternate = ptWiseExtremaIndex;
		bigX_initial2 = 0.99*xVals(ptWiseExtremaIndex) + 0.01*xVals(ptWiseExtremaIndex+1);
	elseif ( leftIndex_initial == ptWiseExtremaIndex )
		leftIndex_alternate = ptWiseExtremaIndex-1;
		bigX_initial2 = 0.99*xVals(ptWiseExtremaIndex) + 0.01*xVals(ptWiseExtremaIndex-1);
	elseif ( leftIndex_initial > ptWiseExtremaIndex )
		msg( thisFile, __LINE__, sprintf( ...
		  "bigX_initial is too far to right. Moving to right side of exterma.") );
		bigX_initial = 0.99*xVals(ptWiseExtremaIndex+1) + 0.01*xVals(ptWiseExtremaIndex);
		leftIndex_initial = ptWiseExtremaIndex;
		leftIndex_alternate = ptWiseExtremaIndex-1;
		bigX_initial2 = 0.99*xVals(ptWiseExtremaIndex) + 0.01*xVals(ptWiseExtremaIndex-1);
	else
		msg( thisFile, __LINE__, sprintf( ...
		  "bigX_initial is too far to left. Moving to left side of exterma.") );
		bigX_initial = 0.99*xVals(ptWiseExtremaIndex-1) + 0.01*xVals(ptWiseExtremaIndex);
		leftIndex_initial = ptWiseExtremaIndex;
		leftIndex_alternate = ptWiseExtremaIndex+1;
		bigX_initial2 = 0.99*xVals(ptWiseExtremaIndex) + 0.01*xVals(ptWiseExtremaIndex+1);
	end
	bigX_altIntBoundary = xVals(ptWiseExtremaIndex);
elseif ( leftIndex_initial == 0 )
	leftIndex_alternate = 1;
	bigX_altIntBoundary = xVals(1);
	bigX_initial2 = 0.99*xVals(1) + 0.01*xVals(2);
elseif ( leftIndex_initial == 1 )
	leftIndex_alternate = 0;
	bigX_altIntBoundary = xVals(1);
	bigX_initial2 = 1.01*xVals(1) - 0.01*xVals(2);
elseif ( leftIndex_initial == numPts-1 )
	leftIndex_alternate = numPts;
	bigX_altIntBoundary = xVals(numPts);
	bigX_initial2 = 1.01*xVals(numPts) - 0.01*xVals(numPts-1);
elseif ( leftIndex_initial == numPts )
	leftIndex_alternate = numPts-1;
	bigX_altIntBoundary = xVals(numPts);
	bigX_initial2 = 0.99*xVals(numPts) + 0.01*xVals(numPts-1);
elseif ( bigX_initial < 0.5*(xVals(2)+xVals(nmPts-1)) )
	msg( thisFile, __LINE__, sprintf( ...
	  "bigX_initial is not close to left edge. Moving towards left edge.") );
	bigX_initial = 0.99*xVals(2) + 0.01*xVals(1);
	leftIndex_initial = 1;
	leftIndex_alternate = 0;
	bigX_altIntBoundary = xVals(1);
	bigX_initial2 = 1.01*xVals(1) - 0.01*xVals(2);
else
	msg( thisFile, __LINE__, sprintf( ...
	  "bigX_initial is not close to right edge. Moving towards right edge.") );
	bigX_initial = 0.99*xVals(numPts-1) + 0.01*xVals(numPts);
	leftIndex_initial = numPts-1;
	leftIndex_alternate = numPts;
	bigX_altIntBoundary = xVals(numPts);
	bigX_initial2 = 1.01*xVals(numPts) - 0.01*xVals(numPts-1);
end
msg( thisFile, __LINE__, sprintf( "leftIndex_alternate = %d.", leftIndex_alternate ) );
haveVisitedAlternateInterval = false; % For now. Not counting BT trials.


matK_initial = [ ones(numPts,1), abs(xVals'-bigX_initial).^bigP_initial ];
vecF = fVals';
vecBigF_initial = matK_initial \ vecF;
bigF0_initial = vecBigF_initial(1);
bigF1_initial = vecBigF_initial(2);
rhoVals_initial = findExt__rho( xVals, fVals, bigX_initial, bigP_initial );
res_initial = 0.5*sum(rhoVals_initial.^2);

%
%
% Note that res is 0.5*sum( rhoVals.^2 ).
resTarget = numPts * eps^1.5;
resAccept = 1E2*resTarget;
maxLoopCount = 1000;
minResDecrease = 1E-8;

loopCount = 0;
btCount = 0;
vecDelta0 = zeros(2,1);
vecDelta = zeros(2,1);
bigX = bigX_initial;
bigP = bigP_initial;
res = res_initial;
findExt__loop; thisFile = "findExt";

if ( res > resAccept )
if (haveVisitedAlternateInterval)
	msg( thisFile, __LINE__, "Have already visited alternate interval." );
else
	haveVisitedAlternateInterval = true;
	% Re-do from "alternate" interval.
	% We assume our initial guess was at least close.
	% Note that the loop might not prevent jumping between intervals.
	msg( thisFile, __LINE__, "Examining alternate interval." );
	bigX_firstSolve = bigX;
	bigP_firstSolve = bigP;
	res_firstSolve = res;
	%
	bigP_initial2 = bigP_initial;
	matK_initial2 = [ ones(numPts,1), abs(xVals'-bigX_initial2).^bigP_initial2 ];
	vecF = fVals';
	vecBigF_initial2 = matK_initial2 \ vecF;
	bigF0_initial2 = vecBigF_initial2(1);
	bigF1_initial2 = vecBigF_initial2(2);
	rhoVals_initial2 = findExt__rho( xVals, fVals, bigX_initial2, bigP_initial2 );
	res_initial2 = 0.5*sum(rhoVals_initial2.^2);
	%
	bigX = bigX_initial2;
	bigP = bigP_initial2;
	res = res_initial2;
	findExt__loop; thisFile = "findExt";
	%
	bigX_secondSolve = bigX;
	bigP_secondSolve = bigP;
	res_secondSolve = res;
	if ( res_secondSolve >= res_firstSolve )
		msg( thisFile, __LINE__, "Alternate interval was WORSE." );
		bigX = bigX_firstSolve;
		bigP = bigP_firstSolve;
		res = res_firstSolve;
	end
end
end

% Finally: Find corresponding F0 and F1
matK = [ ones(numPts,1), abs(xVals'-bigX).^bigP ];
vecF = fVals';
vecBigF = matK \ vecF;
bigF0 = vecBigF(1);
bigF1 = vecBigF(2);
return;
