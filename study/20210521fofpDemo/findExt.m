thisFile = "findExt";

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
bigX_initial = super_vecC(1);
bigP_initial = super_vecC(2);
matK_initial = [ ones(numPts,1), abs(xVals'-bigX_initial).^bigP_initial ];
vecF = fVals';
vecBigF_initial = matK_initial \ vecF;
bigF0_initial = vecBigF_initial(1);
bigF1_initial = vecBigF_initial(2);
rhoVals_initial = findExt__rho( xVals, fVals, bigX_initial, bigP_initial );
res_initial = 0.5*sum(rhoVals_initial.^2);
%
%
resTarget = numPts * eps^1.5;
resAccept = 1E2*resTarget;
maxLoopCount = 1000;

loopCount = 0;
btCount = 0;
vecDelta0 = zeros(2,1);
vecDelta = zeros(2,1);
bigX = bigX_initial;
bigP = bigP_initial;
res = res_initial;
findExt__loop; thisFile = "findExt";

if ( res > resAccept )
	% Re-do from "alternate" interval.
	% We assume our initial guess was at least close.
	% Note that the loop might not prevent jumping between intervals.
	msg( thisFile, __LINE__, "Examining alternate interval." );
	bigX_firstSolve = bigX;
	bigP_firstSolve = bigP;
	res_firstSolve = res;
	%
	msg( thisFile, __LINE__, "HACK!!!" );
	msg( thisFile, __LINE__, "HACK!!!" );
	msg( thisFile, __LINE__, "PLACEHOLDER HACK: FORCING bigX_initial2 = 1.0." );
	bigX_initial2 = 1.0;
	msg( thisFile, __LINE__, "HACK!!!" );
	msg( thisFile, __LINE__, "HACK!!!" );
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
	if ( res_secondSolve < res )
		msg( thisFile, __LINE__, "Alternate interval was WORSE." );
		bigX = bigX_firstSovle;
		bigP = bigP_firstSovle;
		res = res_firstSolve;
	end
end

% Finally: Find corresponding F0 and F1
matK = [ ones(numPts,1), abs(xVals'-bigX).^bigP ];
vecF = fVals';
vecBigF = matK \ vecF;
bigF0 = vecBigF(1);
bigF1 = vecBigF(2);
return;
