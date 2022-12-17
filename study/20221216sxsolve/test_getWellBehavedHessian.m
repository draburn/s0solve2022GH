printf( "\n\n\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n" );
msg( __FILE__, __LINE__, "kthxhai" )
clear;
setprngstates(0);
%
sz = 5;
vecXCrit = randn(sz,1);
vecXIG = zeros(sz,1); %T = "temporary"?
%
for cid = [ 110, 120, 130, 140, 150, 210, 220, 230, 240, 250, 310, 320, 330, 340, 350, 410, 420, 430, 440, 450, 510, 520, 530, 540, 550, 560, 570, 580, 590 ]
echo__cid = cid
switch (cid)
% 100~
% Actually pos-def cases with fMin >= 0.0.
case 110
	fCrit = 0.0;
	matH = mtm(randn(sz,sz));
case 120
	fCrit = 0.0;
	matH = mtm(randn(sz-1,sz)) + sqrt(eps)*mtm(randn(sz,sz));
case 130
	fCrit = 0.0;
	matH = mtm(randn(sz-1,sz)) + sqrt(eps)*mtm(randn(1,sz));
case 140
	fCrit = 0.0;
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(sz,sz));
case 150
	fCrit = 0.0;
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(sz-1,sz));
% 200~
% Actually pos-semi-def cases with fMin >= 0.0.
case 210
	% ... Assuming sz >= 3.
	fCrit = 0.0;
	matH = mtm(randn(sz-2,sz)) + sqrt(eps)*mtm(randn(1,sz));
case 220
	fCrit = 0.0;
	matH = mtm(randn(sz-1,sz));
case 230
	% ... Assuming sz >= 3.
	fCrit = 0.0;
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(sz-2,sz));
case 240
	fCrit = 0.0;
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(1,sz));
case 250
	fCrit = 0.0;
	matH = mtm(randn(1,sz));
% 300~
% Atually pos-def cases with fMin < 0.0.
% Copy from "actual pos-def cases with fMin >= 0.0" except "fCrit = []" instead.
case 310
	fCrit = [];
	matH = mtm(randn(sz,sz));
case 320
	fCrit = [];
	matH = mtm(randn(sz-1,sz)) + sqrt(eps)*mtm(randn(sz,sz));
case 330
	fCrit = [];
	matH = mtm(randn(sz-1,sz)) + sqrt(eps)*mtm(randn(1,sz));
case 340
	fCrit = [];
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(sz,sz));
case 350
	fCrit = [];
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(sz-1,sz));
% 400~
% Actual pos-semi-def cases with fMin < 0.0.
% Copy from "actual pos-semi-def cases with fMin >= 0.0" except "fCrit = []" instead.
case 410
	% ... Assuming sz >= 3.
	fCrit = [];
	matH = mtm(randn(sz-2,sz)) + sqrt(eps)*mtm(randn(1,sz));
case 420
	fCrit = [];
	matH = mtm(randn(sz-1,sz));
case 430
	% ... Assuming sz >= 3.
	fCrit = [];
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(sz-2,sz));
case 440
	fCrit = [];
	matH = mtm(randn(1,sz)) + sqrt(eps)*mtm(randn(1,sz));
case 450
	fCrit = [];
	matH = mtm(randn(1,sz));
% 500~
% Actual has-neg cases with fInitialGuess > 0.0.
case 510
	fCrit = [];
	matH = mtm(randn(sz-1,sz)) - sqrt(eps)*mtm(randn(1,sz));
case 520
	fCrit = [];
	matH = mtm(randn(1,sz)) - sqrt(eps)*mtm(randn(sz-1,sz));
case 530
	fCrit = [];
	matH = mtm(randn(sz-1,sz)) - mtm(randn(1,sz));
case 540
	fCrit = [];
	matH = mtm(randn(1,sz)) - mtm(randn(sz-1,sz));
case 550
	fCrit = [];
	matH = sqrt(eps)*mtm(randn(sz-1,sz)) - mtm(randn(1,sz));
case 560
	fCrit = [];
	matH = sqrt(eps)*mtm(randn(1,sz)) - mtm(randn(sz-1,sz));
case 570
	fCrit = [];
	matH = -mtm(randn(1,sz));
case 580
	fCrit = [];
	matH = -mtm(randn(sz-1,sz));
case 590
	fCrit = [];
	matH = -mtm(randn(sz,sz));
otherwise
	error( "Invalid case." );
endswitch
%
vecDIG = vecXIG - vecXCrit;
vecGIG = matH*vecDIG;
if ( isempty(fCrit) )
	%msg( __FILE__, __LINE__, "THIS (BRANCH) IS NOT (FULLY) FUNCTIONAL; expect exceptions." );
	% Make it so that fIG > 0 but fCrit < 0.
	foo = vecDIG'*vecGIG/2.0;
	fCrit = 1.0 - foo;
endif
fIG = fCrit + vecDIG'*vecGIG/2.0;
%
prm = [];
prm.f0 = fIG;
prm.vecG0 = vecGIG;
fMinAllowed = 0.0;
prm.fMinAllowed = fMinAllowed;
%msg( __FILE__, __LINE__, "Calling getWellBehavedHessian()..." );
[ matHWB, gwbhDat ] = getWellBehavedHessian( matH, prm );
%msg( __FILE__, __LINE__, "Back from getWellBehavedHessian()." );
vecXNewt = vecXIG - matHWB\vecGIG;
vecDNewt = vecXNewt - vecXCrit;
vecGNewt = matH*vecDNewt;
afterVals = [ gwbhDat.muTot, min(eig(matHWB)), fIG, fCrit + vecDNewt'*vecGNewt/2.0 ]
endfor %cid
msg( __FILE__, __LINE__, "kthxbai" )
