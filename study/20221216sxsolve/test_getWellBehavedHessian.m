printf( "\n\n\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n" );
msg( __FILE__, __LINE__, "kthxhai" )
clear;
setprngstates(0);
%
sz = 5;
vecXCrit = randn(sz,1);
vecXIG = zeros(sz,1); %T = "temporary"?
%
for cid = [ 0, 10, 100, 200, 300 ]
echo__cid = cid
switch (cid)
case 0
	fCrit = 0.0;
	matH = mtm(randn(sz,sz));
case 10
	fCrit = 0.0;
	matH = mtm(randn(sz-1,sz)) + sqrt(eps)*mtm(randn(1,sz));
case 100
	fCrit = 0.0;
	matH = mtm(randn(1,sz));
case 200
	fCrit = -0.5;
	matH = mtm(randn(sz,sz));
case 300
	fCrit = 0.5;
	matH = -mtm(randn(sz,sz));
otherwise
	fCrit = 0.0;
	matH = -mtm(randn(sz,sz));
endswitch
%
vecDIG = vecXIG - vecXCrit;
vecGIG = matH*vecDIG;
fIG = fCrit + vecDIG'*vecGIG/2.0;
if ( fIG < 0.0 )
	fCrit += 2.0*abs(fIG);
	fIG = fCrit + vecDIG'*vecGIG/2.0;
endif
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
muTot = gwbhDat.muTot
fNewt = fCrit + vecDNewt'*vecGNewt/2.0
endfor %cid
msg( __FILE__, __LINE__, "kthxbai" )
