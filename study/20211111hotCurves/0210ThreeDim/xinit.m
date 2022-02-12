ax = [];
caseNum = -2;
msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
switch (caseNum)
case -20
	sizeX = 3;
	sizeF = 3;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[]); % Calls setprngstates.
	testFuncPrm.ary3K *= 0.0;
	vecX0 = randn(sizeX,1);
case -2
	sizeX = 10;
	sizeF = 10;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[]); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case -1
	sizeX = 3;
	sizeF = 3;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[]); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 0
	sizeX = 3;
	sizeF = 3;
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = ones(sizeX,1);
	testFuncPrm.vecFE = zeros(sizeF,1);
	testFuncPrm.matJ = eye(sizeF,sizeX);
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	vecX0 = zeros(sizeX,1);
case 1
	sizeX = 3;
	sizeF = 3;
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = [1:sizeX]';
	testFuncPrm.vecFE = zeros(sizeF,1);
	testFuncPrm.matJ = diag([1:sizeX]);
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	vecX0 = zeros(sizeX,1);
case 100
	sizeX = 3;
	sizeF = 3;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 18153952
	% Full steps are good; grad segment doesn't go far enough.
	sizeX = 3;
	sizeF = 3;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,18153952); % Calls setprngstates.
	testFuncPrm.ary3K *= 0.0;
	vecX0 = randn(sizeX,1);
case 82177120
	sizeX = 3;
	sizeF = 3;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,82177120); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
otherwise
	error( "Invalid value of switch." );
end
%
%
%
funchF = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
funchOmega = @(dummyX)( testfunc2021_funcOmega(dummyX,testFuncPrm) );
%
[ vecF0, matJ0 ] = funchF( vecX0 );
[ omega0, vecG0, matH0_full ] = funchOmega( vecX0 );
%
matH0_jtj = matJ0'*matJ0;
[ vecPhi0_full, vecEta0_full ] = calcOCQTerms( vecX0, funchF, matH0_full );
[ vecPhi0_jtj,  vecEta0_jtj  ] = calcOCQTerms( vecX0, funchF, matH0_jtj  );
%
vecXE = testFuncPrm.vecXE;
[ vecFE, matJE ] = funchF( vecXE );
[ omegaE, vecGrad, matHE_full ] = funchOmega( vecXE );
%
%
%
doExtCheck = true;
if (doExtCheck)
	assert( sumsq(vecGrad) < sqrt(eps)*( sumsq(vecFE) + sum(sumsq(matJE)) + sum(sumsq(matHE_full)) ) );
end
%
%
doKSymCheck = true;
if (doKSymCheck)
for k=1:testFuncPrm.sizeF
for m=1:testFuncPrm.sizeX
for n=1:testFuncPrm.sizeX
	assert( fleq( testFuncPrm.ary3K(m,n,k), testFuncPrm.ary3K(n,m,k) ) );
end
end
end
end
%
%
%
% Set these to false for now.
doGrad_alytG = false;
doLev_dispena = false;
doLev_minford = false;
doFOCQ_minXi0_jtj = false;
doFOCQ_minXi0_fullish = false;
doFOCQ_L0jtj = false;
doFOCQ_C0jtj = false;
doFOCQ_R0jtj = false;
doGrad_cnstH_jtj = false;
doLev_cnstH_jtj = false;
doGradSeg_cnstH_jtj = false;
doGrad_cnstH_full = false;
doLev_cnstH_full = false;
doGradSeg_cnstH_full = false;
