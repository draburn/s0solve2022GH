ax = [];
caseNum = 33215776;
msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
switch (caseNum)
case -200
	sizeX = 2;
	sizeF = 2;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[]); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case -20
	sizeX = 3;
	sizeF = 3;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[]); % Calls setprngstates.
	testFuncPrm.ary3K *= 0.0;
	vecX0 = randn(sizeX,1);
case -3
	sizeX = 10;
	sizeF = 10;
	tfpPrm.matJPreMod = ones(sizeF,sizeX);
	tfpPrm.matJPreMod(1,1) = 100.0;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[],true,true,true,tfpPrm); % Calls setprngstates.
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
case 10
	sizeX = 3;
	sizeF = 3;
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = [1:sizeX]';
	testFuncPrm.vecFE = zeros(sizeF,1);
	testFuncPrm.matJ = diag([1:sizeX].^4);
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	vecX0 = zeros(sizeX,1);
case 20
	sizeX = 10;
	sizeF = 10;
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = [1:sizeX]';
	testFuncPrm.vecFE = zeros(sizeF,1);
	testFuncPrm.matJ = diag([1:sizeX].^4);
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	vecX0 = zeros(sizeX,1);
case 21
	% A rather more harshly scaled perfectly linear case.
	sizeX = 10;
	sizeF = 10;
	testFuncPrm.sizeX = sizeX;
	testFuncPrm.sizeF = sizeF;
	testFuncPrm.vecXE = [1:sizeX]';
	testFuncPrm.vecFE = zeros(sizeF,1);
	testFuncPrm.matJ = diag([1:sizeX].^4);
	testFuncPrm.ary3K = zeros(sizeX,sizeX,sizeF);
	vecX0 = [
	  2.35540053044716e-07
	  1.20593217315640e-04
	  4.63288261115714e-03
	  6.13143311393352e-02
	  4.41288800438645e-01
	  1.99003962965689e+00
	  5.35919717615100e+00
	  7.94276375465571e+00
	  9.00000006011073e+00
	  9.83329957602381e+00 ];
case 100
	sizeX = 3;
	sizeF = 3;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 200
	sizeX = 2;
	sizeF = 2;
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
case 57863392
	% Case where scaled grad is very bad but grad and Lev are okay?!?!
	sizeX = 2;
	sizeF = 2;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,57863392); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 17381488
	% Lev has dramatic knee, but it's not too bad.
	sizeX = 2;
	sizeF = 2;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,17381488); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 93023744
	% Lev has knee, but it's good.
	sizeX = 2;
	sizeF = 2;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,93023744); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 60420592
	% Lev has bad knee.
	sizeX = 2;
	sizeF = 2;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,60420592); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 6042059200
	% Hack of 60420592.
	sizeX = 2;
	sizeF = 2;
	tfpPrm.matJPreMod = ones(sizeF,sizeX);
	tfpPrm.matJPreMod(1,1) = 100.0;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,60420592,true,true,true,tfpPrm); % Calls setprngstates.
	vecX0 = [ -1.55494508772561; -1.15714300019889 ];
% Looking at Lev H_patch.
case 37176128
	% gradSeg >~ Lev patch >> OCQ
	sizeX = 10;
	sizeF = 10;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,37176128); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 55566400
	% gradSeg >> OCQ ~ lev patch
	sizeX = 10;
	sizeF = 10;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,55566400); % Calls setprngstates.
	vecX0 = randn(sizeX,1);
case 33215776
	% OCQ has bad knee.
	sizeX = 2;
	sizeF = 2;
	testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,33215776); % Calls setprngstates.
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
matH0_patch = matH0_jtj + vecPhi0_jtj * ( 2.0*vecF0'*vecEta0_jtj ) * (vecPhi0_jtj');
%
vecXE = testFuncPrm.vecXE;
[ vecFE, matJE ] = funchF( vecXE );
[ omegaE, vecGradE, matHE_full ] = funchOmega( vecXE );
%
vecXVals_temp = calcGradSeg_cnstH( vecX0, omega0, vecG0, matH0_jtj, [] );
vecXG_jtj = vecXVals_temp(:,end);
%[ vecFG_jtj, matJG_jtj ] = funchF( vecXG_jtj );
%[ omegaEG_jtj, vecGradG_jtj, matHG_jtj ] = funchOmega( vecXG_jtj );
vecFG_jtj = funchF( vecXG_jtj );
omegaG_jtj = funchOmega( vecXG_jtj );
modelInaccuracyG_jtj = norm( vecF0 + matJ0*(vecXG_jtj-vecX0) - vecFG_jtj ) / norm(vecF0);
%
%
%
doExtCheck = true;
if (doExtCheck)
	assert( sumsq(vecGradE) < sqrt(eps)*( sumsq(vecFE) + sum(sumsq(matJE)) + sum(sumsq(matHE_full)) ) );
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
doGradSeg_cnstH_jtj_scaled = false;
doGrad_cnstH_patch = false;
doLev_cnstH_patch = false;
doGradSeg_cnstH_patch = false;
doGradSeg_cnstH_patch_scaled = false;
doGrad_cnstH_full = false;
doLev_cnstH_full = false;
doGradSeg_cnstH_full = false;
doGradSeg_cnstH_full_scaled = false;
doPostGradJTJ_minXi0 = false;
doPostGradJTJ_grad_cnstH = false;
doPostGradJTJ_lev_cnstH = false;
doPostGradJTJ_gradSeg_cnstH = false;
