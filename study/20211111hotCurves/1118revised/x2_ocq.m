clear;
thisFile = "x2_ocq";
commondefs;
numFigs = 0;
tic();
%
useAxisEqual = false;
sizeX = 2;
sizeF = 2;
if (0)
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,60525840,true,false,false); % x8 Green bad!
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,45342256,true,false,false); % x8 Curves 2 are short.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,60663152,true,false,false); % x8 red misses.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,91475584,true,false,false); % Dramatic benefit.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,72643728,true,false,false); %Hmm...
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,12186496,true,false,false);
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,45553360,true,false,false); % Lev and Grad cnvg to different!
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,85541360,true,false,false); % Lev bonks.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,81404512,true,false,false); % Lev bonks HARD.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,82178352,true,false,false); % First Lev step is big.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,96042480,true,false,false); % Rattle; First Lev step is funky.
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,53162192,true,false,false); % Lev bonks.
	% 1118, good test cases...
	testFuncPrm = testFunc_genPrm(sizeX,sizeF,71832832,true,false,false);
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,61409680,true,false,false);
	%testFuncPrm = testFunc_genPrm(sizeX,sizeF,80724832,true,false,false); % "Cresent".
	vecX0 = [ 3.0; 3.0 ];
	ax = [ -3.0, 6.0, -3.0, 6.0 ];
elseif (0)
	% One-Component Quadratic.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.1; 0.0 ];
	testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
	%%%testFuncPrm.ary3K(:,:,1) = [ 0.01, 0.0; 0.0, 0.0 ];
	%%%testFuncPrm.ary3K(:,:,2) = 10*[ 0.05, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.1, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 1.0, 0.0; 0.0, 0.0 ];
	%vecX0 = [ 5.0; 0.0 ];
	vecX0 = [ 2.2; -1.3 ];
	ax = [ -1.0, 3.0, -2, 0.5 ];
elseif (0)
	% Narrow valley.
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 0.0; 0.0 ];
	testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
	%%%testFuncPrm.ary3K(:,:,1) = [ 0.01, 0.0; 0.0, 0.0 ];
	%%%testFuncPrm.ary3K(:,:,2) = 10*[ 0.05, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,1) = [ 0.1, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.5, 0.0; 0.0, 0.0 ];
	%vecX0 = [ 5.0; 0.0 ];
	vecX0 = [ 2.2; -1.3 ];
	ax = [ -1.0, 5.5, -2, 0.5 ];
elseif (0)
	% Dbl sepx.
	testFuncPrm = testFunc_genPrm(sizeX,sizeF,25094192,true,false,false);
	vecX0 = [ 4.0; -0.11742 ]; % Some go to UR, one goes to LL.
	%vecX0 = [ 4.0017; -0.11742 ]; % All dead end at LX.
	%vecX0 = [ 4.002; -0.11742 ]; % Ones goes to LR, two dead end at LX.
	%vecX0 = [ -1.413; -0.405 ]; % Start near outer sepx.
	%vecX0 = [ -0.33; -0.5 ]; % Get split.
	%vecX0 = [ -2.2; -0.22 ]; % Dead ends at UX.
	%
	ax = [ -5.0, 5.0, -5.0, 5.0 ];
	%ax = [ -1.414 -1.412 -0.406 -0.404 ]; % Zoom on double sepratrix.
	%ax = [ -2.5 0 -0.6 -0.2 ];
	%ax = [ -1.0, -0.7, -1.5, -1.3 ];
elseif (1)
	testFuncPrm.sizeX = 2;
	testFuncPrm.sizeF = 2;
	testFuncPrm.vecXE = [ 0.0; 0.0 ];
	testFuncPrm.vecFE = [ 1.0; 0.0 ];
	testFuncPrm.matJ = [ 0.0, 0.0; 0.0, 1.0 ];
	testFuncPrm.ary3K(:,:,1) = [ -0.1, 0.0; 0.0, 0.0 ];
	testFuncPrm.ary3K(:,:,2) = [ 0.1, 0.0; 0.0, 0.0 ];
	vecX0 = [ 0.01; 3.0 ];
	ax = [ -3.0, 6.0, -3.0, 6.0 ];
end

doOCQLevDev = true;
if (doOCQLevDev)
	%vecX0 = [ 2; -1.35 ]
	%vecX0 = [ 0.9; -0.4 ]
	%vecX0 = [ 0.2; -0.01 ]
	%vecX0 = [ 0.5; -0.5 ]
	%vecX0 = [ 0.71; -0.26 ]
	%vecX0 = [ -0.014; 0.001 ]
	%ax = [ -0.5, 0.5, -0.05, 0.05 ];
	%ax = 0.1*[-1,1,-1,1];
	%vecX0 = [ 1.81; -1.446 ]
	%vecX0 = [ 1.0; -0.4 ]
	%vecX0 = [ 0.4; -0.025 ]
	%vecX0 = [ 0.1; 0.01 ]
	%ax = 0.5*[1,-1,0.1,-0.1]
	vecX0 = [ 1.5; 0.0 ]
	%
	msg( thisFile, __LINE__, "Note: OCQ Lev params may be generated from exact values and not finite-difference." );
	[ omega0, vecF0, matJ0 ] = testFunc_evalDeriv( vecX0, testFuncPrm );
	switch (1)
	case {1}
		% Use eig min.
		msg( thisFile, __LINE__, "Using vecPhi per eigMin." );
		matH0 = matJ0'*matJ0
		[ matPsi, matLam ] = eig( matH0 );
		[ lambdaMin, nOfMin ] = min(diag(matLam))
		assert( lambdaMin >= 0.0 );
		vecPhi = matPsi(:,nOfMin)
	case {2}
		% Use grad.
		msg( thisFile, __LINE__, "Using vecPhi per grad." );
		vecG0 = matJ0'*vecF0;
		vecPhi = vecG0/norm(vecG0)
	case {3}
		% Use deltaN.
		msg( thisFile, __LINE__, "Using vecPhi per deltaN." );
		matH0 = matJ0'*matJ0;
		vecG0 = matJ0'*vecF0;
		vecDeltaN = -matH0\vecG0;
		vecPhi = vecDeltaN/norm(vecDeltaN)
	case {4}
		% Use forced.
		msg( thisFile, __LINE__, "Using vecPhi HACK." );
		vecPhi = [ 1.0; 0.0 ]
		%vecPhi = [ 1.0; 1.0 ] % So, so wrong.
		vecPhi /= norm(vecPhi);
	otherwise
		msg( thisFile, __LINE__, "Invalid value of switch for vecPhi." );
		error( "Invalid value of switch for vecPhi." );
	end
	myEps = eps025;
	vecFP = testFunc_eval( vecX0 + myEps*vecPhi, testFuncPrm );
	vecFM = testFunc_eval( vecX0 - myEps*vecPhi, testFuncPrm );
	vecEta = ( vecFP + vecFM - 2.0*vecF0 ) / ( 2.0 * (myEps^2) )
	% Under 2021.11.24 def, vecEta is half of the relevant ary3K.
	%
	callAsScript = false;
	if (callAsScript)
		matJ = matJ0;
		calcOCQLevCurve
		matY_ocqLev = matX;
	else
		[ matY_ocqLev, ocqLevDat ] = calcOCQLevCurve( vecF0, matJ0, vecEta, vecPhi, vecX0 );
	end
	%return;
end

%
funchF = @(x)( testFunc_eval(x,testFuncPrm) );
funchOmegaG = @(x)( testFunc_evalOmegaG(x,testFuncPrm) );
%
tempPrm = [];
tempPrm.useScaling = false;
matY_grad_old = OLD_calcHOTGradCurveRK4( testFuncPrm, vecX0, tempPrm );
%matY_grad_simple = calcHOTGradCurve_simple( funchOmegaG, vecX0 );
%
tempPrm = [];
%tempPrm.maxStepSize = 1000.0;
%%%matY_grad_simple = OLD_calcHOTGradCurve2( testFuncPrm, vecX0, tempPrm );
matY_grad_simple = OLD_calcHOTGradCurveNewt( testFuncPrm, vecX0, tempPrm );
%
tempPrm = [];
matY_grad1 = OLD_calcHOTGradCurve2Heun( testFuncPrm, vecX0 );
%matY_grad1 = calcHOTGradCurve2( funchOmegaG, vecX0 );
%
matD_grad_old = matY_grad_old - repmat( vecX0, 1, size(matY_grad_old,2) );
matD_grad_simple = matY_grad_simple - repmat( vecX0, 1, size(matY_grad_simple,2) );
matD_grad1 = matY_grad1 - repmat( vecX0, 1, size(matY_grad1,2) );
matD_ocqLev = matY_ocqLev - repmat( vecX0, 1, size(matY_ocqLev,2) );
%
matF_grad_old = funchF(matY_grad_old,testFuncPrm);
matF_grad_simple = funchF(matY_grad_simple,testFuncPrm);
matF_grad1 = funchF(matY_grad1,testFuncPrm);
matF_ocqLev = funchF(matY_ocqLev,testFuncPrm);
%%%%%%%%%%%%%%%%%%%%%%
%
%
p1;
v1;
%
toc();
