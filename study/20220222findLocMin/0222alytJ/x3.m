	clear;
	numFigs = 0;
	%
	%caseNum = 104; % kupd with inter is best; Newt with TR and kupd sans inter are okay.
	%caseNum = 200; % Newt with TR is best (before), kupd with inter is best-ish (after TR tweak);
	%   kupd with inter and kupd sans inter is okay.
	%caseNum = 300; % blind Newt regu is best; everything except blind Newt is okay.
	%caseNum = 38104560; % kupd winter appreciably better than newt with tr.
	caseNum = 999;
	%caseNum = 9392336; % CAUSED AN ERROR.
	%caseNum = 99041968; % Causes chol(matH) okay but chol(matH+mu*matI) (mu~=eps) fails!
	%caseNum = 41765088; % Sans TR is better?!?! (Different roots, ja?)
	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
	switch (caseNum)
	case 0
		sizeX = 2;
		sizeF = 2;
		vecX0 = zeros(sizeX,1);
		funchFJ = @(dummyX)( funcFJ_trivial(dummyX) );
	case 1
		sizeX = 10;
		sizeF = 10;
		vecX0 = zeros(sizeX,1);
		funchFJ = @(dummyX)( funcFJ_easy(dummyX) );
	case 10
		sizeX = 2;
		sizeF = 2;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0); % Calls setprngstates.
		%echo__vecFE = testFuncPrm.vecFE
		%echo__omegaE = sumsq(testFuncPrm.vecFE)/2.0
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 11
		sizeX = 2;
		sizeF = 2;
		tfpPrm.matJPreMod = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0,true,true,true,tfpPrm); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 20
		sizeX = 2;
		sizeF = 2;
		tfpPrm.matJPreMod = ones(sizeF,sizeX);
		tfpPrm.matJPreMod(1,1) = 100.0;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0,true,true,true,tfpPrm); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 100
		msg( __FILE__, __LINE__, "*** WARNING: This is a 'perfectly balanced' case! ***" );
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 0.0; 0.0 ];
		testFuncPrm.vecFE = [ 0.0; 1.0 ];
		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = ones(sizeX,1);
	case 101
		msg( __FILE__, __LINE__, "*** WARNING: This is a 'perfectly balanced' case! ***" );
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 0.0; 0.0 ];
		testFuncPrm.vecFE = [ 0.0; 1.0 ];
		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = 0.01*ones(sizeX,1);
	case 102
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 0.0; 0.0 ];
		testFuncPrm.vecFE = [ 0.0; 1.0 ];
		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = [ 1.0; 0.9 ];
	case 103
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 0.0; 0.0 ];
		testFuncPrm.vecFE = [ 0.0; 1.0 ];
		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = 0.01*[ 1.0; 0.9 ];
	case 104
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 0.0; 0.0 ];
		testFuncPrm.vecFE = [ 0.0; 1.0 ];
		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 2.0, 0.0; 0.0, 2.0 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = [ 1.0; 2.0 ];
	case 200
		sizeX = 20;
		sizeF = 20;
		tfpPrm.matJPreMod = ones(sizeF,sizeX);
		tfpPrm.matJPreMod(1,1) = 10.0;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,1,true,true,true,tfpPrm); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 300
		sizeX = 2;
		sizeF = 2;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,0,false,true,true); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 990
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[],false,true,true); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 41765088
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,41765088,false,true,true); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 999
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[]); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 9392336
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,9392336); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 38104560
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,38104560); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 99041968
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,99041968); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	otherwise
		error( "Invalid caseNum." );
	endswitch
	%
	vecF0 = funchFJ(vecX0);
	omegaE = sumsq(testFuncPrm.vecFE)/2.0 - eps*sumsq(vecF0)/2.0;
	%
	%
	%
	prm_blindNewt.stepType = 10;
	[ vecXF_blindNewt, datOut_blindNewt ] = findLocMin_alytJ( vecX0, funchFJ, prm_blindNewt );
	%
	prm_blindNewtRegu.stepType = 11;
	[ vecXF_blindNewtRegu, datOut_blindNewtRegu ] = findLocMin_alytJ( vecX0, funchFJ, prm_blindNewtRegu );
	%
	prm_blindNewtPatch.stepType = 20;
	[ vecXF_blindNewtPatch, datOut_blindNewtPatch ] = findLocMin_alytJ( vecX0, funchFJ, prm_blindNewtPatch );
	%
	prm_blindNewtPatchRegu.stepType = 21;
	[ vecXF_blindNewtPatchRegu, datOut_blindNewtPatchRegu ] = findLocMin_alytJ( vecX0, funchFJ, prm_blindNewtPatchRegu );
	%
	prm_newtSansBT.stepType = 30;
	[ vecXF_newtSansTR, datOut_newtSansTR ] = findLocMin_alytJ( vecX0, funchFJ, prm_newtSansBT );
	%
	prm_newtWithTR.stepType = 31;
	[ vecXF_newtWithTR, datOut_newtWithTR ] = findLocMin_alytJ( vecX0, funchFJ, prm_newtWithTR );
	%
	prm_kupdSansInter.stepType = 100;
	[ vecXF_kupdSansInter, datOut_kupdSansInter ] = findLocMin_alytJ( vecX0, funchFJ, prm_kupdSansInter );
	%
	prm_kupdWithInter.stepType = 110;
	[ vecXF_kupdWithInter, datOut_kupdWithInter ] = findLocMin_alytJ( vecX0, funchFJ, prm_kupdWithInter );
	%
	numFigs++; figure(numFigs);
	loglog( ...
	  datOut_blindNewt.fevalCountVals, datOut_blindNewt.omegaVals-omegaE, 'o-', 'markersize', 15, ...
	  datOut_blindNewtRegu.fevalCountVals, datOut_blindNewtRegu.omegaVals-omegaE, 'x-', 'markersize', 15, ...
	  datOut_blindNewtPatch.fevalCountVals, datOut_blindNewtPatch.omegaVals-omegaE, 'p-', 'markersize', 15, ...
	  datOut_blindNewtPatchRegu.fevalCountVals, datOut_blindNewtPatchRegu.omegaVals-omegaE, '+-', 'markersize', 15, ...
	  datOut_newtSansTR.fevalCountVals, datOut_newtSansTR.omegaVals-omegaE, '^-', 'markersize', 15, ...
	  datOut_newtWithTR.fevalCountVals, datOut_newtWithTR.omegaVals-omegaE, 'v-', 'markersize', 15, ...
	  datOut_kupdSansInter.fevalCountVals, datOut_kupdSansInter.omegaVals-omegaE, 's-', 'markersize', 15, ...
	  datOut_kupdWithInter.fevalCountVals, datOut_kupdWithInter.omegaVals-omegaE, '*-', 'markersize', 15 );
	grid on;
	xlabel( "feval count" );
	ylabel( "omega" );
	title( "omega vs feval" );
	legend( ...
	  "blind Newt", ...
	  "blind Newt regu", ...
	  "blind Newt patch", ...
	  "blind Newt patch reg", ...
	  "Newt sans TR", ...
	  "Newt with TR", ...
	  "kupd sans inter", ...
	  "kupd with inter", ...
	  "location", "northeast" );
	%
	numFigs++; figure(numFigs);
	loglog( ...
	  datOut_blindNewt.deltaNormVals, 'o-', 'markersize', 15, ...
	  datOut_blindNewtRegu.deltaNormVals, 'x-', 'markersize', 15, ...
	  datOut_blindNewtPatch.deltaNormVals, 'p-', 'markersize', 15, ...
	  datOut_blindNewtPatchRegu.deltaNormVals, '+-', 'markersize', 15, ...
	  datOut_newtSansTR.deltaNormVals, '^-', 'markersize', 15, ...
	  datOut_newtWithTR.deltaNormVals, 'v-', 'markersize', 15, ...
	  datOut_kupdSansInter.deltaNormVals, 's-', 'markersize', 15, ...
	  datOut_kupdWithInter.deltaNormVals, '*-', 'markersize', 15 );
	grid on;
	xlabel( "iter count" );
	ylabel( "deltaNorm" );
	title( "deltaNorm vs iter" );
	legend( ...
	  "blind Newt", ...
	  "blind Newt regu", ...
	  "blind Newt patch", ...
	  "blind Newt patch reg", ...
	  "Newt sans TR", ...
	  "Newt with TR", ...
	  "kupd sans inter", ...
	  "kupd with inter", ...
	  "location", "northeast" );
