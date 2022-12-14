	clear;
	commondefs;
	findLocMin_gnostic_jupdate2__defs;
	setprngstates(0);
	numFigs = 0;
	%
	%caseNum = 104; % Happens to allow matH = matJ'*matJ + mu*matI!
	%caseNum = 110;
	caseNum = 115;
	%caseNum = 200; % Newt with TR is best (before), kupd with inter is best-ish (after TR tweak);
	%   kupd with inter and kupd sans inter is okay.
	%caseNum = 300; % blind Newt regu is best; everything except blind Newt is okay.
	%caseNum = 400;
	%caseNum = 38104560; % kupd winter appreciably better than newt with tr.
	%caseNum = 999;
	%caseNum = 9392336; % CAUSED AN ERROR.
	%caseNum = 99041968; % Causes chol(matH) okay but chol(matH+mu*matI) (mu~=eps) fails!
	%caseNum = 41765088; % Sans TR is better?!?! (Different roots, ja?)
	%caseNum = 6295760;
	%caseNum = 72631264; % Something is very wrong here! Multiple extermum?
	%caseNum = 13175504; % Winter K is rather helpful.
	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
	switch (caseNum)
	case 0
		sizeX = 2;
		sizeF = 2;
		vecX0 = zeros(sizeX,1);
		%funchFJ = @(dummyX)( funcFJ_trivial(dummyX) );
		funchFJ = @(dummyX)( dummyX-ones(size(dummyX)) );
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
	case 110
		% Like case 104, but, make it so that J'*J + mu*I does NOT happen to be the actual Hessian!
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 0.0; 0.0 ];
		testFuncPrm.vecFE = [ 0.0; 1.0 ];
		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 1.0, 0.0; 0.0, 3.0 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = [ 1.0; 2.0 ];
	case 115
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 0.0; 0.0 ];
		testFuncPrm.vecFE = [ 0.0; 1.0 ];
		testFuncPrm.matJ = [ 1.0, 1.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 10.0, 0.0; 0.0, 30.0 ];
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
	case 400
		% Near linear case; kupd bad!
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 1.0; 1.0 ];
		testFuncPrm.vecFE = [ 0.0; 0.0 ];
		testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 1.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.0, 0.0; 0.0, 0.0 ];
		testFuncPrm.ary3K(:,:,2) = [ 0.1, 0.0; 0.0, 0.0 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = [ 0.0; 0.0 ];
	case 401
		% Near linear case; kupd bad!
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 1.0; 1.0 ];
		testFuncPrm.vecFE = [ 0.0; 0.0 ];
		testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 1.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.1, 0.0; 0.0, 0.2 ];
		testFuncPrm.ary3K(:,:,2) = [ 0.3, 0.0; 0.0, 0.3 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = [ 0.0; 0.0 ];
	case 405
		% Near linear case; kupd bad!
		sizeX = 2;
		sizeF = 2;
		testFuncPrm.sizeX = 2;
		testFuncPrm.sizeF = 2;
		testFuncPrm.vecXE = [ 1.0; 1.0 ];
		testFuncPrm.vecFE = [ 0.0; 0.0 ];
		testFuncPrm.matJ = [ 1.0, 0.0; 0.0, 1.0 ];
		testFuncPrm.ary3K(:,:,1) = [ 0.4, 0.2; 0.2, 0.5 ];
		testFuncPrm.ary3K(:,:,2) = [ 0.3, 0.1; 0.1, 0.3 ];
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = [ 0.0; 0.0 ];
	case 990
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,[],false,true,true); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 72631264
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,72631264,false,true,true); % Calls setprngstates.
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
	case 13175504 % Winter K is rather helpful.
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,13175504); % Calls setprngstates.
		funchFJ = @(dummyX)( testfunc2021_funcF(dummyX,testFuncPrm) );
		vecX0 = zeros(sizeX,1);
	case 6295760
		% Bline Newt regu is good early, but terminates.
		sizeX = 20;
		sizeF = 20;
		testFuncPrm = testfunc2021_genPrm(sizeX,sizeF,6295760); % Calls setprngstates.
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
	omega0 = sumsq(vecF0)/2.0;
	%omegaE = sumsq(testFuncPrm.vecFE)/2.0 - eps*sumsq(vecF0)/2.0;
	%
	doFSolveGnostic = true;
	if (doFSolveGnostic)
		for n = 1 : 16
			tolFun = omega0*10.0^(1-n);
			options_fsolve = optimset( 'Updating', 'on', 'Jacobian', 'off', 'GradObj', 'off', 'TolFun', tolFun );
			[  vecXF_fsolve, vecFF_fsolve, info_fsolve, output_fsolve, matJF_fsolve ] = fsolve( funchFJ, vecX0, options_fsolve );
			%output_fsolve
			omegaVals_fsolve(n) = sumsq(vecFF_fsolve,1)/2.0;
		fevalCountVals_fsolve(n) = output_fsolve.funcCount;
		endfor
		semilogy( fevalCountVals_fsolve, omegaVals_fsolve, 'o-' );
		grid on;
		[  vecXF_fsolve, vecFF_fsolve, info_fsolve, output_fsolve, matJF_fsolve ] = fsolve( funchFJ, vecX0, options_fsolve )
		return
	endif
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_broyd_scan = [];
	prm_broyd_scan.jupdateType = JUPDATE_TYPE__BROYDEN;
	prm_broyd_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_broyd_scan, datOut_broyd_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_broyd_scan );
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + STEP_TYPE__SCAN_LEV_MIN_KUPDATE ~~~ " );
	prm_broyd_scank = [];
	prm_broyd_scank.jupdateType = JUPDATE_TYPE__BROYDEN;
	prm_broyd_scank.stepType = STEP_TYPE__SCAN_LEV_MIN_KUPDATE;
	[ vecXF_broyd_scank, datOut_broyd_scank ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_broyd_scank );
	%
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_recalc_scan = [];
	prm_recalc_scan.jupdateType = JUPDATE_TYPE__RECALC;
	prm_recalc_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_recalc_scan, datOut_recalc_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_recalc_scan );
	%
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC_KUPDATE + STEP_TYPE__SCAN_LEV_MIN_KUPDATE  ~~~ " );
	prm_recalck_scank = [];
	prm_recalck_scank.jupdateType = JUPDATE_TYPE__RECALC_KUPDATE ;
	prm_recalck_scank.stepType = STEP_TYPE__SCAN_LEV_MIN_KUPDATE ;
	[ vecXF_recalck_scank, datOut_recalck_scank ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_recalck_scank );
	%
	return
	%
	%
	if (0)
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC + STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH ~~~ " );
	prm_recalc_scanfp = [];
	prm_recalc_scanfp.jupdateType = JUPDATE_TYPE__RECALC;
	prm_recalc_scanfp.stepType = STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH;
	[ vecXF_recalc_scanfp, datOut_recalc_scanfp ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_recalc_scanfp );
	endif
	%return
	%
	%
	%
	if (0)
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__LESQUJ_PRIMAL + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_lesquj_scan = [];
	prm_lesquj_scan.jupdateType = JUPDATE_TYPE__LESQUJ_PRIMAL;
	prm_lesquj_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_lesquj_scan, datOut_lesquj_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_lesquj_scan );
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + default ~~~ " );
	prm_broyd_def = [];
	prm_broyd_def.jupdateType = JUPDATE_TYPE__BROYDEN;
	[ vecXF_broyd_def, datOut_broyd_def ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_broyd_def );
	endif
	%
	%
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__REORTHONORM_POOL + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_rpool_scan = [];
	prm_rpool_scan.jupdateType = JUPDATE_TYPE__REORTHONORM_POOL;
	prm_rpool_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_rpool_scan, datOut_rpool_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_rpool_scan );
	%
	if (0)
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__REORTHONORM_POOL + STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH ~~~ " );
	prm_rpool_scanfp = [];
	prm_rpool_scanfp.jupdateType = JUPDATE_TYPE__REORTHONORM_POOL;
	prm_rpool_scanfp.stepType = STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH;
	[ vecXF_rpool_scanfp, datOut_rpool_scanfp ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_rpool_scanfp );
	endif
	%
	%
	%
	if (0)
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH ~~~ " );
	prm_broyd_scanfp = [];
	prm_broyd_scanfp.jupdateType = JUPDATE_TYPE__BROYDEN;
	prm_broyd_scanfp.stepType = STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH;
	[ vecXF_broyd_scanfp, datOut_broyd_scanfp ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_broyd_scanfp );
	endif
	%
	%
	%
	if (1)
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__NONE + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_none_scan = [];
	prm_none_scan.jupdateType = JUPDATE_TYPE__NONE;
	prm_none_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_none_scan, datOut_none_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_none_scan );
	endif
	%
	if (0)
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__NONE + STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH ~~~ " );
	prm_none_scanfp = [];
	prm_none_scanfp.jupdateType = JUPDATE_TYPE__NONE;
	prm_none_scanfp.stepType = STEP_TYPE__SCAN_LEV_MIN_FORCE_PATCH;
	[ vecXF_none_scanfp, datOut_none_scanfp ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_none_scanfp );
	endif
