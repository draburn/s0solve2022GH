	clear;
	%
	caseNum = 104;
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
	otherwise
		error( "Invalid caseNum." );
	endswitch
	%
	[ vecXF, datOut ] = findLocMin_alytJ( vecX0, funchFJ );
	%
	%
	figure();
	semilogy( ...
	  datOut.fevalCountVals, datOut.omegaVals-sumsq(testFuncPrm.vecFE)/2.0+eps, 'o-', ...
	  datOut.fevalCountVals(2:end), datOut.deltaNormVals+eps, 'x-' );
	grid on;
	xlabel( "feval count" );
	legend( "omega", "||delta||", "location", "northeast" );
