	clear;
	commondefs;
	findLocMin_gnostic_jupdate2__defs;
	setprngstates(0);
	numFigs = 0;
	%
	function [ vecF, matJ ] = funcFJ_cubyDiagTest( vecX, c )
		sizeX = size(vecX,1);
		vecXE = (1:sizeX)';
		vecF = (vecX-vecXE) + c*(vecX-vecXE).^3;
		if ( nargout >= 2 )
			matJ = eye(sizeX) + c*3.0*diag((vecX-vecXE).^2);
		endif
	endfunction
	%
	caseNum = 100;
	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
	switch (caseNum)
	case 0
		sizeX = 2;
		c_cuby = 0.0;
	case 10
		sizeX = 15;
		c_cuby = 0.0;
	case 20
		sizeX = 2;
		c_cuby = 0.01;
	case 25
		sizeX = 15;
		c_cuby = 0.001;
	case 30
		sizeX = 15;
		c_cuby = 0.01;
	case 40
		sizeX = 15;
		c_cuby = 0.1;
	case 100
		sizeX = 15;
		c_cuby = 1.0;
	case 200
		sizeX = 20;
		c_cuby = 1.0;
	otherwise
		error( "Ivalid caseNum." );
	endswitch
	funchFJ = @(dummyX)( funcFJ_cubyDiagTest( dummyX, c_cuby ) );
	vecX0 = zeros(sizeX,1);
	sizeF = sizeX;
	%
	if (0)
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__LESQUJ_PRIMAL + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_lesquj_scan = [];
	prm_lesquj_scan.jupdateType = JUPDATE_TYPE__LESQUJ_PRIMAL;
	prm_lesquj_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_lesquj_scan, datOut_lesquj_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_lesquj_scan );
	return
	end
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__REORTHONORM_POOL + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_rpool_scan = [];
	prm_rpool_scan.jupdateType = JUPDATE_TYPE__REORTHONORM_POOL;
	prm_rpool_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_rpool_scan, datOut_rpool_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_rpool_scan );
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + default ~~~ " );
	prm_broyd_def = [];
	prm_broyd_def.jupdateType = JUPDATE_TYPE__BROYDEN;
	[ vecXF_broyd_def, datOut_broyd_def ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_broyd_def );
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_broyd_scan = [];
	prm_broyd_scan.jupdateType = JUPDATE_TYPE__BROYDEN;
	prm_broyd_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_broyd_scan, datOut_broyd_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_broyd_scan );
	%
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC + STEP_TYPE__SCAN_LEV_MIN ~~~ " );
	prm_recalc_scan = [];
	prm_recalc_scan.jupdateType = JUPDATE_TYPE__RECALC;
	prm_recalc_scan.stepType = STEP_TYPE__SCAN_LEV_MIN;
	[ vecXF_recalc_scan, datOut_recalc_scan ] = findLocMin_gnostic_jupdate2( vecX0, funchFJ, prm_recalc_scan );
