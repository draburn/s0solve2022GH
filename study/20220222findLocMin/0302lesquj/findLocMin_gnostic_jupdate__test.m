	clear;
	commondefs;
	findLocMin_gnostic_jupdate__defs;
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
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__NONE ~~~ " );
	prm = [];
	prm.jupdateType = JUPDATE_TYPE__NONE;
	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
	vecX0 = zeros(sizeX,1);
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN ~~~ " );
	prm = [];
	prm.jupdateType = JUPDATE_TYPE__BROYDEN;
	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__LESQUJ_PRIMAL ~~~ " );
	prm = [];
	prm.jupdateType = JUPDATE_TYPE__LESQUJ_PRIMAL;
	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
	vecX0 = zeros(sizeX,1);
	%
	msg( __FILE__, __LINE__, "" );
	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC ~~~ " );
	prm = [];
	prm.jupdateType = JUPDATE_TYPE__RECALC;
	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
