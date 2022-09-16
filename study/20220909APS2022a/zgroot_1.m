function grootXDatOut = zgroot_1( zcdo, probIndex, algoIndex, prm=[] );
	backup_prngStateDat = getprngstatedat();
	mydefs;
	startTime = time();
	%
	assert( isposintscalar(algoIndex) );
	assert( 1 <= algoIndex );
	assert( algoIndex <= zcdo.algoSetPrm.n );
	%
	% This could easily get broken. But, we'll cross that bridge if we come to it.
	modified_zcdo = zcdo;
	modified_zcdo.algoSetPrm.n = 1;
	modified_zcdo.algoSetPrm.s = [];
	modified_zcdo.algoSetPrm.s(1) = zcdo.algoSetPrm.s(algoIndex);
	modified_zcdo.algoSetPrm.s(1).p.verbLev = VERBLEV__DETAILS;
	modified_zcdo.algoSetPrm.s(1).p.valdLev = VALDLEV__VERY_HIGH;
	%
	grootXDatOut = zgroot_x( modified_zcdo, probIndex, prm );
return;
endfunction
