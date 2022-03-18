function [ matJ, datOut ] = findZero_baseline__jupdate( vecX, vecF, vecX_prev, vecF_prev, matJ_prev, funchF, prm );
	datOut = [];
	datOut.fevalCount = 0;
	jUpdateType = mygetfield( prm, "jUpdateType", "full" );
	switch (jUpdateType)
	case {"none"}
		% Nothing to do.
	case {"broyden"}
		error( "Not implemented" );
	case {"pool"}
		error( "Not implemented" );
	case {"full"}
		cfdj_prm = mygetfield( prm, "cdfj_prm", [] );
		cfdj_prm.vecF0 = vecF;
		[ matJ, cfdj_datOut ] = calcFDJ( vecX, funchF, cfdj_prm );
		datOut.fevalCount += cfdj_datOut.fevalCount;
	otherwise
		error( "Invalid value of jUpdateType." );
	endswitch
return;
endfunction
