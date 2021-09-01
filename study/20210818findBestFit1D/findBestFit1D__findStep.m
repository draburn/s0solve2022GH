function [ vecDelta, retCode, datOut ] = findBestFit1D__findStep( funchRho, rhoArgs, vecZ, prm )
	%
	% Init
	commondefs;
	thisFile = "findBestFit1D__findStep";
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
	%
	vecDelta = [];
	retCode = RETCODE__NOT_SET;
	datOut = [];
	%
	%
	%
	sizeZ = max(size(vecZ));
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecZ,[sizeZ,1]) );
	end
	%
	[ omega, vecG, matH, retCode, datOut_calcLocalModel ] = findBestFit1D__calcLocalModel( ...
	  funchRho, rhoArgs, vecZ, prm );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__calcLocalModel() returned %s.", retcode2str(retCode) ) );
		% Leave retCode as is.
		return;
	end
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega) );
		assert( 0.0 <= omega );
		assert( isrealarray(vecG,[sizeZ,1]) );
		assert( isrealarray(matH,[sizeZ,sizeZ]) );
	end
	datOut.datOut_calcLocalModel = datOut_calcLocalModel;
	%
	%
	% TODO: proper BT loop.
	[ vecDelta, retCode, datOut_calcDelta ] = findBestFit1D__findStep__calcDelta( omega, vecG, matH, 0.0, prm );
	if ( RETCODE__SUCCESS ~= retCode )
		msg_error( verbLev, thisFile, __LINE__, sprintf( ...
		  "__calcDelta() returned %s.", retcode2str(retCode) ) );
		% Leave retCode as is.
		return;
	end
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecDelta,[sizeZ,1]) );
	end
	datOut.datOut_calcDelta = datOut_calcDelta;
	%
	retCode = RETCODE__SUCCESS;
	return;
end
