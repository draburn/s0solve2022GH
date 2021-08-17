%function [ xExt, fExt, retCode, datOut ] = findExt_sym_basic( xVals=[-1,0,1,2,3], fVals=[1,0,1,8,27], ...
%function [ xExt, fExt, retCode, datOut ] = findExt_sym_basic( xVals=[-1,0,1,2], fVals=[1,0,1,4], ...
function [ xExt, fExt, retCode, datOut ] = findExt_sym_basic( xVals=[0,1,2], fVals=[0,1,4], ...
%function [ xExt, fExt, retCode, datOut ] = findExt_sym_basic( xVals, fVals, ...
  xExt_initial=[], p_initial=2.0, prm=[], datIn=[] )
	%
	%
	% Init.
	findExt_sym_basic__init;
	thisFile = "findExt_sym_basic";
	%
	[ errFlag, rhoVals_initial, omega_initial, fExt_initial, f1_initial ] = findExt_sym__calcRes( ...
	  xExt_initial, p_initial, xVals, fVals, dVals );
	if ( errFlag )
		msg_error( verbLev, thisFile, __LINE__, "Initial guess produced invalid results." );
		retCode = RETCODE__BAD_INPUT;
		return;
	end
	if ( valLev >= VALLEV__MEDIUM )
		assert( isrealscalar(omega_initial) );
		assert( 0.0 <= omega_initial );
	end
	%
	%
	%
	xExt = xExt_initial;
	p = p_initial;
	omgea = omega_initial;
	numIter = 0;
	while (1)
		prm_findStep = mygetfield( prm, "prm_findStep", [] );
		[ xExt_new, p_new, retCode, datOut_findStep ] = findExt_sym_basic__findStep( ...
		  xExt, p, xVals, fVals, ...
		  xExtMin, xExtMax, pMin, pMax, dVals, ...
		  prm_findStep );
		if ( RETCODE__SUCCESS != retCode )
			msg_main( verbLev, thisFile, __LINE__, sprintf( "__findStep() was unsuccessful (%d).", retCode ) );
			% Leave retCode as-is (for now).
			break;
		end
		retCode = RETCODE__NOT_SET;
		%
		if ( abs(xExt_new-xExt) <= xExtTol  &&  abs(p_new-p) <= pTol )
			msg_main( verbLev, thisFile, __LINE__, "Converged-ish(?)." );
			xExt = xExt_new;
			p = p_new;
			retCode = RETCODE__SUCCESS;
			break;
		end
		if ( numIter >= maxNumIter )
			msg_main( verbLev, thisFile, __LINE__, "Failed." );
			retCode = RETCODE__IMPOSED_STOP;
			break;
		end
		%
		numIter++;
		msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
		  "   %3d;   %11.3e, %11.3e;   %11.3e, %11.3e.", ...
		  numIter, xExt_new, p_new, xExt_new-xExt, p_new-p ) );
		xExt = xExt_new;
		p = p_new;
	end
	%
	[ errFlag, rhoVals, omega, fExt, f1 ] = findExt_sym__calcRes( ...
	  xExt, p, xVals, fVals, dVals );
	datOut.dVals = dVals;
	datOut.xExt = xExt;
	datOut.fExt = fExt;
	datOut.p = p;
	datOut.f1 = f1;
	datOut.omega = omega;
	datOut.rhoVals = rhoVals;
	msg_main( verbLev, thisFile, __LINE__, sprintf( "xExt = %g.", xExt ) );
	msg_main( verbLev, thisFile, __LINE__, sprintf( "fExt = %g.", fExt ) );
	msg_main( verbLev, thisFile, __LINE__, sprintf( "p = %g.", p ) );
	msg_main( verbLev, thisFile, __LINE__, sprintf( "f1 = %g.", f1 ) );
	msg_main( verbLev, thisFile, __LINE__, sprintf( "omega = %g.", omega ) );
	%
	plot( ...
	  xVals, fExt+f1*abs(xVals-xExt).^p, 'x', 'markersize', 25, 'linewidth', 3, ...
	  xVals, fVals, 'ko', 'markersize', 25, 'linewidth', 3, ...
	  xVals, fVals, 'k+', 'markersize', 25, 'linewidth', 3 );
	grid on;
end
