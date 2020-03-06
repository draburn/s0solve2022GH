%  Function...
%    msg_retcode( verbLev, fileName, lineNum, functionName, retCode )
%  Overview...
%    Part of my common module.
%    Produces a message about the return code of a function
%     iff verbLev is sufficiently high.
%    See commondefs, msg_thresh, and msg for more information.
function msg_retcode( verbLev, fileName, lineNum, functionName, retCode )
	commondefs;
	if ( RETCODE__SUCCESS == retCode )
		showMsg = ( verbLev >= VERBLEV__COPIOUS );
	elseif ( RETCODE__IMPOSED_STOP == retCode )
		showMsg = ( verbLev >= VERBLEV__PROGRESS );
	elseif ( RETCODE__ALGORITHM_BREAKDOWN == retCode )
		showMsg = ( verbLev >= VERBLEV__PROGRESS );
	else
		showMsg = ( verbLev >= VERBLEV__ERROR );
	end
	if (showMsg)
		msg( fileName, lineNum, sprintf( ...
		  "Function %s() returned %s.", functionName, retcode2str(retCode) ) );
	end
return;
end

%!test
%!	commondefs;
%!	thisFile = "test msg_retcode";
%!	disp( "" );
%!	disp( "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
%!	verbLev = VERBLEV__ERROR;
%!	msg( thisFile, __LINE__, "*** Set verbLev = VERBLEV__ERROR." );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__SUCCESS );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__IMPOSED_STOP );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__ALGORITHM_BREAKDOWN );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__UNSPECIFIC_ERROR );
%!	verbLev = VERBLEV__MAIN;
%!	msg( thisFile, __LINE__, "*** Set verbLev = VERBLEV__MAIN." );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__SUCCESS );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__IMPOSED_STOP );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__ALGORITHM_BREAKDOWN );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__UNSPECIFIC_ERROR );
%!	verbLev = VERBLEV__PROGRESS;
%!	msg( thisFile, __LINE__, "*** Set verbLev = VERBLEV__PROGRESS." );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__SUCCESS );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__IMPOSED_STOP );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__ALGORITHM_BREAKDOWN );
%!	msg_retcode( verbLev, thisFile, __LINE__, "dummy", RETCODE__UNSPECIFIC_ERROR );
%!	disp("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
%!	disp( "*** Is the above text correct?  ***" );
%!	disp( "" );
