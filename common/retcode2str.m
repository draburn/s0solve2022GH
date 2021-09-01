%  Function...
%    str = retcode2str( retCode )
%  Overview...
%    Part of my common module.
%    Generates a string describing the integer return code.
function str = retcode2str( retCode )
	commondefs;
	if (isempty(retCode))
		str = "<< THE RETURN CODE WAS NULL ([]). CHECK ORDER OF ARGUMENTS. >>";
		return;
	elseif (~isrealscalar(retCode))
		str = "<< THE RETURN CODE IS NOT A REAL SCALAR. CHECK ORDER OF ARGUMENTS. >>";
		return;
	elseif ( abs(retCode-round(retCode)) > eps050*abs(retCode) )
		str = sprintf( "<< THE RETURN CODE IS NOT AN INTEGER (%g). CHECK ORDER OF ARGUMENTS. >>", retCode );
		return;
	elseif ( RETCODE__SUCCESS == retCode )
		str0 = "Success";
	elseif ( RETCODE__IMPOSED_STOP == retCode )
		str0 = "Imposed Stop";
	elseif ( RETCODE__ALGORITHM_BREAKDOWN == retCode )
		str0 = "Algorithm Breakdown";
	elseif ( RETCODE__UNSPECIFIC_ERROR == retCode )
		str0 = "UNSPECIFIC ERROR";
	elseif ( RETCODE__BAD_INPUT == retCode )
		str0 = "BAD INPUT";
	elseif ( RETCODE__BAD_DEPENDENCY == retCode )
		str0 = "BAD DEPENDENCY";
	elseif ( RETCODE__INTERNAL_INCONSISTENCY == retCode )
		str0 = "INTERNAL INCONSISTENCY";
	elseif ( RETCODE__NOT_SET == retCode )
		str0 = "RETCODE NOT SET";
	else
		str0 = "UNKNOWN RETURN CODE";
	end
	str = sprintf( "\"%s\" (%d)", str0, retCode );
return;
end

%!test
%!	commondefs;
%!	thisFile = "test retcode2str";
%!	disp( "" );
%!	disp( "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "retcode2str(%d) = \"%s\".", ...
%!	  RETCODE__SUCCESS, ...
%!	  retcode2str(RETCODE__SUCCESS) ) );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "retcode2str(%d) = \"%s\".", ...
%!	  RETCODE__IMPOSED_STOP, ...
%!	  retcode2str(RETCODE__IMPOSED_STOP) ) );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "retcode2str(%d) = \"%s\".", ...
%!	  RETCODE__ALGORITHM_BREAKDOWN, ...
%!	  retcode2str(RETCODE__ALGORITHM_BREAKDOWN) ) );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "retcode2str(%d) = \"%s\".", ...
%!	  RETCODE__UNSPECIFIC_ERROR, ...
%!	  retcode2str(RETCODE__UNSPECIFIC_ERROR) ) );
%!	disp("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
%!	disp( "*** Is the above text correct?  ***" );
%!	disp( "" );
