%  Function...
%    msg_thresh( verbLevRequired, verbLev, fileName, lineNum, msgStr )
%  Overview...
%    Part of my common module.
%    Verbosity level limited wrapper for msg().
%    If (verbLev >= verbLevRequired), then msg( fileName, lineNum, msgStr )
%     is called. Otherwise, nothing happens.
%    See commondefs and msg for more information.
function msg_thresh( verbLevRequired, verbLev, fileName, lineNum, msgStr )
	commondefs;
	if ( verbLev >= verbLevRequired )
		msg( fileName, lineNum, msgStr );
		boo = 1;
	else
		boo = 0;
	end
return;
end

%!test
%!	commondefs;
%!	thisFile = "test msg_thresh";
%!	verbLev = VERBLEV__MAIN;
%!	disp( "" );
%!	disp( "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
%!	msg_thresh( VERBLEV__ERROR, verbLev, thisFile, __LINE__, ...
%!	  "This should be displayed." );
%!	msg_thresh( VERBLEV__MAIN, verbLev, thisFile, __LINE__, ...
%!	  "This should also be displayed." );
%!	msg_thresh( VERBLEV__PROGRESS, verbLev, thisFile, __LINE__, ...
%!	  "*** THIS SHOULD NOT BE DISPLAYED! ***" );
%!	disp("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
%!	disp( "*** Is the above text correct?  ***" );
%!	disp( "" );
