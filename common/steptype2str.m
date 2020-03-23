%  Function...
%    str = steptype2str( stepType )
%  Overview...
%    Part of my common module.
%    Generates a string describing the integer step type.
function str = steptype2str( stepType )
	commondefs;
	if ( STEPTYPE__NEWTON == stepType )
		str = "Newt";
	elseif ( STEPTYPE__PICARD == stepType )
		str = "Pic";
	elseif ( STEPTYPE__PICARD_SCALED == stepType )
		str = "PicS";
	elseif ( STEPTYPE__GRADDIR == stepType )
		str = "Grad";
	elseif ( STEPTYPE__GRADDIR_SCALED == stepType )
		str = "GradS";
	elseif ( STEPTYPE__LEVCURVE == stepType )
		str = "Lev";
	elseif ( STEPTYPE__LEVCURVE_SCALED == stepType )
		str = "LevS";
	elseif ( STEPTYPE__GRADCURVE == stepType )
		str = "GCur";
	elseif ( STEPTYPE__GRADCURVE_SCALED == stepType )
		str = "GCurS";
	else
		str = sprintf("??? (%g)",stepType);
	end
return;
end

%!test
%!	commondefs;
%!	thisFile = "test steptype2str";
%!	disp( "" );
%!	disp( "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "steptype2str(%d) = \"%s\".", ...
%!	  STEPTYPE__NEWTON, ...
%!	  steptype2str(STEPTYPE__NEWTON) ) );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "steptype2str(%d) = \"%s\".", ...
%!	  STEPTYPE__GRADDIR_SCALED, ...
%!	  steptype2str(STEPTYPE__GRADDIR_SCALED) ) );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "steptype2str(%d) = \"%s\".", ...
%!	  STEPTYPE__LEVCURVE, ...
%!	  steptype2str(STEPTYPE__LEVCURVE) ) );
%!	msg( thisFile, __LINE__, sprintf( ...
%!	  "steptype2str(%d) = \"%s\".", ...
%!	  STEPTYPE__GRADCURVE_SCALED, ...
%!	  steptype2str(STEPTYPE__GRADCURVE_SCALED) ) );
%!	disp("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
%!	disp( "*** Is the above text correct?  ***" );
%!	disp( "" );
