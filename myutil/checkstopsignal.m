%  Function...
%    foundStopSignal = checkstopsignal( stopSignalFilePath = "stopsig" )
%  Overview...
%    Part of myutil module.
%    Checks for the file stopSignalFilePath and returns true if it exists,
%     otherwise returns false.
function foundStopSignal = checkstopsignal( stopSignalFilePath = "stopsig" )
	fid = fopen( stopSignalFilePath );
	if ( -1 == fid )
		foundStopSignal = false;
	else
		fclose( fid );
		foundStopSignal = true;
	end
return;
end

%!test
%!	if ( checkstopsignal() )
%!		printf( "DID find stop signal.\n" );
%!	else
%!		printf( "Did NOT find stop signal.\n" );
%!	end
%!	disp( "*** IS ABOVE RESULT CORRECT??? ***" );
