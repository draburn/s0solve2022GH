%  Function...
%    msg( fileName, lineNum, msgStr )
%  Overview...
%    Part of myutil module.
%    My formatted alternative to printf().
%    Adds customized information such as fileName, lineNum, and current time.
%  Input values...
%    fileName: Name of file from which called, ideally __FILE__.
%    lineNum: Line number from which called, ideally __LINE__.
%    msgStr: The message content.
%
%  IMPORTANT NOTE!!!
%   On Octave 3.8.1 and earlier,
%   "__FILE__" is likely to cause Octave to crash;
%   Here, unfortunately, it may be necessary to use some "thisFile",
%   a la, 'thisFile = foo.m'; 'msg( thisFile, __LINE__, "Hello world!" );' instead.
function msg( fileName, lineNum, msgStr )
	clockNow = clock();
	printf( "%02d:%02d:%02d.%03d [%s.%d] %s\n", ...
	  clockNow(4), ...
	  clockNow(5), ...
	  floor(clockNow(6)), ...
	  floor(1000*(clockNow(6) - floor(clockNow(6)) )), ...
	  fileName, ...
	  lineNum, ...
	  msgStr );
return;
end

%!test
%!	thisName = "test msg";
%!	disp( "" );
%!	disp( "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" );
%!	msg( thisName, __LINE__, "Hello world!" );
%!	x = pi;
%!	msg( thisName, __LINE__, ["x was set to pi (" num2str(x) ")."] );
%!	x = e;
%!	msg( thisName, __LINE__, sprintf("x was set to e (%f).",x) );
%!	disp("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
%!	disp( "*** Is the above text correct?  ***" );
%!	disp( "" );
