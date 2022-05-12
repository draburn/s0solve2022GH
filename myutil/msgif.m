%  Function...
%    msgif( cond, fileName, lineNum, msgStr )
%  Overview...
%    Part of myutil module.
%    My formatted alternative to printf().
%    Adds customized information such as fileName, lineNum, and current time.
%  Input values...
%    cond: Condition under which to print message.
%    fileName: Name of file from which called, ideally __FILE__.
%    lineNum: Line number from which called, ideally __LINE__.
%    msgStr: The message content.
function msgif( cond, fileName, lineNum, msgStr )
	if ( 4 ~= nargin )
		msg( __FILE__, __LINE__, "Invalid number of arguments." );
		print_usage();
		return; % Superfluous?
	end
	if (~cond)
		return;
	end
	clockNow = clock();
	printf( "%02d:%02d:%02d.%03d [%s %4d]  %s\n", ...
	  clockNow(4), ...
	  clockNow(5), ...
	  floor(clockNow(6)), ...
	  floor(1000*(clockNow(6) - floor(clockNow(6)) )), ...
	  trimFileName(fileName), ...
	  lineNum, ...
	  msgStr );
return;
end
