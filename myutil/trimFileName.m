function fileName = trimFileName( fullFilePath )
	w = strchr( fullFilePath, '/' );
	if (isempty(w))
		fileName = fullFilePath;
	elseif (w(end)>=length(fullFilePath))
		fileName = fullFilePath;
	else
		fileName = fullFilePath( w(end)+1 : end );
	endif
return;
end

%!test
%!	fullFilePath = "/home/x/Documents/sxsolve/dever/s0solve2022/study/20220326BeyondZero/0328ReturnToZero/zlinsolf100.m"
%!	trimFileName( fullFilePath )
