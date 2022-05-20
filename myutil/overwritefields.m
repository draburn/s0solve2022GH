function s = overwritefields( s, s_over )
	if ( isempty(s_over) )
		return;
	endif
	fieldNames_over = fieldnames(s_over);
	numFields_over = size(fieldNames_over,1);
	for n=1:numFields_over
		s = setfield( s, fieldNames_over{n}, getfield( s_over, fieldNames_over{n} ) );
	endfor
endfunction
% Consider using "rmfield()" if there are fields you want to remove from s_over.

%!test
%!	prm_default.a = sqrt(eps);
%!	prm_default.b = 1.0E-4;
%!	prm_default.myName = "Hello world!";
%!	prm.x = @(x)( x.^2 );
%!	prm.b = 0.1;
%!	prm = overwritefields( prm_default, prm )
%!	prm = overwritefields( prm_default, [] )
