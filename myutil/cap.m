%  Function...
%    aryOut = cap( aryIn, cLo, cHi )
%  Overview...
%    Part of myutil module.
%    Caps each value in aryIn to be >= cLo and <= cHi.
%  Input values...
%    aryIn: Input array.
%  Output values...
%    aryOut: Output array.
%    cLo: Low value cap.
%    cHi: High value cap.
function aryOut = cap( aryIn, cLo, cHi )
	assert( isrealarray(aryIn) );
	assert( isrealscalar(cLo) );
	assert( isrealscalar(cHi) );
	assert( cLo <= cHi );
	aryOut = aryIn + ((cLo-aryIn).*(cLo>aryIn)) + ((cHi-aryIn).*(cHi<aryIn));
return;
end

%!test
%!	aryIn = randn(3,3,2);
%!	cLo = -0.5;
%!	cHi = 0.5;
%!	aryOut = cap(aryIn,cLo,cHi);
%!	assert( cLo <= aryOut );
%!	assert( cHi >= aryOut );
