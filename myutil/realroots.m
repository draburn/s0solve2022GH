% rts = realroots( C )
% Returns real results of built-in roots(C).
function rts = realroots( C )
	rtsFull = roots(C);
	rts = [];
	for n=1:length(rtsFull)
	if (isreal(rtsFull(n)))
		rts = [ rts, rtsFull(n) ];
	endif
	endfor
return;
end
