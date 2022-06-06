function [ f, vecG ] = getL1Norm( x )
	f = norm(x,1);
	vecG = sign(x);
return;
endfunction
