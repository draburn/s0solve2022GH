function [ bigL, vecG ] = funfun( x )
	echo__nargout = nargout
	bigL = (x-2.0).^2;
	vecG = 2.0*(x-2.0);
end
