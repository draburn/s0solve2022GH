function [ f, c, dfdx, dcdx ] = myfzero1216_repackFunc( x, funchFC_arrayForm )
	foo = funchFC_arrayForm(x);
	f = foo(1);
	c = foo(2);
	dfdx = foo(3);
	dcdx = foo(4);
return;
endfunction
