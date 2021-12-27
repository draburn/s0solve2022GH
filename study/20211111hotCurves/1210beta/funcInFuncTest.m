function foo = funcInFuncTest( x )
	function foo2 = funchFoo2( a, b )
		foo2 = a+b;
	end
	foo = x+1
	foo = funchFoo2(x,1)
return;
end
