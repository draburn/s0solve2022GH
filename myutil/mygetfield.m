%  Function...
%    val = mygetfield( s, fieldName, defaultVal = [] )
%  Overview...
%    Part of myutil module.
%    My well-behaved alternative to "getfield()".
%    Calls isfield(s,fieldName), and, if true, returns
%     getfield(s,fieldName). Otherwise, returns defaultVal.
%  Input values...
%    s: A structure.
%    fieldName: The name of the field within s.
%    defaultVal: The return value in case the field is not valid.
%  Output values...
%    val: The obtained value.
function val = mygetfield( s, fieldName, defaultVal = [] )
	if ( isfield(s,fieldName) )
		val = getfield(s,fieldName);
	else
		val = defaultVal;
	end
end

%!test
%!	myA = rand(3,1);
%!	myStr = "Konnichwa, sekaiyo!";
%!	foo.a1 = myA(1);
%!	foo.str = myStr;
%!	foo.foo2.a2 = myA(2);
%!	foo.foo2.a3 = myA(3);
%!	
%!	assert( mygetfield(foo,"a1") == myA(1) );
%!	assert( strcmp(mygetfield(foo,"str"),myStr) );
%!	assert( isempty(mygetfield(foo,"foobarbarganoush")) );
%!	assert( -pi == mygetfield(foo,"barbarTheElephant",-pi) );
%!	assert( i == mygetfield([],"myStr",i) );
%!
%!	assert( -10.0 == mygetfield(mygetfield(foo,"notFoo",-2.0),"a2",-10.0) );
%!	assert( myA(2) == mygetfield(mygetfield(foo,"foo2",-2.0),"a2",-10.0) );
