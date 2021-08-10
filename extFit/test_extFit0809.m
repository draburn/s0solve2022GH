clear;
commondefs;
thisFile = "test_extFit0809";
numFigs = 0;
setprngstates(0);
%
msg( thisFile, __LINE__, "THIS TEST CASE NEEDS ADDRESSING!" );
xVals = [ 2.3306461846735673  -0.0435845986213925   0.3421158369702867  -0.0809152803594539 ];
xVals = [ xVals, 1.0 ];
%fVals = [ 2.6409339442016080   0.0274374833626288   0.2919875291420616   0.0558125665104530 ];
fVals = abs(xVals).^1.14771342552317;
wVals = ones(size(xVals));
%
if (0)
	[ s, p, bigF0, bigF1, retCode, datOut ] = extFit__internal( xVals, fVals );
	msg( thisFile, __LINE__, sprintf("extFit() returned %s.", retcode2str(retCode)) );
	return
end
%
[ s, p, bigF0, bigF1, retCode, datOut ] = extFit( xVals, fVals );
msg( thisFile, __LINE__, sprintf("extFit() returned %s.", retcode2str(retCode)) );
