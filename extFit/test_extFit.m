clear;
commondefs;
thisFile = "test_extFit";
numFigs = 0;
setprngstates(0);
%
xVals = [ 2.3306461846735673  -0.0435845986213925   0.3421158369702867  -0.0809152803594539 ];
%fVals = [ 2.6409339442016080   0.0274374833626288   0.2919875291420616   0.0558125665104530 ];
fVals = abs(xVals).^1.14771342552317;
wVals = ones(size(xVals));
%
if (0)
	[ s, p, retCode, datOut ] = extFit__internal( xVals, fVals, wVals );
	msg( thisFile, __LINE__, sprintf("extFit() returned %s.", retcode2str(retCode)) );
	[ rhoVals, bigF0, bigF1, omega ] = extFit__calcAtPt( s, p, xVals, fVals, wVals );
	omega = omega
	bigF0 = bigF0
	bigF1 = bigF1
	return
end
%
[ s, p, retCode, datOut ] = extFit( xVals, fVals, wVals );
msg( thisFile, __LINE__, sprintf("extFit() returned %s.", retcode2str(retCode)) );
omega = datOut.omega
bigF0 = datOut.bigF0
bigF1 = datOut.bigF1
