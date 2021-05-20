clear;
commondefs;
thisFile = "test_groot1d";
tic();

verbLev = VERBLEV__COPIOUS;
caseNum = 2;

prm = []; % Unless overwritten by case.
datIn = []; % Unless overwritten by case.

switch (caseNum)
case (0)
	msg_main( verbLev, thisFile, __LINE__, ...
	  sprintf("Case: Basic test case.") );
	funchF = @(x)( x.*(10.0-x.^2) + 20.0 );
	x1 = 0.0;
	x2 = 1.0;
case (1)
	msg_main( verbLev, thisFile, __LINE__, ...
	  sprintf("Case: Simple glancing zero, 'complex behavior'.") );
	funchF = @(x)( (x-(pi-1.0)).^20.0 );
	x1 = 1.0;
	x2 = 3.0;
case (2)
	msg_main( verbLev, thisFile, __LINE__, ...
	  sprintf("Case: Very bumpy, slow converging.") );
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
	x1 = 0.0;
	x2 = 1.0;
case (3)
	msg_main( verbLev, thisFile, __LINE__, ...
	  sprintf("Case: Very bumpy, 'complex behavior'.") );
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
	x1 = 0.0;
	x2 = 5.0;
otherwise
	error(sprintf("Unsupported value of caseNum (%d).",caseNum) );
end

[ xFinal, retCode, datOut ] = groot1d( funchF, x1, x2, prm, datIn );
msg_retcode( verbLev, thisFile, __LINE__, retCode );
