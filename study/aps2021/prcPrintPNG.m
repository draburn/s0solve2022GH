thisFile = "prcPrintPNG";
tic();
%
msg( thisFile, __LINE__, "Printing ..." );
figure(1);
print "omega.png";
figure(2);
print "remaining.png";
figure(3);
print "curve.png";
%
toc();
thisFile = [ "RETURN FROM " thisFile ];
return;
