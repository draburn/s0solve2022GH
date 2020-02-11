thisFile = "linsolf__finish";
%
%
msg_main( verbLev, thisFile, __LINE__, sprintf( ...
  "Final:  %4d  %6.2f  %8.2e  %s", ...
   numIter, ...
   time()-startTime, ...
   fracRes, ...
   retcode2str(retCode) )  );
%
%
if ( 0 == numIter )
	vecY = zeros(0,1);
	vecX = zeros(sizeX,1);
else
	vecY = (matR(1:numIter,1:numIter) \ vecZ(1:numIter,1)) * normB;
	vecX = matV(1:sizeX,1:numIter) * vecY(1:numIter,1);
end
%
%
datOut.funchMatAProd = funchMatAProd;
datOut.vecB = vecB;
datOut.prm = prm;
datOut.thisFile = thisFile;
datOut.startTime = startTime;
%
datOut.verbLev = verbLev;
datOut.reportInterval = reportInterval;
%
datOut.sizeF = sizeF;
datOut.sizeX = sizeX;
datOut.normB = normB;
%
datOut.fracResTol = fracResTol;
datOut.exeTimeLimit = exeTimeLimit;
datOut.numIterLimit = numIterLimit;
%
datOut.gsThresh0 = gsThresh0;
datOut.gsThresh1 = gsThresh1;
%
datOut.numIter = numIter;
datOut.matV = matV;
datOut.matW = matW;
datOut.matH = matH;
datOut.matR = matR;
datOut.vecG = vecG;
datOut.fracRes = fracRes;
%
datOut.vecY = vecY;
datOut.vecX = vecX;
%
%
return;
