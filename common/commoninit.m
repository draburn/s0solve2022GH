commondefs;
retCode = RETCODE__NOT_SET;
datOut = [];
startTime = time();
%
verbLev = mygetfield( prm, "verbLev", VERBLEV__NOTIFY );
reportInterval = mygetfield( prm, "reportInterval", 3.0 );
assert( isrealscalar(verbLev) );
assert( isrealscalar(reportInterval) );
assert( 0.0 <= reportInterval );
reportTimePrev = startTime - 0.1;
%
valLev = mygetfield( prm, "valLev", VALLEV__HIGH );
assert( isrealscalar(valLev) );
%
figIndex0 = mygetfield( prm, "figIndex0", 0 );
assert( isposintscalar(figIndex0+1) );
figIndex = figIndex0;
