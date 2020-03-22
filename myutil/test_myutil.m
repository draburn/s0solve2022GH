thisFile = "test_myutil";
figure(1);
title( "This figure intentionally left blank." );
%
msg( thisFile, __LINE__, "Testing borderzoom..." );
test borderzoom;
msg( thisFile, __LINE__, "Finished testing borderzoom.\n" );
%
msg( thisFile, __LINE__, "Testing cap..." );
test cap;
msg( thisFile, __LINE__, "Finished testing cap.\n" );
%
msg( thisFile, __LINE__, "Testing cent..." );
test cent;
msg( thisFile, __LINE__, "Finished testing cent.\n" );
%
msg( thisFile, __LINE__, "Testing fleq..." );
test fleq;
msg( thisFile, __LINE__, "Finished testing fleq.\n" );
%
msg( thisFile, __LINE__, "Testing isposintarray..." );
test isposintarray;
msg( thisFile, __LINE__, "Finished testing isposintarray.\n" );
%
msg( thisFile, __LINE__, "Testing isposintscalar..." );
test isposintscalar;
msg( thisFile, __LINE__, "Finished testing isposintscalar.\n" );
%
msg( thisFile, __LINE__, "Testing isrealarray..." );
test isrealarray;
msg( thisFile, __LINE__, "Finished testing isrealarray.\n" );
%
msg( thisFile, __LINE__, "Testing isrealscalar..." );
test isrealscalar;
msg( thisFile, __LINE__, "Finished testing isrealscalar.\n" );
%
msg( thisFile, __LINE__, "Testing issize..." );
test issize;
msg( thisFile, __LINE__, "Finished testing issize.\n" );
%
msg( thisFile, __LINE__, "Testing msg..." );
test msg;
msg( thisFile, __LINE__, "Finished testing msg.\n" );
%
msg( thisFile, __LINE__, "Testing mygetfield..." );
test mygetfield;
msg( thisFile, __LINE__, "Finished testing mygetfield.\n" );
%
msg( thisFile, __LINE__, "Testing stopsignalpresent..." );
test stopsignalpresent;
msg( thisFile, __LINE__, "Finished testing stopsignalpresent.\n" );
