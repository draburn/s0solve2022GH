page_screen_output(0)
pkg load optim;
addpath( "myutil" );
addpath( "common" );
commondefs;
%
%my_editor = "/usr/bin/codeblocks";
my_editor = "/usr/bin/mousepad";
EDITOR( [my_editor ' %s'] ); % Moot if not USING_GUI; see options within GUI.
%
%msg("init",__LINE__,"Adding path \"study/20211111hotCurves/1210beta\".");
%addpath( "study/20211111hotCurves/1210beta" );
switch 1
case 1
msg("init",__LINE__,"Adding path \"study/20211111hotCurves/0101trine/\".");
addpath( "study/20211111hotCurves/0101trine/" );
msg("init",__LINE__,"Adding path \"study/20211111hotCurves/0101trine/numopt/\".");
addpath( "study/20211111hotCurves/0101trine/numopt/" );
msg("init",__LINE__,"Adding path \"study/20211111hotCurves/0101trine/testfunc2021/\".");
addpath( "study/20211111hotCurves/0101trine/testfunc2021/" );
case 2
msg("init",__LINE__,"Adding path \"study/20211111hotCurves/0103onward/\".");
addpath( "study/20211111hotCurves/0103onward/" );
end