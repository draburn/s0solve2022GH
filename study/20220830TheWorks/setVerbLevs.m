% Verblev likely doesn't apply to built-in error and warning messages.
%  _SILENT: No messages. (May not be fully supported.)
%  _ERROR: only fatal error messages.
%  _FLAGGED: for flagging stuff in development.
%  _WARN: include warnings.
%  _NOTIFY: include information of interest, like _FLAGGED in production.
%  _MAIN: should be limited to one or two lines of the main results.
%  _PROGRESS: include progress information to a varying extent.
%  _COPIOUS: include all progress information and worthwhile internal stuff.
%  _UNLIMITED: anything and everything.
VERBLEV__SILENT = 0;
VERBLEV__ERROR = 100;
VERBLEV__FLAGGED = 200;
VERBLEV__WARNING = 300;
VERBLEV__NOTICE = 400;
VERBLEV__MAIN = 500;
VERBLEV__PROGRESS = 600;
VERBLEV__COPIOUS = 700;
VERBLEV__UNLIMITED = 800;
% "warn" and "notify" are verbs, the rest aren't...
VERBLEV__WARN = 300;
VERBLEV__NOTIFY = 400;
%
VERBLEV__INFO = 550;
VERBLEV__PROG = 600;
