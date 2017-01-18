#!/bin/tcsh

if ($#argv != 0) then
   echo "Usage: $0"
   exit -1
endif

set PROGRAM = 'DsServerMgr'
#echo PROGRAM = $PROGRAM

set APPCHK = `ps aux | grep -c $PROGRAM`
#echo APPCHK = $APPCHK

if ($APPCHK <= 1) then
   echo $PROGRAM is not running - Restarting
   /bin/rm /tmp/DsServerMgr.log
   DsServerMgr >& /tmp/DsServerMgr.log &
else
   echo $PROGRAM is alrady running
endif
