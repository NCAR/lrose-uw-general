#! /bin/csh 


cd $PROJ_DIR/system/params

running "DsFileDist -params DsFileDist.matlab"
if ($status == 1) then
  DsFileDist -params DsFileDist.matlab |& \
	LogFilter -d $ERRORS_LOG_DIR -p DsFileDist -i matlab >& /dev/null &
endif

