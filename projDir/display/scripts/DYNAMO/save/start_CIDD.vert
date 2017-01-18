#! /bin/csh

cd $PROJ_DIR/display/params

setenv DATA_HOST pgen

CIDD -p CIDD.vert -i vert  |& \
    LogFilter -d $ERRORS_LOG_DIR -p CIDD -i vert &

