#!/bin/bash

echo $1"namelist.input"
SOURCE_PATH=$1"/namelist.input"
DEST_PATH='.'
cp $SOURCE_PATH $DEST_PATH
#find KWN code
LOC_KWN=which 'KWN-Deform'
cp $SOURCE_PATH $DEST_PATH
#execute KWN code
$LOC_KWN
#rm namelist.input 
