#!/bin/bash

echo $1"namelist.input"
SOURCE_PATH=$1"namelist.input"
DEST_PATH='.'
cp $SOURCE_PATH $DEST_PATH
/home/mbignon/.local/bin/KWN-Deform
#rm namelist.input
