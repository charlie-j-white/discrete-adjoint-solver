#!/bin/bash
# chmod +x display.sh


mv ~/Downloads/TapenadeResults.zip ./DEV
unzip ./DEV/TapenadeResults.zip -d ./DEV

rm ./DEV/TapenadeResults.zip
rm ./DEV/*-all.f


echo "move stuff to ./DEV/"
