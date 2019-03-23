#! /bin/bash

eval `scramv1 runtime -sh`

echo "[$0] remove *~ and plots/"
rm *~
rm -rf plots

echo "[$0] run data"
./showerPlots config/config.ini strip_digis

echo "[$0] run powheg"
./showerPlots config/config_powheg.ini strip_digis_powheg

echo "[$0] run makePlot.py"
for file in `ls | grep json`
do 
    python makePlot.py $file
done

echo "[$0] run publishDir.py"
python publishDir.py plots/ index.php 

echo "[$0] publish files on EOS webpage:"
echo "[$0] https://battilan.web.cern.ch/battilan/HighPt/high_pt_showers_2017/"
rm -rf /eos/user/b/battilan/www/HighPt/high_pt_showers_2017/*
cp -r plots/* /eos/user/b/battilan/www/HighPt/high_pt_showers_2017/



