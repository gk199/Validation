#!/bin/bash
# list all files intersted in, change name in efficiency.C, run in root
FileList="LLP_mh125_pl3000 LLP_mh125_pl30000 LLP_mh250_pl500 LLP_mh250_pl1000 LLP_mh250_pl10000 LLP_mh350_pl500 LLP_mh350_pl1000 LLP_mh350_pl10000 LLP_mh1000_pl10000 LLP_mh1000_pl100000 QCD"
for mh_pl in $FileList
do
    echo $mh_pl
    sed -i "s/file_type = \".*/file_type = \"${mh_pl}\";/" efficiency_jetPt.C
    sed -i "s/rates_new_cond_.*/rates_new_cond_${mh_pl}\.root\");/" efficiency_jetPt.C
    sed -i "s/for .*/for ${mh_pl}\");/" efficiency_jetPt.C
    sed -i "s/L1Jetefficiency.*/L1Jetefficiency_${mh_pl}\.pdf\");/" efficiency_jetPt.C
    root <<EOF
.L efficiency_jetPt.C++
.x efficiency_jetPt.C
.q
EOF
done
