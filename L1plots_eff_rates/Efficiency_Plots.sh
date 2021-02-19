#!/bin/bash
# list all files intersted in, change name in efficiency.C, run in root
FileList="LLP_mh125_pl3000 LLP_mh125_pl30000 LLP_mh250_pl500 LLP_mh250_pl1000 LLP_mh250_pl10000 LLP_mh350_pl500 LLP_mh350_pl1000 LLP_mh350_pl10000 LLP_mh1000_pl1000 LLP_mh1000_pl10000 QCD"
for mh_pl in $FileList
do
    echo $mh_pl
    sed -i "s/file_type = \".*/file_type = \"${mh_pl}\";/" efficiency.C
    sed -i "s/rates_new_cond_.*/rates_new_cond_${mh_pl}\.root\");/" efficiency.C
    sed -i "s/for .*/for ${mh_pl}\");/" efficiency.C
    sed -i "s/L1efficiency.*/L1efficiency_${mh_pl}\.pdf\");/" efficiency.C
    root <<EOF
.L efficiency.C++
.x efficiency.C
.q
EOF
done
