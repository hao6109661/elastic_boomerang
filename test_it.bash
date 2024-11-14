#! /bin/bash


make reparametrise_beam_test

rm -rf RESLT RESLT_old RESLT_new

mkdir RESLT
./reparametrise_beam_test --q 0.3
mv RESLT RESLT_new

mkdir RESLT
./reparametrise_beam_test --q 0.3 --old_version
mv RESLT RESLT_old

