#!/bin/sh

mkdir corpus
cd corpus/

wget http://pizzachili.dcc.uchile.cl/texts/code/sources.gz
wget http://pizzachili.dcc.uchile.cl/texts/dna/dna.gz
wget http://pizzachili.dcc.uchile.cl/repcorpus/artificial/fib41.gz
wget http://pizzachili.dcc.uchile.cl/repcorpus/real/kernel.gz

gunzip *.gz

cd -
