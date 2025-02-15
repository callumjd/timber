#!/bin/bash
gfortran ./readit.f 
grep "Atomic Center " $1 > a
grep "ESP Fit" $1 > b
grep "Fit    " $1 > c
./a.out 
rm -f a b c a.out readit.o
