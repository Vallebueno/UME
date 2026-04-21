#!/bin/bash

input_file=$1
n=$2

x=input_file
x=${x/lst/lstR}
y=${x/lstR/lstS}
shuf "$input_file" | head -n "$n" > $x
sort -n ${x} > ${y}
