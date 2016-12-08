#! /bin/bash

head=$1
tail=$2
l=$3 
x=$4
true_tree=$5

d_head=`nw_distance $true_tree -sl $head`
d_tail=`nw_distance $true_tree -sl $tail`

d2trueRoot.py $d_head $d_tail $l $x
