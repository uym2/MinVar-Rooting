#! /bin/bash

head=$1
tail=$2
l=$3 
x=$4
rand_root=$5
true_tree=$6


if [ "$l" == "None" ]; then
	d_head=0
	d_tail=0
else
	d_head=`nw_distance $true_tree -sl $head`
	if [ "$tail" == "None" ]; then
		if [ "$head" == "$rand_root"]; then
			d_tail=$(($d_head-$l))
		else
			d_tail=$(($d_head+$l))
		fi
	fi
	d_tail=`nw_distance $true_tree -sl $tail`
fi

d2trueRoot.py $d_head $d_tail $l $x
