#! /bin/bash

#intree=$1
#method=$2
#outtree=$3
#outinfo=$4

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -i|--input)
    intree="$2"
    shift # past argument
    ;;
    -m|--method)
    method="$2"
    shift # past argument
    ;;
    -o|--output)
    outtree="$2"
    shift # past argument
    ;;
    -f|--info)
    outinfo="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done


[ -e $outtree ] && rm $outtree
#[ -e $outinfo ] && rm $outinfo

temp_in=`mktemp tmp-XXXXX`
temp_out=`mktemp tmp-XXXXX`
temp_info=`mktemp tmp-XXXXX`


while read l; do 
	echo $l > $temp_in
	FastRoot.py -i $temp_in -m $method -o $temp_out -f $temp_info
	cat $temp_out >> $outtree
	cat $temp_info >> $outinfo
	rm $temp_in $temp_out $temp_info
done < $intree
