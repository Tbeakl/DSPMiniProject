#!/bin/bash
echo "Recompressing the given image"
echo "Converting to pnm"
echo $(pngtopnm $1 > temp.ppm)
echo "Finished conversion"

for COMPRESSION in {1..73}
do
    echo "Compressing with all $COMPRESSION"
    echo $(cat temp.ppm | cjpeg -grayscale -qtables ./constant_matrices/constant_$COMPRESSION.txt | djpeg | cjpeg -grayscale -q $3 | djpeg | pnmtopng > $2/test_$COMPRESSION.png)
done