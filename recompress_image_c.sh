#!/bin/bash
echo "Recompressing the given image"
echo "Converting to pnm"
echo $(pngtopnm $1 > temp.ppm)
echo "Finished conversion"

for COMPRESSION in {1..100}
do
    echo "Compressing with Q$COMPRESSION"
    echo $(cat temp.ppm | cjpeg -quality $COMPRESSION | djpeg | pnmtopng > $2/test_$COMPRESSION.png)
done