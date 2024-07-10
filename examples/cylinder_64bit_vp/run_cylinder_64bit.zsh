#!/usr/bin/zsh
#
echo "!  N_interior   N_neighbours     Drag" > cylinder_64bit_summary.txt
for n in 400 800 1600 3200 4800 6400 9600 12800; 
do
  for p in 9 16 25 36; 
  do
    echo "40" > input.txt
    echo "120" >> input.txt
    echo $n >> input.txt
    echo $p >> input.txt
    ./cylinder_64bit < input.txt > output.txt
    drag=`tail -1 output.txt | awk '{print $8}'`
    echo $n $p $drag >> cylinder_64bit_summary.txt
  done
done
