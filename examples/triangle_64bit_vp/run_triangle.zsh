#!/usr/bin/zsh
#
echo "!  N_interior   N_neighbours     S" > triangle_summary.txt
for n in 400 800 1600 3200 4800 6400 9600 12800 19200; 
do
  for p in 9 16 25; 
  do
    echo $n > input.txt
    echo $p >> input.txt
    ./triangle < input.txt > output.txt
    S=`tail -1 output.txt | awk '{print $2}'`
    echo $n $p $S >> triangle_summary.txt
  done
done
