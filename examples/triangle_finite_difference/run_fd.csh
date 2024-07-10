#!/bin/tcsh
#
foreach name (`seq 50 50 400`)
echo $name > input.txt
./finite_difference < input.txt > output.txt
set S=`tail -1 output.txt | awk '{print $2}'`
set n_int=`tail -1 output.txt | awk '{print $3}'`
echo $name $S $n_int
end
