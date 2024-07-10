#!zsh
#
for e2 in 64 16 4 1 0.25 6.25e-2 1.56e-2 3.91e-3 9.77e-4;
do
  echo $e2 > input.txt
  ./poisson < input.txt >& output.txt
  rms=`tail -1 output.txt | awk '{if ($2~"epsilon2" ) {print "invalid"} else {print $2}}'`
  echo $e2 $rms
done
rm input.txt
rm output.txt
