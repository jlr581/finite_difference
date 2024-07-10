#!zsh
#
echo "With Taylor series expansion" 
for e2 in 64 16 4 1 0.25 6.25e-2 1.56250000E-02 3.90625000E-03 9.76562500E-04 2.44140625E-04 6.10351562E-05 1.52587891E-05 3.81469727E-06 9.53674316E-07 2.38418579E-07 5.96046448E-08;
do
  echo $e2 > input.txt
  echo T >> input.txt
  ./poisson < input.txt >& output.txt
  rms=`tail -1 output.txt | awk '{if ($2~"epsilon2" ) {print "invalid"} else {print $2}}'`
  echo $e2 $rms
done

echo "Without Taylor series expansion" 
for e2 in 64 16 4 1 0.25 6.25e-2 1.56250000E-02 3.90625000E-03 9.76562500E-04 2.44140625E-04 6.10351562E-05 1.52587891E-05 3.81469727E-06 9.53674316E-07 2.38418579E-07 5.96046448E-08;
do
  echo $e2 > input.txt
  echo F >> input.txt
  ./poisson < input.txt >& output.txt
  rms=`tail -1 output.txt | awk '{if ($2~"epsilon2" ) {print "invalid"} else {print $2}}'`
  echo $e2 $rms
done

rm input.txt
rm output.txt
