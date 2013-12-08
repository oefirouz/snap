echo "Making SNAP Push-Relabel"
make all

echo "Making BOOST Push-Relabel"
#g++ -O3 -I /usr/local/Cellar/boost/1.54.0/ boost_flow.cpp -o boost_flow


for filename in "rmflong_8_64" "rmflong_10_91" "rmflong_11_128" "rmflong_13_181" \
"rmflong_16_256" "rmflong_19_362" "rmflong_23_512" \
"rmfwide_28_5" "rmfwide_37_6" "rmfwide_49_7" "rmfwide_64_8" #"rmfwide_85_9"
#"rmfwide_111_10" "rmfwide_147_12" "rmfwide_194_14"
#"rmflong_16_256" "rmflong_19_362" "rmflong_23_512" "rmflong_30_724" \
do
    printf $filename
    printf "\t"
    ./pushrelabel Tests/${filename}.txt
    cat Tests/${filename}.txt | ./boost_flow
    printf "\n"
done

#for filename in "rmflong_8_64" "rmflong_10_91" "rmflong_11_128" "rmflong_13_181" \
#"rmflong_19_362"
#do
#    printf $filename
##    printf "\t"
#    cat Tests/${filename}.txt | ./pseudo_fifo
#done
#./pushrelabel Tests/rmflong_8_64.txt
#cat Tests/rmflong_8_64.txt | ./boost_flow
#./pushrelabel Tests/rmflong_10_91.txt
#./pushrelabel Tests/rmflong_11_128.txt
#./pushrelabel Tests/rmflong_13_181.txt
#./pushrelabel Tests/rmflong_16_256.txt
#./pushrelabel Tests/rmflong_19_362.txt
#./pushrelabel Tests/rmflong_23_512.txt

#./pushrelabel Tests/rmfwide_28_5.txt
#./pushrelabel Tests/rmfwide_37_6.txt
#./pushrelabel Tests/rmfwide_49_7.txt
#./pushrelabel Tests/rmfwide_64_8.txt
