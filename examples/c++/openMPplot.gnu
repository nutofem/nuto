set terminal postscript eps enhanced color font 'Helvetica,15'
set format x "%g"
set format y "%g"

#set termoption dashed

set key bottom right

set output 'SpeedUp.eps'
set title "Speedup"
plot \
     "./SpeedUp.txt" using 1:3 with lines lt 1 lc rgb "red" title "static load balancing, high effort in loop",\
     "./SpeedUp.txt" using 1:5 with lines lt 1 lc rgb "blue" title "dynamic load balancing, high effort in loop",\
     "./SpeedUp.txt" using 1:7 with lines lt 1 lc rgb "green" title "guided load balancing, high effort in loop",\
     "./SpeedUp.txt" using 1:9 with lines lt 2 lc rgb "red" title "static load balancing, low effort in loop",\
     "./SpeedUp.txt" using 1:11 with lines lt 2 lc rgb "blue" title "dynamic load balancing, low effort in loop",\
     "./SpeedUp.txt" using 1:13 with lines lt 2 lc rgb "green" title "guided load balancing, low effort in loop"

set key top right
set output 'Time.eps'
set title "Time"
plot \
     "./SpeedUp.txt" using 1:2 with lines lt 1 lc rgb "red" title "static load balancing, high effort in loop",\
     "./SpeedUp.txt" using 1:4 with lines lt 1 lc rgb "blue" title "dynamic load balancing, high effort in loop",\
     "./SpeedUp.txt" using 1:6 with lines lt 1 lc rgb "green" title "guided load balancing, high effort in loop",\
     "./SpeedUp.txt" using 1:8 with lines lt 2 lc rgb "red" title "static load balancing, low effort in loop",\
     "./SpeedUp.txt" using 1:10 with lines lt 2 lc rgb "blue" title "dynamic load balancing, low effort in loop",\
     "./SpeedUp.txt" using 1:12 with lines lt 2 lc rgb "green" title "guided load balancing, low effort in loop"
