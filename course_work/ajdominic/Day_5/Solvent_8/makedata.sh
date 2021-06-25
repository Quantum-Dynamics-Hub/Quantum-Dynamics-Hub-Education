#!/bin/sh

sed -i -n -e '2,$p' temperature.out
awk '{print $1}' temperature.out > time.txt
awk '{print $2}' temperature.out > temp.txt

sed -i -n -e '2,$p' energy-ev.out
awk '{print $1}' energy-ev.out > E_time.txt
awk '{print $3}' energy-ev.out > kinetic.txt
awk '{print $5}' energy-ev.out > potential.txt
awk '{print $7}' energy-ev.out > total_E.txt
