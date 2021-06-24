#!/bin/sh

awk '{print $1}' temperature.out > time.txt
awk '{print $2}' temperature.out > temp.txt
