#!/bin/bash
for i in `seq $1 $2`; do
   cd chunk_$i
      sbatch launch.sh
   cd ..
done
