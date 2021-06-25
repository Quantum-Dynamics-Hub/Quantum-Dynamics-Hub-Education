INPUT=cof

echo 'ngrids cutoff rel       energy 	    memory     time    grid1   	 grid2 	   grid3     grid4     grid5   	 grid6'
#for ngrids in {4,6,8,10,12,14}
#do
    ngrids=6
#    for cutoff in {50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100} 
#    do
        cutoff=1000
        for relcutoff in {50,60,70,80,90,100,110,120,130,140,150}
        do
#        relcutoff=150
        logfile=$INPUT-$ngrids-$cutoff-$relcutoff.log
        energy=`grep 'Total energy:' $logfile|awk '{print $3}'`
        mem=`grep peak $logfile|awk '{print $NF}'`
        time=`grep -A4 'SCF WAVEFUNC' $logfile|tail -n1|awk '{print $4}'`
        printf "%5d %5d %5d %15.10f %5d %10.1f" "$ngrids" "$cutoff" "$relcutoff" "$energy" "$mem" "$time"
            for ((i=1; i <= ngrids; i++))
            do
                grid=`grep 'count for grid' $logfile|sed -n "$i"p|awk '{print $5}'`
                printf "%10d" "$grid"
            done
        #done
        printf "\n"

    done
#done



    

