#!/usr/bin/env bash

#Setup FACTOR if unset with 5
if [[ -z $FACTOR ]]
then
     FACTOR=5
fi


if (( $(bc <<< "$norms_end > 0") ))
then
	i=8
	norms[i+1]=$norms_end

	for (( c=i; c>=0; c--))
	do
		diff=$(echo "scale=6; $norms_end-$norms_init" | bc -l | sed 's/^\./0./')
		mid=$(echo "scale=6; ${diff}/${FACTOR}" | bc -l | sed 's/^\./0./')
		norms[c]=$(echo "scale=6; $norms_init+$mid" | bc -l | sed 's/^\./0./')
		norms_end=${mid}
#		echo ${norms[c]}
	done
fi
#unset FACTOR to don't use the same in a different dataset
FACTOR=

declare -a epsilons=(0.0625 0.125 0.25 0.5 1.0)
declare -a TOLs=(0.01 0.001 0.0001 0.00001 0.000001 0.0000001)
declare -a etas=(0.00390625 0.015625 0.0625 0.25 1.0 4.0 16.0)




usage()
{
cat << EOF
usage: $0 options

This script run the test1 or test2 over a machine.

OPTIONS:
   -h      Show this message
   -t      tolerance
   -e      epsilon
   -v      Verbose
EOF
}

OPTIND=1
OP=
while getopts "ht:T:e:E:n:N:o:v" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         t)
             TOLs=(${TOLs[${OPTARG}]})
             echo "set t: ${TOLs}"
             ;;
         T)
             TOLs=(${OPTARG})
             echo "set t: ${TOLs}"
             ;;
         e)
             epsilons=(${epsilons[${OPTARG}]})
             echo "set e: $epsilons"
             ;;
         E)
             epsilons=($OPTARG)
             echo "set e: $OPTARG"
             ;;
         n)
             norms=(${norms[${OPTARG}]})
             echo "set n: ${norms}"
             ;;
         N)
             norms=(${OPTARG})
             echo "set nv: ${norms}"
             ;;
         o)
             OP=${OPTARG}
             echo "set op: ${OPTARG}"
             ;;
         v)
             VERBOSE=1
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [[ -z $OP ]]
then
     usage
fi


#initialized and executing sparcifier
if [ "$OP" == "optimize2" ]; then
    . ../issvm/experiments/job_pool.sh
    job_pool_init 7 0
	
    for EPSILON in "${epsilons[@]}"
    do
	for TOL in "${TOLs[@]}"
	do
	    for NORM in "${norms[@]}"
	    do
            for ETA in "${etas[@]}"
            do
                if [ ! -f $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init ]; then
                if [ "$EPSILON" == "INF" ]; then
                    #echo "issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predicted -A $NORM -A $ETA"
                    issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predicted -A $NORM -A $ETA #-A $EPSILON
                else
                    #echo "issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predicted -A $NORM -A $ETA -A $EPSILON"
                    issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predicted -A $NORM -A $ETA -A $EPSILON
                fi
                fi
            done

            SECONDS=0
            start='date +%s'
            for ETA in "${etas[@]}"
            do
                if [ ! -f $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL} ]; then
                    echo "issvm_optimize2 -i $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -o $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL} -s ${METHOD}/${DATASET}_sparsifier_stats.txt -t ${TOL}"
                    job_pool_run issvm_optimize2 -i $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -o $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL} -s ${METHOD}/${DATASET}_sparsifier_stats.txt -t ${TOL}
                fi
            done
            #wait
            duration=$SECONDS
            now=$(date +'%Y-%m-%d %H:%M:%S')
#            echo " ***  [${now}] sparcifier with norm ${NORM} completed in $(($duration / 60))m:$(($duration % 60))s *** "
         done
	done
	done
	job_pool_shutdown

    # check the $job_pool_nerrors for the number of jobs that exited non-zero
    echo "job_pool_nerrors: ${job_pool_nerrors}"
fi

#initialized and executing optimize sparsifier
if [ "$OP" == "optimize" ]; then
    EPSILON=0.5
    for NORM in "${norms[@]}"
    do
        for ETA in "${etas[@]}"
        do
            issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predicted -A $NORM -A $ETA -A $EPSILON
            echo "issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predicted -A $NORM -A $ETA -A $EPSILON"
        done

        SECONDS=0
        start='date +%s'
        for ETA in "${etas[@]}"
        do
#            if [ ! -f $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL} ]; then
                echo "issvm_optimize -i $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -o $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS} -n ${ITERATIONS}"
                issvm_optimize -i $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -o $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS} -n ${ITERATIONS} &
#            fi
        done
        wait
        duration=$SECONDS
        now=$(date +'%Y-%m-%d %H:%M:%S')
        echo " ***  [${now}] sparcifier with norm ${NORM} completed in $(($duration / 60))m:$(($duration % 60))s *** "
     done
fi


if [ "$OP" == "test" ]; then
	for NORM in "${norms[@]}"
    do
	sum_test_error=0
	sum_sv=0
	for ((i=1;i<=10;i++)); do
		BEST_ERROR=1
		#for TOL in "${TOLs[@]}"
	    #do
            for ETA in "${etas[@]}"
            do
                #echo "issvm_evaluate -f $dataset_dir/$TRAIN_DATA -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS} -o $TEST_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS}_train.predited"
                OUTPUT="$(issvm_test -f $dataset_dir/${VAL_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS})"
                #echo $OUTPUT
                arrIN=(${OUTPUT})
                #echo ${arrIN[2]} $ETA
                if (( $(bc <<< "${arrIN[2]} < $BEST_ERROR") ))
                then
                    BEST_ERROR=${arrIN[2]}
                    BEST_ETA=$ETA
                    #BEST_TOL=$TOL
                fi
            done
        #done
		OUTPUT="$(issvm_test -f $dataset_dir/${TEST_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${BEST_ETA}_EP-${EPSILON}_BIASED_${ITERATIONS})"
		arrIN=(${OUTPUT})
		test_error=${arrIN[2]}
		echo "${arrIN[1]} $test_error ${BEST_ETA}"
		sum_test_error=$(echo "scale=5; $test_error+$sum_test_error" | bc -l)
		sum_sv=$(echo "${arrIN[1]}+$sum_sv" | bc -l)
	done
	mean_sv=$(echo "scale=1; ${sum_sv}/10" | bc -l)
	mean_test_error=$(echo "scale=5; ${sum_test_error}/10" | bc -l)
	echo $mean_sv $mean_test_error "0 0"
	done
fi


if [ "$OP" == "ttol" ]; then
    epsilonLen=${#epsilons[@]}
    for TOL in "${TOLs[@]}"
	    do
        for NORM in "${norms[@]}"
        do
        sum_test_error=0
        sum_sv=0
        for ((i=1;i<=10;i++)); do
            BEST_ERROR=1
            BEST_SV=0
            for EPSILON in "${epsilons[@]}"
            do
                for ETA in "${etas[@]}"
                do
                    #echo "issvm_test -f $dataset_dir/${VAL_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL}"
                    OUTPUT="$(issvm_test -f $dataset_dir/${VAL_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL})"
                    #echo $OUTPUT
                    arrIN=(${OUTPUT})
                    #echo ${arrIN[2]} $ETA
                    if (( $(bc <<< "${arrIN[2]} < $BEST_ERROR") ))
                    then
                        BEST_ERROR=${arrIN[2]}
                        BEST_ETA=${ETA}
                        BEST_EPSILON=${EPSILON}
                        BEST_SV=${arrIN[1]}
                    fi
                    if [ ${epsilonLen} == 1 ]; then
#                        echo "${METHOD}/${DATASET}_EP-${EPSILON}_TOL-${TOL}_validation.txt"
                        echo "${NORM};${ETA};${EPSILON};${i};${arrIN[1]};${arrIN[2]}" >> ${METHOD}/${DATASET}_EP-${EPSILON}_TOL-${TOL}_validation.txt
                    else
#                        echo "${METHOD}/${DATASET}_TOL-${TOL}_AGGRESSIVE_validation.txt"
                        echo "${NORM};${ETA};${EPSILON};${i};${arrIN[1]};${arrIN[2]}" >> ${METHOD}/${DATASET}_TOL-${TOL}_AGGRESSIVE_validation.txt
                    fi
                done
            done
            if [ ${epsilonLen} == 1 ]; then
                echo "${NORM};${BEST_ETA};${BEST_EPSILON};${i};${BEST_SV};${BEST_ERROR};****" >> ${METHOD}/${DATASET}_EP-${EPSILON}_TOL-${TOL}_validation.txt
            else
                echo "${NORM};${BEST_ETA};${BEST_EPSILON};${i};${BEST_SV};${BEST_ERROR};****" >> ${METHOD}/${DATASET}_TOL-${TOL}_AGGRESSIVE_validation.txt
            fi
            OUTPUT="$(issvm_test -f $dataset_dir/${TEST_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${BEST_ETA}_EP-${BEST_EPSILON}_BIASED_${TOL})"
            arrIN=(${OUTPUT})
            test_error=${arrIN[2]}
            #echo "${arrIN[1]} $test_error ${BEST_ETA} ${BEST_TOL}"
            sum_test_error=$(echo "scale=5; $test_error+$sum_test_error" | bc -l)
            sum_sv=$(echo "${arrIN[1]}+$sum_sv" | bc -l)
        done
        mean_sv=$(echo "scale=1; ${sum_sv}/10" | bc -l)
        mean_test_error=$(echo "scale=5; ${sum_test_error}/10" | bc -l)
        #echo $mean_sv $mean_test_error "0 0"

        if [ ${epsilonLen} == 1 ]; then
#            echo "${METHOD}/${DATASET}_EP-${EPSILON}_TOL-${TOL}_test.txt"
            echo "${NORM};${mean_sv};${mean_test_error}" >> ${METHOD}/${DATASET}_EP-${EPSILON}_TOL-${TOL}_test.txt
        else
#            echo "${METHOD}/${DATASET}_TOL-${TOL}_AGGRESSIVE_test.txt"
            echo "${NORM};${mean_sv};${mean_test_error}" >> ${METHOD}/${DATASET}_TOL-${TOL}_AGGRESSIVE_test.txt
        fi
        done
    done
fi
