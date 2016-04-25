#!/usr/bin/env bash

if (( $(bc <<< "$norms_end > 0") ))
then
	i=8
	#norms[0]=$norms_init
	norms[i+1]=$norms_end
	#while ((i>=2)); do
	#	diff=$(echo "scale=2; $norms_end-$norms_init" | bc -l)
	#	mid2=$(echo "$mid/2" | bc -l)
	#	norms[i]=$(echo "$norms_init+$mid+$mid2" | bc -l)
	#	norms[i-1]=$(echo "$norms_init+$mid" | bc -l)
	#	i=$((i-2))
	#	norms_end=$mid
	#done

	for (( c=i; c>=0; c--))
	do
		diff=$(echo "scale=2; $norms_end-$norms_init" | bc -l)
		mid=$(echo "scale=2; $diff/2" | bc -l)
		norms[c]=$(echo "scale=2; $norms_init+$mid" | bc -l)
		norms_end=$mid
	done
fi

declare -a epsilons=(0.0625 0.125 0.25 0.5 1.0)
declare -a TOLs=(0.01 0.001 0.0001 0.00001)
declare -a etas=(0.00390625 0.015625 0.0625 0.25 1.0 4.0 16.0)

. ../issvm/experiments/job_pool.sh



EPSILON=INF
#TOL=0.0001

#initialized and executing sparcifier
if [ "$1" == "optimize" ]; then
    job_pool_init 7 0
	if [ "$2" != "" ]; then
		norms=(${norms[${2}]})
	fi
	
	for NORM in "${norms[@]}"
	do
	    for TOL in "${TOLs[@]}"
	    do
            for ETA in "${etas[@]}"
            do
                if [ ! -f $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init ]; then
                    #echo "issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predited -A $NORM -A $ETA -A $EPSILON"
                    issvm_initialize -f $dataset_dir/$TRAIN_DATA -o $INIT_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED.init -k gaussian -K $K -b 1 -a $METHOD -A issvm_error/${DATASET}_SVM_SMO_BIASED_100000_train.predicted -A $NORM -A $ETA #-A $EPSILON
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
	job_pool_shutdown

    # check the $job_pool_nerrors for the number of jobs that exited non-zero
    echo "job_pool_nerrors: ${job_pool_nerrors}"
fi


if [ "$1" == "test" ]; then
	if [ "$2" != "" ]; then
		norms=(${norms[${2}]})
	fi

	for NORM in "${norms[@]}"
    do
	#NORM=${norms[${2}]}
	sum_test_error=0
	sum_sv=0
	for ((i=1;i<=10;i++)); do
		BEST_ERROR=1
		for TOL in "${TOLs[@]}"
	    do
            for ETA in "${etas[@]}"
            do
                #echo "issvm_evaluate -f $dataset_dir/$TRAIN_DATA -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS} -o $TEST_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS}_train.predited"
                OUTPUT="$(issvm_test -f $dataset_dir/${VAL_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL})"
                #echo $OUTPUT
                arrIN=(${OUTPUT})
                #echo ${arrIN[2]} $ETA
                if (( $(bc <<< "${arrIN[2]} < $BEST_ERROR") ))
                then
                    BEST_ERROR=${arrIN[2]}
                    BEST_ETA=$ETA
                    BEST_TOL=$TOL
                fi
            done
        done
		OUTPUT="$(issvm_test -f $dataset_dir/${TEST_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${BEST_ETA}_EP-${EPSILON}_BIASED_${BEST_TOL})"
		arrIN=(${OUTPUT})
		test_error=${arrIN[2]}
		#echo "${arrIN[1]} $test_error ${BEST_ETA}"
		sum_test_error=$(echo "scale=5; $test_error+$sum_test_error" | bc -l)
		sum_sv=$(echo "${arrIN[1]}+$sum_sv" | bc -l)
	done
	mean_sv=$(echo "scale=1; ${sum_sv}/10" | bc -l)
	mean_test_error=$(echo "scale=5; ${sum_test_error}/10" | bc -l)
	echo $mean_sv $mean_test_error "0 0"
	done
fi


if [ "$1" == "ttol" ]; then
	if [ "$2" != "" ]; then
		norms=(${TOLs[${2}]})
	fi

    for TOL in "${TOLs[@]}"
	    do
        for NORM in "${norms[@]}"
        do
        #NORM=${norms[${2}]}
        sum_test_error=0
        sum_sv=0
        for ((i=1;i<=10;i++)); do
            BEST_ERROR=1
                for ETA in "${etas[@]}"
                do
                    #echo "issvm_evaluate -f $dataset_dir/$TRAIN_DATA -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS} -o $TEST_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${ITERATIONS}_train.predited"
                    OUTPUT="$(issvm_test -f $dataset_dir/${VAL_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${ETA}_EP-${EPSILON}_BIASED_${TOL})"
                    #echo $OUTPUT
                    arrIN=(${OUTPUT})
                    #echo ${arrIN[2]} $ETA
                    echo "${NORM};${ETA};${i};${arrIN[1]};${arrIN[2]}" >> ${METHOD}/${DATASET}_${TOL}_validation.txt
                    if (( $(bc <<< "${arrIN[2]} < $BEST_ERROR") ))
                    then
                        BEST_ERROR=${arrIN[2]}
                        BEST_ETA=$ETA
                        BEST_TOL=$TOL
                    fi
                done

            OUTPUT="$(issvm_test -f $dataset_dir/${TEST_DATA}${i} -i $MODEL_DIR/${DATASET}_SVM_${METHOD}_NORM-${NORM}_ETA-${BEST_ETA}_EP-${EPSILON}_BIASED_${BEST_TOL})"
            arrIN=(${OUTPUT})
            test_error=${arrIN[2]}
            #echo "${arrIN[1]} $test_error ${BEST_ETA} ${BEST_TOL}"
            sum_test_error=$(echo "scale=5; $test_error+$sum_test_error" | bc -l)
            sum_sv=$(echo "${arrIN[1]}+$sum_sv" | bc -l)
        done
        mean_sv=$(echo "scale=1; ${sum_sv}/10" | bc -l)
        mean_test_error=$(echo "scale=5; ${sum_test_error}/10" | bc -l)
        #echo $mean_sv $mean_test_error "0 0"
        echo "${NORM};${mean_sv};${mean_test_error}" >> ${DATASET}_${TOL}_test.txt
        done
    done
fi
