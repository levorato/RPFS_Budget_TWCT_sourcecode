array1=("2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20")
array2=("2" "3" "4" "5" "6" "7" "8" "9" "10")
array3=("0.5" "0.8")
for n in "${array1[@]}"; do   # The quotes are necessary here
  for m in "${array2[@]}"; do
    for c in "${array3[@]}"; do 
        cmd="python generate_random_instances.py -n $n -m $m -c $c"
        echo $cmd
        eval "$cmd"
    done
  done
done

