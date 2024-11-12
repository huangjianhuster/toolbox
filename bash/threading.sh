# set up the maximun core to be used
MAX_THREAD=6

counter=1
while read p;
do
echo "processing: $p"

# run scripts here
# python Convert_xtc2dcd.py eq.pdb $p &

core=$(echo "$counter % $MAX_THREAD" | bc)  # # also can use: core=$(expr $counter % $MAX_THREAD)

if [ $core -eq 0 ]; 
then
wait
fi

counter=$(expr $counter + 1)

done    <   xtc.list

