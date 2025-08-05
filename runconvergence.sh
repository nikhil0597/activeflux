NCELL=$1
FILE=error.txt
rm -rf $FILE && touch $FILE
for ncell in $NCELL
do 
   echo "ncell = $ncell"
   python af.py -ic gauss -Tf 0.1 -cfl 0.01 -nc $ncell -pde linear -compute_error yes \
          -plot_freq 0 -time_scheme ssprk3 >log.txt
   tail -n 1 log.txt
   tail -n 1 log.txt >> $FILE
done
echo "Wrote file $FILE"