NCELL=$1
FILE_AVG=error_avg.txt
FILE_INT=error_int.txt
rm -f $FILE_AVG $FILE_INT
touch $FILE_AVG $FILE_INT

for ncell in $NCELL
do 
   echo "ncell = $ncell"
   python af.py -ic gauss -Tf 0.1 -cfl 0.01 -nc $ncell -pde linear -compute_error yes \
          -plot_freq 0 -time_scheme ssprk3 >log.txt
   line=$(tail -n 1 log.txt)
   echo $line

   # Assume output format: h er1 er2 er3 er4 (space-separated)
   # Extract fields using awk or cut

   h=$(echo $line | awk '{print $1}')
   er1=$(echo $line | awk '{print $2}')
   er2=$(echo $line | awk '{print $3}')
   er3=$(echo $line | awk '{print $4}')
   er4=$(echo $line | awk '{print $5}')

   echo "$h $er1 $er2" >> $FILE_AVG
   echo "$h $er3 $er4" >> $FILE_INT
done

echo "Wrote files $FILE_AVG and $FILE_INT"
