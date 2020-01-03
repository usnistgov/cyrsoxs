VTIDIR=. #Directory of VTI Files
CSVDIR=$VTIDIR/CSV/ #Directory of CSV Files
CODEDIR=/home/maksbh/Documents/TalyFem/prsoxs/scripts #Directory of code 



for i in $( ls $VTIDIR/*.vti ); do
        filename=$(basename -- "$i")
        filename="${filename%%.*}"".csv"
        echo $i $CSVDIR$filename
        pvpython $CODEDIR/vtiToCSV.py $i $CSVDIR$filename
        echo Done
done