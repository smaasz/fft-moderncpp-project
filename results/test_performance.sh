ALGORITHM=(1 3)
RADIXALGORITHM=3
RADIXTHRESHOLD=16
TEXT=("iterative" "FFTW")
GNUPLOTSCRIPT=" set title 'Performance: Powers-of-2';
                set title font 'Helvetica,14';
                set xlabel 'Input Length';
                set ylabel 'Time [ms]';
                set key left top;
                set style line 1 \
                linetype 1 linewidth 1 \
                pointtype 7 pointsize 1.5;
                set style line 2 \
                linetype 2 linewidth 1 \
                pointtype 1 pointsize 1.5;
                set style line 3 \
                linetype 3 linewidth 1 \
                pointtype 1 pointsize 1.5;
                set style line 4 \
                linetype 4 linewidth 1 \
                pointtype 1 pointsize 1.5;
                set logscale x 2;
                set logscale y 2;
                plot"

for i in 0 1
do
    ../testit -a ${ALGORITHM[$i]} -g $RADIXALGORITHM -r $RADIXTHRESHOLD -p ${TEXT[$i]} |
    awk 'BEGIN{OFS=" "}
         NR > 2 {print $1,$4}' > "${TEXT[$i]}"
    GNUPLOTSCRIPT="${GNUPLOTSCRIPT} '${TEXT[$i]}' title '${TEXT[$i]}' with linespoints linestyle $((i+1)),"
done
    
gnuplot -p -e "$GNUPLOTSCRIPT"

ALGORITHM=(1 3)
RADIXALGORITHM=3
RADIXTHRESHOLD=32
TEXT=("iterative" "FFTW")
GNUPLOTSCRIPT=" set title 'Performance: Non-Powers-of-2';
                set title font 'Helvetica,14';
                set xlabel 'Input Length';
                set ylabel 'Time [ms]';
                set key left top;
                set style line 1 \
                linetype 1 linewidth 1 \
                pointtype 7 pointsize 1.5;
                set style line 2 \
                linetype 2 linewidth 1 \
                pointtype 1 pointsize 1.5;
                set style line 3 \
                linetype 3 linewidth 1 \
                pointtype 1 pointsize 1.5;
                set style line 4 \
                linetype 4 linewidth 1 \
                pointtype 1 pointsize 1.5;
                set logscale x 2;
                set logscale y 2;
                plot"

for i in 0 1
do
    ../testit -a ${ALGORITHM[$i]} -g $RADIXALGORITHM -r $RADIXTHRESHOLD -n -p ${TEXT[$i]} |
    awk 'BEGIN{OFS=" "}
         NR > 2 {print $1,$4}' > "${TEXT[$i]},non-powers-of-2"
    GNUPLOTSCRIPT="${GNUPLOTSCRIPT} '${TEXT[$i]},non-powers-of-2' title '${TEXT[$i]}' with linespoints linestyle $((i+1)),"
done
    
gnuplot -p -e "$GNUPLOTSCRIPT"





