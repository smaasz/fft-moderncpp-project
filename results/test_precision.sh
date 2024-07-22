ALGORITHM=(1 2 1 2)
TEXT=("single precision, iterative" "single precision, recursive" "double precision, iterative" "double precision, recursive")
GNUPLOTSCRIPT=" set title 'Single vs. Double Precision';
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
    ../testit -a ${ALGORITHM[$i]} -g 1 -s -p ${TEXT[$i]} |
    awk 'BEGIN{OFS=" "}
         NR > 2 {print $1,$4}' > "${TEXT[$i]}"
    GNUPLOTSCRIPT="${GNUPLOTSCRIPT} '${TEXT[$i]}' title '${TEXT[$i]}' with linespoints linestyle $i,"
done

for i in 2 3
do
    ../testit -a ${ALGORITHM[$i]} -g 1 -p ${TEXT[$i]} |
    awk 'BEGIN{OFS=" "}
         NR > 2 {print $1,$4}' > "${TEXT[$i]}"
    GNUPLOTSCRIPT="${GNUPLOTSCRIPT} '${TEXT[$i]}' title '${TEXT[$i]}' with linespoints linestyle $i,"
done
    
gnuplot -p -e "$GNUPLOTSCRIPT"


