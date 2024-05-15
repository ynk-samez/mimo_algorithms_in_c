# Font Size Settings
set xtics font " Times new roman,15";
set ytics font " Times new roman,15"
set xlabel font "Arial,15";
set ylabel font "Arial,15";
set key font "Arial,15";

set format y "10^{%L}"
set xlabel "Eb/N0 [ dB ]";
set ylabel "BER ";
set grid xtics ytics mxtics mytics linewidth 0.5 linetype 4 linecolor 0;

set logscale y;
set border lw 1.5;
set ytics scale 2.0;
set key right top;
set key box;
#colors
set colorsequence classic

#file0
file0= "'MMSE_2x2.d'";
title0="'MMSE (2x2)'";

#file1
file1="'ZF_2x2.d'";
title1="'ZF (2x2)'";

#file2
file2="'MMSE-SIC_2x8.d'";
title2="'MMSE-SIC (2x8)'";

#styles
style1="lp ps 2 lw 1";


plot0=sprintf("%s with %s title %s", file0, style1,title0);
plot1=sprintf("%s with %s  title %s", file1, style1,title1);
plot2=sprintf("%s with %s title %s", file2, style1,title2);

plot[-10:20]\
 @plot0 , @plot1 , @plot2 ;
 
