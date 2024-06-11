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

#unset key;
set key box ls 1 lc rgb "black";
set key spacing 1.5          # 行間を広げる（適宜調整）
set key width 3              # キーボックスの幅を調整（適宜調整）
#colors
set colorsequence classic



#file6
file5= "'EDC_1x1.d'";
title5="'pure  SISO'";

#file0
file0= "'MMSE-SIC_2x2.d'";
title0="'MMSE-SIC (2x2)'";

#file3
file3="'MMSE-SQRD_2x2.d'";
title3="'MMSE-SQRD (2x2)'";


#file1
file1="'ZF_2x2.d'";
title1="'ZF (2x2)'";

#file2
file2="'MLD_2x2.d'";
title2="'MLD (2x2)'";



#styles
style1="lp ps 3 lw 1";


plot0=sprintf("%s with %s title %s", file0, style1,title0);
plot1=sprintf("%s with %s  title %s", file1, style1,title1);
plot2=sprintf("%s with %s title %s", file2, style1,title2);

plot3=sprintf("%s with %s title %s", file3, style1,title3);
#plot4=sprintf("%s with %s title %s", file4, style1,title4);
plot5=sprintf("%s with %s title %s", file5, style1,title5);
plot[-10:10]\
 @plot5 ;