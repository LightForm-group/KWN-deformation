clear all
close all
%%% to plot
font='Times New Roman';
%font='Helvetica';
set(0,'DefaultAxesXGrid','on', 'DefaultAxesYGrid','on');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontName',font)
set(groot,'defaultLegendFontName',font)
set(groot,'defaultTextFontName',font)
set(0, 'DefaultLineLineWidth', 2);
%set(groot,'defaultAxesTickLabelFontName','Times New Roman')
set(0,'defaultAxesFontSize',24);
kb=physconst('boltzman')



%xw_al=0.021;
xw_fe=0.01;
xw_mo=0.04;
xw_w=0.04;
xw_co=0.13;
xw_cr=0.16;
xw_ti=0.037;
xw_nb=0.007;
xw_mg=0.022;
xw_zn=0.083;
xw_cu=0.019
xw_al=1-(xw_mg+xw_zn+xw_cu);
%molar weight elements
w_al=26.98;
w_fe=55.84;
w_mo=95.95;
w_w=183.84;
w_co=58.93;
w_cr=52.00;
w_ti=47.87;
w_nb=92.90;
w_ni=58.69;
w_zn=65.38;
w_mg=24.305;
w_cu=63.546;

x_al=xw_al/w_al/(xw_al/w_al+xw_zn/w_zn+xw_mg/w_mg+xw_cu/w_cu)
x_mg=xw_mg/w_mg/(xw_al/w_al+xw_zn/w_zn+xw_mg/w_mg+xw_cu/w_cu)
x_cu=xw_cu/w_cu/(xw_al/w_al+xw_zn/w_zn+xw_mg/w_mg+xw_cu/w_cu)
x_zn=xw_zn/w_zn/(xw_al/w_al+xw_zn/w_zn+xw_mg/w_mg+xw_cu/w_cu)
