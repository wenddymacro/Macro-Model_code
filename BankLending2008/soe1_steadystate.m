function [ys,check] = soe1_steadystate(ys, junk)

global M_;

growth_us_ss = get_param_by_name('growth_us_ss'); 
rr_us_bar_ss = get_param_by_name('rr_us_bar_ss'); 
pietar_us_ss = get_param_by_name('pietar_us_ss');


% because of that crapy (i am sorry) Koopman&durbin stuff, which doesn't
% take into account cointagration of unit-root variables, we initialize
% all unit root variables to the first observation

data;

RR_US = rr_us_bar_ss;
RR_US_BAR = rr_us_bar_ss;
UNR_US = UNR_US(1);
UNR_US_GAP = 0;
UNR_US_BAR = UNR_US;
PIE_US     = pietar_us_ss;
PIE_US4    = pietar_us_ss;
Y_US       = 0;
LGDP_US    = LGDP_US(1);
LGDP_US_BAR = LGDP_US(1);
RS_US      = rr_us_bar_ss+pietar_us_ss;
G_US  = growth_us_ss;
LCPI_US    = LCPI_US(1);
E4_PIE_US4 = pietar_us_ss;
E1_Y_US    = 0;
E1_PIE_US = pietar_us_ss;
UNR_G_US = 0;
GROWTH_US = growth_us_ss;
GROWTH4_US = growth_us_ss;
GROWTH4_US_BAR = growth_us_ss;
BLT_US = BLT_US(1);
BLT_US_BAR = BLT_US(1);
E=0;
E2=0;


check = 0;

ys = [
RR_US 
RR_US_BAR
UNR_US 
UNR_US_GAP 
UNR_US_BAR
PIE_US 
PIE_US4 
Y_US 
LGDP_US 
LGDP_US_BAR 
RS_US 
G_US 
LCPI_US 
E4_PIE_US4 
E1_Y_US 
E1_PIE_US 
UNR_G_US
GROWTH_US
GROWTH4_US
GROWTH4_US_BAR
BLT_US
BLT_US_BAR
E
E2
    ];
