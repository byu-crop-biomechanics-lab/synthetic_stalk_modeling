function [normalized_sensitivity] = get_sensitivity(yref,ynew,pct_change)



normalized_sensitivity = ((ynew - yref)/yref)/pct_change;




end