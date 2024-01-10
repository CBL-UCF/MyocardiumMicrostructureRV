clc; clear;

%% Inputs

% Dimensions of RUC

global RUC_w RUC_t RUC_l % L_SLT W_SLT T_SLT min_Y min_Z

RUC_w = 1000;  
RUC_t = 1000;
RUC_l = 500;

% Approximate sheetlet width(250 microm) and thickness(100 microm) used to
% create voronoi patterns - Sheetlet_voronoi.m

Sheetlet_voronoi

sheetlet_face1 = Shrinked_after_inf;
stlt_count = length(sheetlet_face1);

for stlt = 1:stlt_count
    stlt_v = sheetlet_face1{1,stlt};
    min_Y = min(stlt_v(:,1));
    min_Z = min(stlt_v(:,2));
    max_Y = max(stlt_v(:,1));
    max_Z = max(stlt_v(:,2));

    SLT_L = RUC_l;        %length in micarometers(Along fibre dirrection)
    SLT_W = max_Y-min_Y;  %Width in micrometers (Along Sheetlet dirrection)
    SLT_T = max_Z-min_Z;  %Thickness in micrometers (Along normal dirrection)
    
    Sheetlet_innitial
end