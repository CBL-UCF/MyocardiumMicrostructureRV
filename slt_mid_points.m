function mid_point_list = slt_mid_points(RUC_w,RUC_t,Slt_w, Slt_t)

mid_point_list = {};

Ns = round(RUC_w/(Slt_w*1.1))+2;
Nn = round(RUC_t/(Slt_t*1.1))+2;

for pts = 1:Ns
    for ptn = 1:Nn
        y=(pts+0.25+0.5*rand())*Slt_w;
        z=(ptn+0.45+0.1*rand())*Slt_t;
        point = [y,z];
        mid_point_list{end+1} = point;
    end
end
        