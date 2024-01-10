function Vf = fibreVolfrac(con_list,fibre_r,RUC_l,RUC_w,RUC_t)  %Calculating total fibre volume fraction in the RUC acoording to connectivity table

total_fibre_length = 0;


for bp_vol = 1:length(con_list)
    for line = 1:length(con_list{bp_vol})
        total_fibre_length = total_fibre_length + norm(con_list{bp_vol}{line}{2}-con_list{bp_vol}{line}{1});
    end
end

total_fibre_volume = total_fibre_length*pi*fibre_r^2;
RUC_total_vol = RUC_l*RUC_w*RUC_t;
Vf = (total_fibre_volume/RUC_total_vol)*100;

end 




        
    

