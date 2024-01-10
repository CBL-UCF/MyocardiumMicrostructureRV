function [graph_size, largest_component] = largest_con_group(v1_list,v2_list)

graph_size = [];
largest_component =[];
v1_now = [];
v2_now = [];

for myo = 1: length(v1_list)
    v1_now(end+1) = v1_list(myo);
    v2_now(end+1) = v2_list(myo);
    G_now = graph(v1_now,v2_now);
    
    component_list = conncomp(G_now);
    
    [comp_name,freq] = mode(component_list);
    
    graph_size(end+1) = length(v1_now);
    largest_component(end+1) = freq;
    
end  

