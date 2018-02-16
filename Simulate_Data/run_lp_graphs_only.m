

graphs_only = 1
%remote_to_server = 1

xsplit_graph_meanweight = [10,80]
xlim_graph_meanweight = [-10,100]

xsplit_graph_thetag = [-70,30]
xlim_graph_thetag = [-90,50]

run_lp_cs_1thetac_basicmoments
run_lp_cs_3thetacs_basicmoments
run_lp_cs_9thetacs_basicmoments


xsplit_graph_meanweight = [35,55]
xlim_graph_meanweight = [-10,100]


run_lp_cs_1thetac_interactedmoments
run_lp_cs_3thetacs_interactedmoments
run_lp_cs_9thetacs_interactedmoments

make_power_table

use_zero_cutoff = 0;
make_excess_length_table;
use_zero_cutoff =1;
make_excess_length_table;

make_size_table
