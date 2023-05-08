disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];
ports = [];

% figure refresh rate
plotting_step = 100;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 20e6;
frequency_domain.end = 20e9;
frequency_domain.step = 20e6;

% define sampled voltages
sampled_voltages(1).min_x = 14*dx;
sampled_voltages(1).min_y = 10*dy;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 20*dx;
sampled_voltages(1).max_y = 10*dy;
sampled_voltages(1).max_z = 3*dz;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = false;

sampled_voltages(2).min_x = 30*dx;
sampled_voltages(2).min_y = 36*dy;
sampled_voltages(2).min_z = 0;
sampled_voltages(2).max_x = 36*dx;
sampled_voltages(2).max_y = 36*dy;
sampled_voltages(2).max_z = 3*dz;
sampled_voltages(2).direction = 'zp';
sampled_voltages(2).display_plot = false;

% define sampled currents
sampled_currents(1).min_x = 14*dx;
sampled_currents(1).min_y = 10*dy;
sampled_currents(1).min_z = 3*dz;
sampled_currents(1).max_x = 20*dx;
sampled_currents(1).max_y = 10*dy;
sampled_currents(1).max_z = 3*dz;
sampled_currents(1).direction = 'yp';
sampled_currents(1).display_plot = false;

sampled_currents(2).min_x = 30*dx;
sampled_currents(2).min_y = 36*dy;
sampled_currents(2).min_z = 3*dz;
sampled_currents(2).max_x = 36*dx;
sampled_currents(2).max_y = 36*dy;
sampled_currents(2).max_z = 3*dz;
sampled_currents(2).direction = 'yn';
sampled_currents(2).display_plot = false;

% define ports
ports(1).sampled_voltage_index = 1;
ports(1).sampled_current_index = 1;
ports(1).impedance = 50;
ports(1).is_source_port = true;

ports(2).sampled_voltage_index = 2;
ports(2).sampled_current_index = 2;
ports(2).impedance = 50;
ports(2).is_source_port = false;
