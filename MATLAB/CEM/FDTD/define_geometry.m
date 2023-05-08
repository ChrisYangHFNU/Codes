disp('defining the problem geometry');

bricks = [];
spheres = [];

% define a substrate
bricks(1).min_x = 0;
bricks(1).min_y = 0;
bricks(1).min_z = 0;
bricks(1).max_x = 50*dx;
bricks(1).max_y = 46*dy;
bricks(1).max_z = 3*dz;
bricks(1).material_type = 4;

% define a PEC plate
bricks(2).min_x = 14*dx;
bricks(2).min_y = 0;
bricks(2).min_z = 3*dz;
bricks(2).max_x = 20*dx;
bricks(2).max_y = 20*dy;
bricks(2).max_z = 3*dz;
bricks(2).material_type = 2;

% define a PEC plate
bricks(3).min_x = 30*dx;
bricks(3).min_y = 26*dy;
bricks(3).min_z = 3*dz;
bricks(3).max_x = 36*dx;
bricks(3).max_y = 46*dy;
bricks(3).max_z = 3*dz;
bricks(3).material_type = 2;

% define a PEC plate
bricks(4).min_x = 0;
bricks(4).min_y = 20*dy;
bricks(4).min_z = 3*dz;
bricks(4).max_x = 50*dx;
bricks(4).max_y = 26*dy;
bricks(4).max_z = 3*dz;
bricks(4).material_type = 2;

% define a PEC plate
bricks(5).min_x = 0;
bricks(5).min_y = 0;
bricks(5).min_z = 0;
bricks(5).max_x = 50*dx;
bricks(5).max_y = 46*dy;
bricks(5).max_z = 0;
bricks(5).material_type = 2;
