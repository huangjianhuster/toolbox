bg_color white
set antialias_shader, 2
set surface_quality, 2
set cartoon_sampling, 50
set cartoon_discrete_colors, on

set cartoon_oval_length, 1.0
set cartoon_fancy_sheets, 1

set max_ups, 0
set use_shaders, 1
set cartoon_use_shader, 1
set cgo_use_shader, 1
set dash_use_shader, 1
set dot_use_shader, 1
set line_use_shader, 1
set mesh_use_shader, 1
set nb_spheres_use_shader, 1
set nonbonded_use_shader, 1
set ribbon_use_shader, 1
set sphere_use_shader, 1
set stick_use_shader, 1
set surface_use_shader, 1


alias macaron, set_color macaron1, [250,191,209]; set_color macaron2, [250,209,179]; set_color macaron3, [250,230,173]; set_color macaron4, [199,240,209]; set_color macaron5, [191,217,250]; spectrum count, macaron1 macaron2 macaron3 macaron4 macaron5; set specular, 0; set reflect, 0; set ambient, 1; set depth_cue, 0; set ray_shadow, off; set ray_trace_mode, 1; set ray_trace_depth_factor, 1; set ray_trace_disco_factor, 1; set bg_rgb, [0.98,0.98,0.95]; ray



alias macaron2, remove solvent; remove inorganic; bg_color white; set depth_cue, 0; set ray_opaque_background, 1; set ray_shadow, 0; set antialias, 2; set ray_trace_mode, 1; set ray_trace_gain, 0.005; set ray_trace_depth_factor, 1; set ray_trace_disco_factor, 1; set cartoon_flat_sheets, 1; set cartoon_tube_radius, 0.8; set cartoon_oval_length, 1.5; set cartoon_oval_width, 0.2; set cartoon_loop_radius, 0.2; set_color color1, [0.69, 0.64, 0.82]; set_color color2, [0.60, 0.78, 0.88]; set_color color3, [0.65, 0.82, 0.71]; set_color color4, [0.87, 0.86, 0.55]; set_color color5, [0.98, 0.87, 0.55]; set_color color6, [0.97, 0.69, 0.49]; set_color color7, [0.96, 0.64, 0.49]; spectrum count, color1 color2 color3 color4 color5 color6 color7, polymer; set ambient, 1.0; set direct, 0.0; set reflect, 0.0; set specular, 0; set specular_intensity, 0; set spec_reflect, 0; set shininess, 0; ray


alias gradient, set_color grad1, [194,219,245]; set_color grad2, [158,191,235]; set_color grad3, [204,166,224]; spectrum count, grad1 grad2 grad3; set ambient, 0.85; set direct, 0.15; set specular, 0.1; set reflect, 0.5; set depth_cue, 1; set ray_shadow, on; set ray_trace_mode, 0; set antialias, 2

alias surfocus, show surface, polymer; set transparency, 0.6; set surface_color, gray90; color green, polymer; show sticks, ligand; set stick_radius, 0.3; color cyan, ligand and elem C

