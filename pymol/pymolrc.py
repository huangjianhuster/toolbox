from pymol import cmd, preset

def stylized_bs(selection="all"):

    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.show("spheres", selection)

    cmd.set("stick_radius", 0.07)
    cmd.set("sphere_scale", 0.18)
    cmd.set("sphere_scale", 0.13, "elem H")

    cmd.set("stick_quality", 50)
    cmd.set("sphere_quality", 4)

    cmd.color("gray85", "elem C")
    cmd.color("red", "elem O")
    cmd.color("slate", "elem N")
    cmd.color("gray98", "elem H")

    cmd.set("stick_color", "black")

    cmd.bg_color("white")

    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_texture", 2)
    cmd.set("antialias", 3)
    cmd.set("ambient", 0.5)
    cmd.set("spec_count", 5)
    cmd.set("shininess", 50)
    cmd.set("specular", 1)
    cmd.set("reflect", 0.1)

    preset.ball_and_stick(selection)

cmd.extend("stylized_bs", stylized_bs)
# example: stylized_bs org
# <-- show ball and stick for organic molecules

def suitup(selection="sele"):
    """
    Stylized protein cartoon with transparent surface
    and a highlighted region.

    Usage:
        suitup
        suitup sele
        suitup chain A
    """

    # Clean up
    cmd.remove("solvent")
    cmd.remove("metals")

    # Background and rendering
    cmd.bg_color("white")
    cmd.set("depth_cue", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_transparency_oblique", 0.1)
    cmd.set("ray_trace_depth_factor", 1)
    cmd.set("ray_trace_disco_factor", 1)
    cmd.set("ray_trace_gain", 0.008)
    cmd.set("ray_shadow", 0)

    # Cartoon representation
    cmd.show("cartoon", "all")
    cmd.cartoon("automatic", "all")

    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("cartoon_tube_radius", 0.8)
    cmd.set("cartoon_oval_length", 1.5)
    cmd.set("cartoon_oval_width", 0.2)
    cmd.set("cartoon_loop_radius", 0.2)

    # Custom colors (RGB: 0-255)
    cmd.set_color("custom_bluegray", [138, 159, 178])
    cmd.set_color("custom_red", [159, 16, 18])

    # Base coloring
    cmd.color("custom_bluegray", "all")
    cmd.color("custom_red", selection)

    # Cartoon-specific colors
    cmd.set("cartoon_color", "custom_bluegray", "all")
    cmd.set("cartoon_color", "custom_red", selection)

    # Surface representation
    cmd.show("surface", "all")
    cmd.set("surface_mode", 3)
    cmd.set("surface_color", "custom_bluegray")
    cmd.set("transparency", 0.9)

    # Lighting and shading
    cmd.set("ambient_occlusion_mode", 0)
    cmd.set("ray_transparency_contrast", 8)
    cmd.set("light_count", 6)
    cmd.set("ambient", 0.25)
    cmd.set("direct", 6)
    cmd.set("reflect", 0.0)
    cmd.set("specular", 0)

    # Render
    cmd.ray()


# Register PyMOL command
cmd.extend("suitup", suitup)
# example: suitup resi 1-10


