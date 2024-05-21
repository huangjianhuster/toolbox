proc change_transparency {new_alpha} {
        # This will always get the correct colors even if VMD
        # gains new definitions in the future
        set color_start [colorinfo num]
        set color_end [expr $color_start * 2]
        # Go through the list of colors (by index) and 
        # change their transp. value
        for {set color $color_start} {$color < $color_end} {incr color} {
                color change alpha $color $new_alpha
        }
}

user add menu Transp
for {set i 0} {$i < 10} {incr i} {
   set f [expr $i / 10.0]
   user add subitem Transp $f "change_transparency $f"
}

