proc morph_linear {t N} {
  return [expr {double($t) / double($N)}]
}
proc morph_cycle {t N} {
  global M_PI
  return [expr {(1.0 - cos( $M_PI * double($t) / ($N + 1.0)))/2.0}]
}
proc morph_sin2 {t N} {
  global M_PI
  return [expr {sqrt(sin( $M_PI * double($t) / double($N) / 2.0))}]
}



proc morph {molid N {morph_type morph_linear}} {
    # make sure there are only two animation frames
    if {[molinfo $molid get numframes] != 2} {
	error "Molecule $molid must have 2 animation frames"
    }
    # workaround for the 'animate dup' bug; this will translate
    # 'top' to a number, if needed
    set molid [molinfo $molid get id]

    # Do some error checking on N
    if {$N != int($N)} {
	  error "Need an integer number for the number of frames"
    }
    if {$N <= 2} {
	  error "The number of frames must be greater than 2"
    }

    # Get the coordinates of the first and last frames (there are only 2)
    set sel1 [atomselect $molid "all" frame 0]
    set sel2 [atomselect $molid "all" frame 1]
    set x1 [$sel1 get x]
    set y1 [$sel1 get y]
    set z1 [$sel1 get z]
    set x2 [$sel2 get x]
    set y2 [$sel2 get y]
    set z2 [$sel2 get z]

    # Make N-2 new frames (copied from the last frame)
    for {set i 2} {$i < $N} {incr i} {
	  animate dup frame 1 $molid
    }
    # there are now N frames

    # Do the linear interpolation in steps of 1/N so
    # f(0) = 0.0 and f(N-1) = 1.0
    for {set t 0} {$t < [expr $N -1]} {incr t} {
	  # Here's the call to the user-defined morph function
	  set f [$morph_type $t $N]
	  # calculate the linear interpolation for each coordinate
	  # go to the given frame and set the coordinates
	  $sel1 frame $t
      $sel1 set x [vecadd [vecscale [expr {1.0 - $f}] $x1] [vecscale $f $x2]]
      $sel1 set y [vecadd [vecscale [expr {1.0 - $f}] $y1] [vecscale $f $y2]]
      $sel1 set z [vecadd [vecscale [expr {1.0 - $f}] $z1] [vecscale $f $z2]]
   } 
}
