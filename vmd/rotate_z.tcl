proc rotate_z {angle} {
    set sel [atomselect top "all"]
    set com [measure center $sel weight mass]
    set matrix [transaxis z $angle]

    $sel moveby [vecscale -1.0 $com]
    $sel move $matrix
    $sel moveby $com 
}
