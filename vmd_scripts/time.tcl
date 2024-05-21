## for time


proc enabletrace {} {
  global vmd_frame;
  trace variable vmd_frame([molinfo top]) w drawcounter
}

proc disabletrace {} {
  global vmd_frame;
  trace vdelete vmd_frame([molinfo top]) w drawcounter
}

proc drawcounter { name element op } {
  global vmd_frame;

  draw delete all
  # puts "callback!"
  draw color red
  set timeperframe 0.05
  set time [format "%8.2f ns" [expr $vmd_frame([molinfo top]) * $timeperframe]]
  draw text {50 65 130} "$time" size 2 thickness 4
}

enabletrace

Type:   enabletrace

