## for time
set text_time {}

proc enabletrace {} {
  global vmd_frame
  global text_time
  trace variable vmd_frame([molinfo top]) w drawcounter
}

proc disabletrace {} {
  global vmd_frame
  trace vdelete vmd_frame([molinfo top]) w drawcounter
}

proc drawcounter { name element op } {
  global vmd_frame
  global text_time

  # Delete old labels
  foreach id $text_time {
      graphics top delete $id
  }
  set text_time {}

  draw color red
  set timeperframe 0.05
  set time [format "%8.2f ns" [expr $vmd_frame([molinfo top]) * $timeperframe]]
  set id [graphics top text {50 65 130} "$time" size 2 thickness 4]
  lappend text_time $id
}

enabletrace
