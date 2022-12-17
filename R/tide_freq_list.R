

tide_freq_radians_hr <- function()
{
  hours_per_cycle <- c(
    4*2.*12.420612,
    23.93447213,  #K1
    12.4206012,
    25.81933871,
    12.,
    26.868350,    #O1
    12.65834751,
    12.19162085,  #0.0820235525 S2
    6.210300601,
    8.177140247,
    8.38630296301,
    4.140200401,
    4.930880215,
    12.4206012*4,   # half cycle
    12.4206012*2,
    12.4206012*2/3.,
    12.4206012*2/5.,
    12.4206012*2/8.,
    25.5,
    19.# 19 is 1.26cpd
  ) #0.2028035475

  freq_names <- c("LOW","K1","M2","O1","S2","Q1","N2","L2","M4",
                  "MK3","MO3","M6","MK5","M05","M1","M3","M5","M8","f","f12")
  names(hours_per_cycle)  <- freq_names
  radians_per_hour <- 2.*pi/hours_per_cycle
  radians_per_hour
}

tide_freq_radhr <- tide_freq_radians_hr()
usethis::use_data(tide_freq_radhr,overwrite=TRUE)


