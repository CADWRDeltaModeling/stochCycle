

tide_freq_radians_hr <- function()
{
  hours_per_cycle <- c(
    4*2.*12.420612,
    23.93447213,  #K1
    12.4206012,   #M2
    25.81933871,  #O1
    28.00622,     #2Q1
    24.06589,     #P1
    23.09848,     #J1
    22.30608,     #OO1
    12.90537,     # 2N2
    12.62600,     # NU2
    11.96723,     # K2
    12.,          #S2
    26.868350,    #Q1
    12.65834751,  #N2
    12.19162085,  #L2
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
  )

  freq_names <- c("LOW","K1","M2","O1","2Q1","P1","J1","OO1","2N2","NU2","K2",
                  "S2","Q1","N2","L2","M4",
                  "MK3","MO3","M6","MK5","M05",
                  "M1","M3","M5","M8","f","f12")
  names(hours_per_cycle)  <- freq_names
  radians_per_hour <- 2.*pi/hours_per_cycle
  radians_per_hour
}

tide_freq_radhr <- tide_freq_radians_hr()
usethis::use_data(tide_freq_radhr,overwrite=TRUE)

