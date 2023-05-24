library("demic")
library("usethis")

#print(sessioninfo::session_info())
#print(lsf.str('package:demic'))
#print(help(package = demic))

### Main
args = commandArgs(trailingOnly=TRUE)

# Argument number and default settings
if(length(args)==1) {
  demic::estGrowthRate(args[1], paste(args[1],".eptr",sep=""))
} else if (length(args)==2) {
  demic::estGrowthRate(args[1], args[2])
} else if (length(args)==3) {
  demic::estGrowthRate(args[1], args[2], args[3])
} else if (length(args)==4) {
  demic::estGrowthRate(args[1], args[2], args[3])
}