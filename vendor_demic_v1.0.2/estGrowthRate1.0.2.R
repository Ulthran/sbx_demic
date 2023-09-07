library("demic")
library("usethis")

#print(sessioninfo::session_info())
#print(lsf.str('package:demic'))
#print(help(package = demic))

### Main
args = commandArgs(trailingOnly=TRUE)

# Argument number and default settings
if(length(args)==1) {
  demic::estPTR(args[1], output=paste(args[1],".eptr",sep=""))
} else if (length(args)==2) {
  demic::estPTR(args[1], output=args[2])
} else if (length(args)==3 or length(args)==4) {
  demic::estPTR(args[1], output=args[2], max_candidate_iter=args[3])
}