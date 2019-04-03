##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            Bratwurst requirements                          ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

#options(unzip = "internal")
#Sys.setenv(TAR = "/bin/tar")

devtools::install_github('andquintero/bratwurst', ref='dev_hdsu')


writeLines("complete installation", ".snakemake/completeLibrary.txt")
