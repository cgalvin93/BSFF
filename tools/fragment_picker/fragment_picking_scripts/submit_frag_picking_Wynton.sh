#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o ../job_outputs                        #-- output directory (fill in)
#$ -e ../job_outputs                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=40G                  #-- submits on nodes with enough free memory (required)
#$ -l scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
##$ -t 1                        #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section)

$1 -verbose -id inpuA -frag_sizes 9 -n_frags 200 -n_candidates 1000 $2
