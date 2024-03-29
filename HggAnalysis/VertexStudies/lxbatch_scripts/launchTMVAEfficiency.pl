#!/usr/bin/perl

print "reading ".$ARGV[0]."\n" ;

open (USERCONFIG,$ARGV[0]) ;

while (<USERCONFIG>)
{
    chomp; 
    s/#.*//;                # no comments
    s/^\s+//;               # no leading white
    s/\s+$//;               # no trailing white
#    next unless length;     # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences{$var} = $value;
}


$DATASETPATH      = $User_Preferences{"DATASETPATH"} ;
$CFGTEMPLATE      = $User_Preferences{"CFGTEMPLATE"} ;
$JOBDIR           = $User_Preferences{"JOBDIR"} ;
$OUTPUTDIR        = $User_Preferences{"OUTPUTDIR"} ;
$OUTPUTFILENAME   = $User_Preferences{"OUTPUTFILENAME"} ;
$JOBMODULO        = $User_Preferences{"JOBMODULO"} ;
$QUEUE            = $User_Preferences{"QUEUE"} ;



$LISTFile = "./list.txt" ;
print ("cmsLs -R ".$DATASETPATH." | grep root | awk '{print \"root://eoscms//eos/cms\" \$5}' > ".$LISTFile."\n") ;
system ("cmsLs -R ".$DATASETPATH." | grep root | awk '{print \"root://eoscms//eos/cms\" \$5}' > ".$LISTFile."\n") ;
#system ("cmsLs -R ".$DATASETPATH." | egrep '2012A|2012B|PromptReco-v1|198934_202016|202017_203002' | grep root | awk '{print \"root://eoscms//eos/cms\" \$5}' > ".$LISTFile."\n") ;
#system ("cmsLs -R ".$DATASETPATH." | grep 2012C | grep root | awk '{print \"root://eoscms//eos/cms\" \$5}' > ".$LISTFile."\n") ;
#system ("cmsLs -R ".$DATASETPATH." | grep 13Jul | grep root | awk '{print \"root://eoscms//eos/cms\" \$5}' > ".$LISTFile."\n") ;

$totNumber = 0;
$jobNumber = 0;

open (LISTFile,$LISTFile) ;
while (<LISTFile>)
{
##print "File = ".$_;
  ++$totNumber;
}

$jobNumber = int($totNumber/$JOBMODULO);
if( $totNumber%$JOBMODULO != 0)
{
  $jobNumber = $jobNumber+1;
}

print "NumberOfFiles = ".$totNumber."\n";
print "NumberOfJobs = ".$jobNumber."\n";



###############################
# create folders for job output
##############################    

print ("mkdir ".$JOBDIR."\n") ;
system ("mkdir ".$JOBDIR."\n") ;


################
# loop over jobs 
################

for($jobIt = 1; $jobIt <= $jobNumber; ++$jobIt)
{ 
  $currDir = `pwd`;
  system ("echo ".$currDir." \n");
  chomp ($currDir) ;
  system ("cd ".$currDir." \n");
  system ("mkdir ".$JOBDIR."/JOB_".$jobIt." \n") ;
  $currDir = $currDir."/".$JOBDIR."/JOB_".$jobIt ;
 
  ## set output file name in template cfg
  $tempo1 = "./tempo1" ;
  $tempfilename = $OUTPUTFILENAME."_".$jobIt.".root" ;
  system ("cat ".$CFGTEMPLATE." | sed -e s%OUTPUTFILENAME%".$tempfilename."%g > ".$tempo1."\n") ;
  

  ## make list of files
  $listOfFiles;   
  $it = 0;
  open (LISTFile2,$LISTFile) ;
  while (<LISTFile2>)
  {
      chomp; 
      s/#.*//;                # no comments
      s/^\s+//;               # no leading white
      s/\s+$//;               # no trailing white
      $file = $_ ;
      
      if( ($it >= ($jobIt - 1)*$JOBMODULO) && ($it < ($jobIt)*$JOBMODULO) )
      { 
	  $listOfFiles = $listOfFiles.$file."   " ; 
      }
      ++$it;
  }
  
  chop($listOfFiles) ;
  #  print $listOfFiles."\n" ;

  
  $LISTNAME = $currDir."/input_".$jobIt.".txt" ;
  system ("echo ".$listOfFiles." > ".$LISTNAME."\n") ;

  $tempo2 = "./tempo2" ;    
  system ("cat ".$tempo1." | sed -e s%INPUTLIST%".$LISTNAME."%g > ".$tempo2."\n") ;
  $listOfFiles = "" ;

  $tempo3 = "./tempo3" ;    
  system ("cat ".$tempo2." | sed -e s%OUTPUTDIR%".$OUTPUTDIR."%g > ".$tempo3."\n") ;
     
  print("echo ".$currDir);
  $CFGFILE = $currDir."/".$jobIt.".cfg" ;
  system ("mv ".$tempo3." ".$CFGFILE) ;
  system ("rm ./tempo*") ;
  
  $tempBjob = $currDir."/bjob_".$jobIt.".sh" ;
  $command = "touch ".$tempBjob ;
  system ($command) ;
  $command = "chmod 777 ".$tempBjob ;
  system ($command) ;
    
  $command = "cd ".$currDir ;
  system ("echo ".$command." > ".$tempBjob) ;
  
  $command = "export SCRAM_ARCH=slc5_amd64_gcc462" ;
  system ("echo ".$command." >> ".$tempBjob) ;
  
  $command = "eval \\` scramv1 runtime -sh \\`" ;
  system ("echo ".$command." >> ".$tempBjob) ;

  $command = "source ../../../scripts/setup.sh" ;
  system ("echo ".$command." >> ".$tempBjob) ;

  $command = "" ;
  system ("echo ".$command." >> ".$tempBjob) ;
  
  $command = "MyVertexAnalysis.exe ".$CFGFILE ;
  system ("echo ".$command." >> ".$tempBjob) ;
  
  
  
######################
# copy files
######################    
  
  system ("echo \\` pwd \\` ");

  $command = "mkdir ".$OUTPUTDIR; 
  system ("echo ".$command." >> ".$tempBjob) ;

  $command = "mv ./".$OUTPUTFILENAME."_".$jobIt.".root ".$OUTPUTDIR; 
  system ("echo ".$command." >> ".$tempBjob) ;
  
############
# submit job
############
  
  system ("cd ".$currDir." \n") ;
  print ("bsub -q ".$QUEUE." -cwd ".$currDir." ".$tempBjob."\n") ;
  system ("bsub -q ".$QUEUE." -cwd ".$currDir." ".$tempBjob."\n") ;
  
  print "\n" ;
}
  

 

   
