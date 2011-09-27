#!/usr/bin/perl

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------

use Env;

#PG lettura dei parametri da cfg file
#PG --------------------------------
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

$BASEDir          = $User_Preferences{"BASEDir"};
$LISTOFSamples    = $User_Preferences{"LISTOFSamples"} ;
$EXEName          = $User_Preferences{"EXEName"} ;
$JOBCfgTemplate   = $User_Preferences{"JOBCfgTemplate"} ;
$INPUTSAVEPath    = $User_Preferences{"INPUTSAVEPath"} ;
$OUTPUTSAVEPath   = $User_Preferences{"OUTPUTSAVEPath"} ;
$OUTPUTFILEName   = $User_Preferences{"OUTPUTFILEName"} ;
$JOBModulo        = $User_Preferences{"JOBModulo"} ;



print "BASEDir = "          .$BASEDir."\n" ;
print "LISTOFSamples = "    .$LISTOFSamples."\n" ;
print "EXEName = "          .$EXEName."\n" ;
print "JOBCfgTemplate = "   .$JOBCfgTemplate."\n" ;
print "INPUTSAVEPath = "    .$INPUTSAVEPath."\n" ;
print "OUTPUTSAVEPath = "   .$OUTPUTSAVEPath."\n" ;
print "OUTPUTFILEName = "   .$OUTPUTFILEName."\n" ;
print "JOBModulo = "        .$JOBModulo."\n\n" ;






$sampleJobListFile = "./lancia.sh";
open(SAMPLEJOBLISTFILE, ">", $sampleJobListFile);


#####################################################
# PG prepare the array containing the root files list
#####################################################


open (LISTOFSamples,$LISTOFSamples) ;
while (<LISTOFSamples>)
{
  system("cd ".$BASEDir."\n");
  
  chomp($_);
  
  ($sample,$dummy,$dataFlag,$isDD,$mW,$crossSection,$color,$recoil,$ptHatCut) = split(" ") ;
  $subsample = substr($sample,0,1);
  if($subsample eq "#")
  {
    next;
  }
  
  system ("rm -r ".$sample."\n") ;

  print("Sample: ".$sample."\n") ;  
  system ("mkdir ".$sample."\n") ;
  
  
  
    
  
  
  $LISTOFFiles = "./list_".$sample.".txt" ;
  system ("ls -ltrh ".$INPUTSAVEPath.$sample." | grep root | grep -v \"MiBiCommonNT.root\" | grep -v \"histo_WprimeTreeAnalysis.root\" | awk '{print \$9}' > ".$LISTOFFiles."\n") ;
  
  
  
  $totNumber = 0;
  $jobNumber = 0;
  
  open (LISTOFFiles,$LISTOFFiles) ;
  while (<LISTOFFiles>)
  {
    #print "File = ".$_;
    ++$totNumber;
  }
  
  $jobNumber = int($totNumber/$JOBModulo);
  if( $totNumber%$JOBModulo != 0)
  {
    $jobNumber = $jobNumber+1;
  }
  
  print "NumberOfJobs = ".$jobNumber."\n";
  
  
  
  
  
  
  ################
  # loop over jobs 
  ################
  
  for($jobIt = 1; $jobIt <= $jobNumber; ++$jobIt)
  { 
    $currDir = `pwd` ;
    chomp ($currDir) ;
    
    $jobDir = $currDir."/".$sample."/JOB_".$jobIt ;
    system ("mkdir ".$jobDir." \n") ;
    
    $tempBjob = $jobDir."/bjob_".$jobIt.".sh" ;
    $command = "touch ".$tempBjob ;
    system ($command) ;
    $command = "chmod 777 ".$tempBjob ;
    system ($command) ;
    


    
    
    $tempo1 = "./tempo1" ;
    system ("cat ".$JOBCfgTemplate."   | sed -e s%MASS%".$mW.
                                   "%g | sed -e s%DATAFLAG%".$dataFlag.
                                   "%g | sed -e s%PTHATCUT%".$ptHatCut.
                                   "%g | sed -e s%RECOIL%".$recoil.
                                   "%g | sed -e s%CROSSSECTION%".$crossSection.
                                   "%g | sed -e s%OUTPUTFILENAME%".$OUTPUTSAVEPath.$sample."/".$OUTPUTFILEName."_".$jobIt.".root".
                                   "%g > ".$tempo1) ;
    
    
    
    $it = 0;
    $JOBLISTOFFiles = $jobDir."/list_".$sample.".txt";
    open (JOBLISTOFFiles, ">", $JOBLISTOFFiles) or die "Can't open file ".$JOBLISTOFFiles;

    open (LISTOFFiles2,$LISTOFFiles) ;
    while (<LISTOFFiles2>)
    {
      chomp; 
      s/#.*//;                # no comments
      s/^\s+//;               # no leading white
      s/\s+$//;               # no trailing white
      $file = $_ ;
      
      if( ($it >= ($jobIt - 1)*$JOBModulo) && ($it < ($jobIt)*$JOBModulo) )
      { 
	print JOBLISTOFFiles $INPUTSAVEPath."/".$sample."/".$file."\n";
      }
      ++$it;
    }
    
    
    
    $tempo2 = "./tempo2" ;    
    system ("cat ".$tempo1." | sed -e s%INPUTFILELIST%".$JOBLISTOFFiles."%g > ".$tempo2) ;
    
    
    $JOBCfgFile = $jobDir."/selections_cfg.py" ;
    system ("mv ".$tempo2." ".$JOBCfgFile) ;
    system ("rm ./tempo*") ;
    
    
    
    
    
    
    ######################
    # make job files
    ######################    
    
    open (SAMPLEJOBFILE, ">", $tempBjob) or die "Can't open file ".$tempBjob;

    $command = "#!/bin/tcsh" ;
    print SAMPLEJOBFILE $command."\n";

    #$command = "cd /gwteraz/users/deguio/CMSSW_4_2_8/src/" ;
    #print SAMPLEJOBFILE $command."\n";
    
    #$command = "cmsenv" ;
    #print SAMPLEJOBFILE $command."\n";

    #$command = "cd -" ;
    #print SAMPLEJOBFILE $command."\n";

    $command = "cd ".$BASEDir."../" ;
    print SAMPLEJOBFILE $command."\n";

    $command = "source ./scripts/setup.csh" ;
    print SAMPLEJOBFILE $command."\n";
    
    $command = "cd -" ;
    print SAMPLEJOBFILE $command."\n";

    # $command = "cd ".$jobDir ;
    # print SAMPLEJOBFILE $command."\n";
    
    $command = "mkdir ".$OUTPUTSAVEPath."/".$sample;
    print SAMPLEJOBFILE $command."\n";
    
    # $command = "unbuffer ".$EXEName." ".$JOBCfgFile." >> ".$jobDir."/out.txt" ;
    # print SAMPLEJOBFILE $command."\n";    


    $command = "unbuffer ".$EXEName." ".$JOBCfgFile." >> ".$jobDir."/out.txt" ;
    print SAMPLEJOBFILE $command."\n";

    # $command = "mkdir ".$sample;
    #print SAMPLEJOBFILE $command."\n";

    # $command = "mv  ".$OUTPUTFILEName.".root ".$sample;
    #print SAMPLEJOBFILE $command."\n";

    #$command = "scp -r ./".$OUTPUTFILEName."_".$jobIt.".root cmsmi4:".$OUTPUTSAVEPath.$sample;
    #print SAMPLEJOBFILE $command."\n";
    
    
    ############
    # submit job
    ############
    
    $command = "qsub -V -d ".$jobDir." -q production ".$tempBjob."\n" ;  
    print SAMPLEJOBLISTFILE $command."\n";
    
    #print "\n" ;
  }

  system ("rm ".$LISTOFFiles) ;
}  
