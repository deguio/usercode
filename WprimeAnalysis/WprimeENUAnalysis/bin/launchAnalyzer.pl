#!/usr/bin/perl

# ----------------------------------------------------------------------------
#      MAIN PROGRAM
# ----------------------------------------------------------------------------

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
#    next unless length;    # anything left?
    my ($var, $value) = split(/\s*=\s*/, $_, 2);
    $User_Preferences{$var} = $value;

  }


$LISTDatasets     = $User_Preferences{"LISTDatasets"} ;
$CMSSWCfgTemplate = $User_Preferences{"CMSSWCfgTemplate"} ;
$INpath           = $User_Preferences{"INpath"} ;
$INfolderDATA     = $User_Preferences{"INfolderDATA"} ;
$INfolderMC       = $User_Preferences{"INfolderMC"} ;
$OUTpath          = $User_Preferences{"OUTpath"} ;
$OUTname          = $User_Preferences{"OUTname"} ;
$evtListName      = $User_Preferences{"evtListName"} ;
$DATE             = $User_Preferences{"DATE"} ;

print "LISTDatasets     = "  .$LISTDatasets."\n" ;
print "CMSSWCfgTemplate = "  .$CMSSWCfgTemplate."\n" ;
print "INpath           = "  .$INpath."\n" ;
print "INfolderDATA     = "  .$INfolderDATA."\n" ;
print "INfolderMC       = "  .$INfolderMC."\n" ;
print "OUTpath          = "  .$OUTpath."\n" ;
print "OUTname          = "  .$OUTname."\n" ;
print "evtListName      = "  .$evtListName."\n" ;
print "DATE             = "  .$DATE."\n" ;

print "\n";

###################################################
# prepare the array containing the root files list
###################################################


open (LISTFile,$LISTDatasets) ;
while (<LISTFile>)
{
  print "\n";

  chomp($_);
#  $indataset_i =  $INpath.$_."_".$DATE."_merged.root";
  $indataset_i =  $INpath.$_."_merged.root";
  $outdataset_i = $OUTpath.$_."/".$OUTname;
  $outevtlist_i = $OUTpath.$_."/".$evtListName;

  print "File = ".$indataset_i." \n";
  print  "mkdir -p ".$OUTpath.$_."\n";
  system("mkdir -p ".$OUTpath.$_."\n");

  $tempo1 = "./tempo.py" ;
  system ("cat ".$CMSSWCfgTemplate." | sed -e s%LISTOFFILES%".$indataset_i. 
                                 "%g | sed -e s%OUTFILE%".$outdataset_i.
                                 "%g | sed -e s%EVTLIST%".$outevtlist_i.
	                         "%g > ".$tempo1) ;

  print    "WprimeTreeAnalysis_new ".$tempo1."\n";
  system  ("WprimeTreeAnalysis_new ".$tempo1."\n");
}

print "\n";
