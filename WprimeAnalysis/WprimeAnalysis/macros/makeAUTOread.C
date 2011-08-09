///==== modify "AnalysisFileIn" to add entries definitions from "TreeIn" that is in the file "ntupleIn" ====
///==== search the string: +++AUTO-Reader+++ 
///==== e.g.   root makeAUTOread.C\(\"file.root\"\)

TTree* T ;

void makeAUTOread(string ntupleIn="rfio:/castor/cern.ch/user/d/deguio/Wprime/MiBiCommonPAT/Spring11/MC_26042011/WToENu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1_AODSIM_26042011/MiBiCommonNT_84_1_Veh.root", 
		  string TreeIn="MiBiCommonNTOneElectron/SimpleNtuple",
		  string AnalysisFileIn="branches.txt")
{
  using namespace std ;
    
  TFile* rootFileIn =  TFile::Open(ntupleIn.c_str(), "READ") ;
  T = (TTree*) rootFileIn -> Get(TreeIn.c_str()) ;
    
  char toDo[1000];
  sprintf(toDo,"T->Print(\"toponly\") ; >> treeContent.txt") ; 
  gROOT->ProcessLine(toDo);
  

  string filenameIn ;
  filenameIn  = "treeContent.txt" ; 
  char buffer[1000];                 

  ifstream fileIn(filenameIn.c_str());
  
  FILE *out_file;
  out_file = fopen( "insert.txt" , "w");
    
  while (! fileIn.eof() )
  {
    TString temp_name, trash, temp_type ;
    fileIn.getline(buffer,1000);
    istringstream line(buffer);
    
    line >> trash ;
    
    string s_trash = trash.Data() ;
    if ( s_trash != "branch:" ) continue ;

    line >> temp_name ;
    string s_temp_name = temp_name.Data() ;
    
    TBranch* B = (TBranch*) T->GetBranch(s_temp_name.c_str()) ;
    string s_temp_type =  B->GetClassName() ;
    
    
    int pos ;
    
    pos = s_temp_type.find("4D<double") ;  
    if (pos!=string::npos) 
      fprintf ( out_file , " std::vector<ROOT::Math::XYZTVector>* %s = reader.Get4V(\"%s\");\n", s_temp_name.c_str(), s_temp_name.c_str()) ;

    pos = s_temp_type.find("3D<double") ;  
    if (pos!=string::npos) 
      fprintf ( out_file , " std::vector<ROOT::Math::XYZVector>* %s = reader.Get3V(\"%s\");\n", s_temp_name.c_str(), s_temp_name.c_str()) ;

    pos = s_temp_type.find("r<int") ; 
    if (pos!=string::npos) 
      fprintf ( out_file , " std::vector<int>* %s = reader.GetInt(\"%s\");\n", s_temp_name.c_str(), s_temp_name.c_str()) ;
    
    pos = s_temp_type.find("r<float") ; 
    if (pos!=string::npos)
      fprintf ( out_file , " std::vector<float>* %s = reader.GetFloat(\"%s\");\n", s_temp_name.c_str(), s_temp_name.c_str()) ;

    pos = s_temp_type.find("r<double") ; 
    if (pos!=string::npos)
      fprintf ( out_file , " std::vector<double>* %s = reader.GetDouble(\"%s\");\n", s_temp_name.c_str(), s_temp_name.c_str()) ;
  }
  
  fclose(out_file) ;

  TString Command2Line = Form("rm treeContent.txt");
  gSystem->Exec(Command2Line);
  
  TString Command2Line = Form(("sed  '/\\+\\+\\+AUTO\\+\\+\\+/r insert.txt' < " + AnalysisFileIn + " > tmp.txt").c_str());
  gSystem->Exec(Command2Line);
  TString Command2Line = Form(("mv tmp.txt " + AnalysisFileIn).c_str()); 
  gSystem->Exec(Command2Line);

  TString Command2Line = Form("rm insert.txt");
  gSystem->Exec(Command2Line);

  gApplication->Terminate(0) ;


}
