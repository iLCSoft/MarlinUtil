#include "MCTree.h"



using namespace lcio ;
using namespace std;
using namespace IMPL;
using namespace EVENT;




// ##########################################
// #####                                #####
// #####   Constructor and Destructor   #####
// #####                                #####
// ##########################################

//=============================================================================

MCTree::MCTree(LCCollection* col): _col(col) {}


//=============================================================================

MCTree::~MCTree() {


}

//=============================================================================






// ##########################################
// #####                                #####
// #####        public methods          #####
// #####                                #####
// ##########################################

//=============================================================================

void MCTree::print(int opt) {

  if (_col->getTypeName() == LCIO::MCPARTICLE ) {
    
    int nParticles =  _col->getNumberOfElements() ;

    struct ispis{
      int level;
      int index;
      int line ;
    };
    // fill map with particle pointers and collection indices
    typedef std::map< MCParticle*, int > PointerToIndexMap ;
    PointerToIndexMap p2i_map ;
    std::vector<MCParticle*> moms ;
    PointerToIndexMap id2line_map;
    vector <int > number;
    // maximal size of PDG name
    unsigned int max_size_of_pdg=0;
    
    std::vector<ispis> test(nParticles+1);
    
    for( int k=0; k<nParticles; k++)
      {
	unsigned int  size_of_pdg=0;
	MCParticle* part =  dynamic_cast<MCParticle*>( _col->getElementAt( k ) ) ;
	p2i_map[ part ] = k ;
	int kakogod =part->getPDG();
	size_of_pdg=length_of_int(kakogod);
	if ( size_of_pdg>max_size_of_pdg) max_size_of_pdg=size_of_pdg;
	test[k].index=0;
	test[k].level=0;
	moms.push_back( part ) ;
	number.push_back(k);
      }
    
    for ( int i=0 ; i< nParticles;++i) {  test[i].line=i;}
   
    int lastlevel=0;
    int max_index=0;
    int current_line=0;
    
    // loop over collection - preserve order
    for(  int index = 0 ; index < nParticles ; index++)
      {
#ifdef CLHEP
	MCParticle4V part( _col->getElementAt( index ) ) ;
#else
	MCParticle* part =  dynamic_cast<MCParticle*>( _col->getElementAt( index ) ) ;
#endif
	for(unsigned int k=0;k<part->getParents().size();k++)
	  {

	    if(k>0)
	      {
		test[index].level=lastlevel;
		test[index].index=0;
		//   if ( index!=nParticles-1){
		test[index].line=current_line+1;
		current_line++;
		// }else{
		
		
	      }else {
	      
	      if ( index==(p2i_map[(part->getParents()[k])->getDaughters()[0]]))
		{
		  test[index].index=test[p2i_map[(part->getParents()[k])]].index+1;
		  test[index].level=test[p2i_map[(part->getParents()[k])]].level;
		  lastlevel= test[index].level;
		  test[index].line=test[p2i_map[(part->getParents()[k])]].line+1;
		  // test[index].line=current_line+1;
		  current_line++;
		  
		  if ( test[index].index>max_index) max_index=test[index].index;
		}else{
		test[index].index=test[p2i_map[(part->getParents()[k])]].index+1;
		if ( test[index].index>max_index) max_index=test[index].index;
		test[index].level=lastlevel+1;
		lastlevel++;
		for(unsigned int ii=0; ii< (part->getParents()[0])->getDaughters().size(); ++ii)
		  {
		    if ( index==p2i_map[((part->getParents()[0])->getDaughters()[ii])])
		      {
			
			test[index].line=test[p2i_map[(part->getParents()[k])]].line+ii;
			break;
		      }
		    
		  }
		current_line++;
	      }
	    }
	  }
      }// over n particles
    
    std::vector<int> new_line(nParticles);
    for ( int i=0 ; i< nParticles;++i) new_line[i]=0;
    
    
    current_line=0;	 stack <int> lines;
       for ( int i=0 ; i< nParticles;++i)
	 {
	   
	   
	   if( test[i].index==0 && moms[i]->getDaughters().size()==0)
	     {
	       new_line[i]=current_line;
	       current_line++;
             }
	   if (test[i].index==0 && moms[i]->getDaughters().size()!=0)
	     {
	       int kkk=0;
	       lines.push(i);
	       
	       while( !lines.empty())
		 {
		   
		   kkk=lines.top();
		   
		   
		   lines.pop();
		   
		   new_line[kkk]=current_line;
		   current_line++;
		   if ( (moms[kkk])->getDaughters().size() > 1)
		     {
		       for(unsigned int ii=0; ii< (moms[kkk])->getDaughters().size();++ii)
			 {
			   lines.push(p2i_map[((moms[kkk])->getDaughters()[ii])]);
			   
			 }
		     }
		   if ( (moms[kkk])->getDaughters().size()== 1)
		     {
		       lines.push(p2i_map[((moms[kkk])->getDaughters()[0])]);
		       
		     }
		   
		 }
	     }
	   
	 }
       
       
       
       
       
       current_line=0;
       
       // strings for output of the MC particle tree
       string  decayed(" -x +");
       string  lives_on(" -x-> ");
       string  daughter_of(" \\-> ");
       string  interaction(" x--> ");
       string  died(" --+ ");
       string  blank(" ");

       
       
       string  maximal_pdg_length;
       string max_blank;
       
       cout << endl
	    << "------------------------------------------------------------------ "
	    << endl << " Number of  MC particles = "<<nParticles << endl;
       int tmp_shift_head= max_size_of_pdg*(max_index+1)+(5*(max_index))+max_index+6;
       switch(opt)
	 {
	 case -1:
	   
	   cout << endl;
	   break;
	 case 0:
	   if ( max_size_of_pdg>=3)
	     {
	       cout << adjust_position1(tmp_shift_head) << "| "<< adjust_position1(max_size_of_pdg-3)
                    << "PDG | GS|CIS| B |INE| DT| DC| LD| S | "<< endl;
	     }else{
	     cout << adjust_position1(tmp_shift_head) << "|"
		  << "PDG| GS|CIS| B |INE| DT| DC| LD| S | "<< endl;
	   }
	   break;
	 case 1:
	   if ( max_size_of_pdg>=3)
	     {
	       cout << adjust_position1(tmp_shift_head) << "| "<< adjust_position1(max_size_of_pdg-3)
                    << "PDG |       Px     |       Py     |       Pz     |    Energy   |     Mass    |"
                    << endl;
	     }else{
	     cout << adjust_position1(tmp_shift_head) << "|"
		  << "PDG |       Px     |       Py     |       Pz     |    Energy   |     Mass    |"
                    << endl;
	   }
	   break;
	 case 2:
	   
	   if ( max_size_of_pdg>=3)
	     {
	       cout << adjust_position1(tmp_shift_head) << "| "
		    << adjust_position1(max_size_of_pdg-3)
                    << "PDG |    Energy   |     Mass    |  "<< endl;
	     }else{
	     cout << adjust_position1(tmp_shift_head) << "| "
		  << "PDG |    Energy   |     Mass    |  "<< endl;
	   }
	   break;
	 default:
	   cout << adjust_position1(tmp_shift_head) << "| Free your mind support LDC !" << endl;
	   break;
	 }
       
       
       while(maximal_pdg_length.size()<max_size_of_pdg)
	 {
	   maximal_pdg_length+=" ";
	   }
       
       int  novi_red=0 ;
       

       
       for ( int index=0; index<nParticles; ++index)
	 {
	   if (index !=0) cout << endl;
// ispisvanje redukovanog sadrzaja mc_particlea
	   int index_n=0;
	   while(new_line[index_n]!=index)
	     {
		     index_n++;
	     }
	   
	   bool first_in_row=true;
	   novi_red=index_n-1;
	      if ( first_in_row && test[index_n].index==0)
                {
		  int pdgi= moms[index_n]->getPDG();
		  string pdg= pdg_to_string(pdgi,max_size_of_pdg);
		  int tmp_shift=5;
		  if(moms[index_n]->getDaughters().size()!=0)
                  {
		    if( !(moms[index_n]->getDaughters()[0])->vertexIsNotEndpointOfParent())
		      {
                        if ( max_index!=0) tmp_shift= max_size_of_pdg*(max_index-1)+(5*(max_index-1))+max_index+tmp_shift-1;
			cout << pdg << interaction << pdg
			     << adjust_position1(tmp_shift) ;
			printShortMCInfo(moms[index_n],max_size_of_pdg,opt);
		    }else{
		      if ( max_index!=0) tmp_shift= max_size_of_pdg*(max_index)+(5*(max_index-1))+max_index+tmp_shift;
		      cout << pdg << decayed
			   << adjust_position1(tmp_shift);
		      printShortMCInfo(moms[index_n],max_size_of_pdg,opt);
		    }
		  }else{
		    if ( max_index!=0)   tmp_shift= max_size_of_pdg*(max_index)+(5*(max_index-1))+max_index+tmp_shift;
		    cout << pdg << died
			 << adjust_position1(tmp_shift);
		    printShortMCInfo(moms[index_n],max_size_of_pdg,opt);
		  }
		  
		  first_in_row=false;
		  if ( test[index_n].level!=novi_red  || index>=nParticles ) continue;
 		}    

              if(first_in_row && test[index_n].index!=0 ) 
                {
		  int pdgi= moms[index_n]->getPDG();
                    string pdg= pdg_to_string(pdgi,max_size_of_pdg);
		    int tmp_shift=5;
		    
		    cout << adjust_position(max_size_of_pdg , test[index_n].index);
		    if(moms[index_n]->getDaughters().size()!=0)
                        {
                          if( !(moms[index_n]->getDaughters()[0])->vertexIsNotEndpointOfParent())
                            {
			      if ( max_index!=0 && test[index_n].index!=max_index)   tmp_shift= max_size_of_pdg*(max_index-test[index_n].index-1)+(5*(max_index-test[index_n].index))+max_index-test[index_n].index+tmp_shift-6;
			      cout<< daughter_of << pdg << interaction<< pdg
				<< adjust_position1(tmp_shift);
			      printShortMCInfo(moms[index_n],max_size_of_pdg,opt);
			    }else{
			    if ( max_index!=0 &&test[index_n].index!=max_index )   tmp_shift= max_size_of_pdg*(max_index-test[index_n].index)+(5*(max_index-test[index_n].index))+max_index-test[index_n].index+tmp_shift-5;
			    cout<< daughter_of << pdg << decayed
				<< adjust_position1(tmp_shift);
			    printShortMCInfo(moms[index_n],max_size_of_pdg,opt);
			  }
			}else{
		      if ( max_index!=0 && test[index_n].index!=max_index)   tmp_shift= max_size_of_pdg*(max_index-test[index_n].index)+(5*(max_index-test[index_n].index))+max_index-test[index_n].index+tmp_shift;
			  cout << daughter_of << pdg
			       << adjust_position1(tmp_shift);
			  printShortMCInfo(moms[index_n],max_size_of_pdg,opt);
		    }
		    
		    first_in_row=false;
		      if ( test[index_n].level!=novi_red || index>=nParticles ) continue;
                }	      
	 }  
       std::cout<<std::endl;     
  }
  else  std::cout << "Not a MC Collection" << std::endl;  
}

//=============================================================================









// ##########################################
// #####                                #####
// #####        private methods         #####
// #####                                #####
// ##########################################

//=============================================================================

// int  MCTree::printShortMCInfo(const MCParticle* part){
    int  MCTree::printShortMCInfo(const MCParticle* part, const unsigned int pdg_size, const int options ){

      int precision_default ;
	switch( options ) {
       	    case -1:

	      break;
	    case 0: // print PDG  and bit feelds
            cout << " | "<< pdg_to_string( part->getPDG(),pdg_size)  << " | "
		 <<  part->getGeneratorStatus() << " | "
		 <<  part->isCreatedInSimulation() << " | "
		 <<  part->isBackscatter() << " | "
		 <<  part->vertexIsNotEndpointOfParent() << " | "
		 <<  part->isDecayedInTracker() << " | "
		 <<  part->isDecayedInCalorimeter() << " | "
		 <<  part->hasLeftDetector() << " | "
		 <<  part->isStopped() << " | ";
	    break;
           case 1: // print PDG  and 5 vector ( px py pz  E  M)
	     precision_default=cout.precision();
	       cout << setprecision(5)
                    << scientific
		    << showpos
	            << " | "<< pdg_to_string( part->getPDG(),pdg_size)  << " | "
                    <<   part->getMomentum()[0]  << " | "
		    <<   part->getMomentum()[1]  << " | "
		    <<   part->getMomentum()[2]  << " | "
                    << noshowpos
		    <<   part->getEnergy()  << " | "
                    <<   part->getMass()    << " | " ;
               cout << setprecision(precision_default) << fixed  ;
               break;
	       case 2: // print PDG  and  Energy   Mass
		 precision_default=cout.precision();
	       cout<< setprecision(5)
                   << scientific
	           << " | "<< pdg_to_string( part->getPDG(),pdg_size)  << " | "
		   <<   part->getEnergy()   << " | "
		   <<   part->getMass()     << " | " ;
	       cout<< setprecision(precision_default) << fixed  ;

               break;
	default :
          cout<< " | Free your mind support LDC !";
	  break;
	};

    return 1 ;

  }

int  MCTree::length_of_int( const int&  to_convert){
     ostringstream dtes;
      dtes<< to_convert;
      int tmp= dtes.str().size();
      return tmp;
}
string  MCTree::pdg_to_string( const int&  to_convert , const unsigned int&  max_size){
     ostringstream dtes;
      dtes<< to_convert;
      string pdg_to_string= dtes.str();
      while( pdg_to_string.size()<max_size)
      {
	  pdg_to_string=" "+pdg_to_string;
      }
      return pdg_to_string;
}

string  MCTree::adjust_position( const unsigned int&  pdg_size, const unsigned int& index){

      string  tmp;
      unsigned int  adjust=pdg_size*(index)+(int)(5*(index-(int)1))+index;


      while(! (tmp.size()==adjust) )
      {
	 tmp+=" ";
      }

      return tmp;
}
string  MCTree::adjust_position1(const  unsigned int& n_blanks){

      string  tmp;

      while(! (tmp.size()==n_blanks) )
      {
	 tmp+=" ";
      }

      return tmp;
}

//=============================================================================
