
// Daniel Pitzl, Simon Spannagel, Simon Schnake (DESY) May-Aug 2015
// event display 4 module planes

// eudaqevd4 -n 99 000073 (for 4028, 4035, 4030, 4023)

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <stdlib.h> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>

#include <TApplication.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLine.h>

class MyMainFrame:public TGMainFrame
{
private:
  TGMainFrame * fMain;
  TRootEmbeddedCanvas *fEcanvas;
public:
  MyMainFrame( const TGWindow * p, UInt_t w, UInt_t h );
  //virtual ~ MyMainFrame(  );
  ~MyMainFrame(  );
  TCanvas *GetCanvas(  );
};

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int adc;
  double cal;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int sumA; // DP
  double charge;
  double col,row;
  bool bigx, bigy;
};

struct track {
  double xA;
  double yA;
  double xD;
  double yD;
};

// globals:

pixel pb[66560]; // global declaration: vector of pixels with hit
int fNHit; // global

//------------------------------------------------------------------------------
MyMainFrame::MyMainFrame( const TGWindow * p, UInt_t w, UInt_t h )
  :TGMainFrame( p, w, h )
{
  cout << "MyMainFrame..." << endl;
  // Create a main frame:
  fMain = new TGMainFrame( p, w, h );

  fMain->SetWMPosition( 99, 0 ); // no effect

  // Create canvas widget:
  fEcanvas = new TRootEmbeddedCanvas( "Ecanvas", fMain, w, h );

  fMain->AddFrame( fEcanvas,
                   new TGLayoutHints( kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1 ) );

  // Set a name to the main frame:
  fMain->SetWindowName( "psi46test" );

  // Map all subwindows of main frame:
  fMain->MapSubwindows(  );

  // Initialize the layout algorithm:
  fMain->Resize( fMain->GetDefaultSize(  ) );

  // Map main frame:
  fMain->MapWindow(  );
}

MyMainFrame::~MyMainFrame(  )
{
  // Clean up used widgets: frames, buttons, layouthints
  fMain->Cleanup(  );
  cout << "MyMainFrame: Cleanup" << endl;
  delete fMain;
  delete fEcanvas;
}

TCanvas *MyMainFrame::GetCanvas(  )
{
  return ( fEcanvas->GetCanvas(  ) );
}

// ----------------------------------------------------------------------
vector<cluster> getClus()
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  const int fCluCut = 1; // clustering: 1 = no gap (15.7.2012)
  //const int fCluCut = 2;

  vector<cluster> v;
  if( fNHit == 0 ) return v;

  int* gone = new int[fNHit];

  for( int i = 0; i < fNHit; ++i )
    gone[i] = 0;

  int seed = 0;

  while( seed < fNHit ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( int i = 0; i < fNHit; ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list o f clusters:

    c.sumA = 0;
    c.charge = 0;
    c.size = 0;
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    c.bigx = 0;
    c.bigy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      c.sumA += p->adc; // Aout
      double Qpix = p->cal; // calibrated [Vcal]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      if( p->col ==  0 ) c.bigx = 1;
      if( p->col == 51 ) c.bigx = 1;
      if( p->row == 79 ) c.bigy = 1;
    }

    c.size = c.vpix.size();

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( ! c.charge == 0 ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with zero charge" << endl;
    }

    v.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( (++seed < fNHit) && gone[seed] );

  } // while over seeds

  // nothing left, return clusters

  delete gone;
  return v;
}

//------------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give file name" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << "run " << run << endl;

  //FileReader reader = FileReader("000073", "run$6R$X");
  FileReader reader = FileReader(runnum.c_str(), "run$6R$X");

  // further arguments:

  int lev = 99;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-n" ) )
      lev = atoi( argv[++i] );

  } // argc
  
  // alignments:

  double alignx[4];
  double aligny[4];
  double fx[4];
  double fy[4];
  double tx[4];
  double ty[4];

  if( run >= 70 ) { // updated fo run 70

    alignx[0] = -1.162; // [mm] same sign as dxAB
    aligny[0] = -0.448; // [mm] same sign as dy
    fx[0] = -0.004; // [rad] rot, same     sign dxvsy
    fy[0] = -0.007; // [rad] rot, opposite sign dyvsx
    tx[0] = 0.0005; // from dxvsx
    ty[0] = 0.0006; // from dyvsy

    alignx[2] = -0.404; // [mm] same sign as dxCB
    aligny[2] =  0.292; // [mm] same sign as dy
    fx[2] = -0.0055; // [rad] rot, same     sign dxvsy
    fy[2] =  0.0048; // [rad] rot, opposite sign dyvsx
    tx[2] = 0.0003; // from dxvsx
    ty[2] =-0.0050; // from dyvsy, same sign

    alignx[3] =  0.288; // [mm] same sign as dxDB
    aligny[3] =  0.588; // [mm] same sign as dy
    fx[3] = -0.0015; // [rad] rot, same     sign dxvsy
    fy[3] =  0.0056; // [rad] rot, opposite sign dyvsx
    tx[3] = 0.0004; // from dxvsx
    ty[3] =-0.0024; // from dyvsy, same sign

  }

  if( run >= 98 ) { // turn 0
    
    alignx[0] =-0.172; // [mm] same sign as dxAB
    aligny[0] =-0.412; // [mm] same sign as dy
    fx[0] =-0.0050; // [rad] rot, same sign dxvsy
    fy[0] =-0.0055; // [rad] rot, opposite sign dyvsx
    tx[0] = 0.0000; // from dxvsx
    ty[0] = 0.0000; // from dyvsy
    
    alignx[2] = 0.185; // [mm] same sign as dxCB
    aligny[2] = 0.390; // [mm] same sign as dy
    fx[2] =-0.0020; // [rad] rot, same sign dxvsy
    fy[2] =-0.0023; // [rad] rot, opposite sign dyvsx
    tx[2] = 0.0000; // from dxvsx
    ty[2] = 0.0030; // from dyvsy, same sign
      
    alignx[3] = 0.228; // [mm] same sign as dxDB
    aligny[3] = 0.253; // [mm] same sign as dy
    fx[3] = 0.0035; // [rad] rot, same sign dxvsy
    fy[3] = 0.0038; // [rad] rot, opposite sign dyvsx
    tx[3] = 0.0000; // from dxvsx
    ty[3] = 0.0005; // from dyvsy, same sign
      
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile* histoFile = new TFile( "evd4.root", "RECREATE" );

  cout << "ROOT application..." << endl;

  TApplication theApp( "comet", &argc, argv );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // set ROOT styles:

  gStyle->SetTextFont( 62 ); // 62 = Helvetica bold
  gStyle->SetTextAlign( 11 );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.01, "y" );
  gStyle->SetTickLength( -0.01, "z" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.013, "y" );
  gStyle->SetLabelOffset( 0.022, "z" );

  gStyle->SetTitleOffset( 1.6, "x" );
  gStyle->SetTitleOffset( 1.6, "y" );
  gStyle->SetTitleOffset( 1.7, "z" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );
  gStyle->SetLabelFont( 62, "z" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );
  gStyle->SetTitleFont( 62, "z" );

  gStyle->SetTitleBorderSize( 0 ); // no frame around global title
  gStyle->SetTitleAlign( 13 ); // 13 = left top align
  gStyle->SetTitleX( 0.12 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetLineWidth( 1 ); // frames
  gStyle->SetHistLineColor( 4 ); // 4=blau
  gStyle->SetHistLineWidth( 3 );
  gStyle->SetHistFillColor( 5 ); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle( 1001 ); // 1001 = solid

  gStyle->SetFrameLineWidth( 2 );

  // statistics box:

  gStyle->SetOptStat( 10 );
  gStyle->SetStatFormat( "8.6g" ); // more digits, default is 6.4g
  gStyle->SetStatFont( 42 ); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize( 1 ); // no 'shadow'

  gStyle->SetStatX( 0.80 );
  gStyle->SetStatY( 0.95 );

  gStyle->SetPalette( 1 ); // rainbow colors

  gStyle->SetHistMinimumZero(  ); // no zero suppression

  gStyle->SetOptDate( 0 );

  cout << "open ROOT window..." << endl;
  MyMainFrame *myMF = new
    MyMainFrame( gClient->GetRoot(  ), 1400, 376 ); // module

  cout << "open Canvas..." << endl;
  TCanvas * c1 = myMF->GetCanvas(  );

  c1->SetBottomMargin( 0.08 );
  c1->SetLeftMargin( 0.03 );
  c1->SetRightMargin( 0.10 );

  gPad->Update(  ); // required


  //------------------------------------------------------------------------------
  //event loop:
  
  size_t event_nr = 0;
  int kev = 0;

  do {
    // Get next event:
    DetectorEvent evt = reader.GetDetectorEvent();

    if (evt.IsBORE()) { eudaq::PluginManager::Initialize(evt); }

    bool ldbg = 0;

    if (event_nr%1000==0)
      ldbg = 1;

    event_nr++;

    if( ldbg ) std::cout<<"Processing event "<< event_nr<<std::endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);
    vector <cluster> cl[4];

    int modlist[4] = { 0, 2, 3, 1 }; // mapping Simon along the beam

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg ) std::cout << "PLANE " << plane.ID() << ": ";

      int mod = modlist[iplane];

      if( mod > 3 ) {
	cout << "illegal plane " << iplane << endl;
	continue;
      }

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix) {


	if( ldbg ) 
	  std::cout << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << " ";

	int xm = plane.GetX(ipix); // global column 0..415
	int ym = plane.GetY(ipix); // global row 0..159
	int adc = plane.GetPixel(ipix); // ADC 0..255

	// leave space for big pixels:

	int roc = xm / 52; // 0..7
	int x = 1 + xm + 2*roc; // 1..52 per ROC
	int y = ym;
	if( ym > 79 ) y += 2;

	// fill pixel block for clustering:
	pb[npx].col = x;
	pb[npx].row = y;
	pb[npx].adc = adc;
	pb[npx].cal = adc;
	++npx;
	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	if( xm%52 == 0 ) {
	  pb[npx].col = x-1; // double
	  pb[npx].row = y;
	  pb[npx].adc = adc;
	  pb[npx].cal = adc;
	  ++npx;
	}

	if( xm%52 == 51 ) {
	  pb[npx].col = x+1; // double
	  pb[npx].row = y;
	  pb[npx].adc = adc;
	  pb[npx].cal = adc;
	  ++npx;
	}

	if( ym == 79 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 80;
	  pb[npx].adc = adc;
	  pb[npx].cal = adc;
	  ++npx;
	}

	if( ym == 80 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 81;
	  pb[npx].adc = adc;
	  pb[npx].cal = adc;
	  ++npx;
	}

      } // pix
      
      if( ldbg ) std::cout << std::endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[mod] = getClus();

      if( ldbg ) cout << "clusters mod " << mod << " " << cl[mod].size() << endl;

    } // planes = mod

    if( cl[0].size() < 1 ) continue; // skip event
    if( cl[1].size() < 1 ) continue;
    if( cl[2].size() < 1 ) continue;
    if( cl[3].size() < 1 ) continue;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // tracks in x:

    int n4x = 0;
    vector <track> xtrks;

    for( vector<cluster>::iterator cA = cl[0].begin(); cA != cl[0].end(); ++cA ) {

      double xA = cA->col*0.15 - alignx[0];
      double yA = cA->row*0.10 - aligny[0];
      double xmid = xA - 32.4;
      double ymid = yA -  8.1;
      double xrot = xmid     - ymid*fx[0] - tx[0]*xmid;
      double yrot = xmid*fy[0] + ymid     - ty[0]*ymid;
      xA = xrot + 32.4;
      yA = yrot +  8.1;

      for( vector<cluster>::iterator cD = cl[3].begin(); cD != cl[3].end(); ++cD ) {

	double xD = cD->col*0.15 - alignx[3];
	double yD = cD->row*0.10 - aligny[3];
	double xmid = xD - 32.4;
	double ymid = yD -  8.1;
	double xrot = xmid     - ymid*fx[3] - tx[3]*xmid;
	double yrot = xmid*fy[3] + ymid     - ty[3]*ymid;
	xD = xrot + 32.4;
	yD = yrot +  8.1;

	double dx2 = xD - xA;

	//if( abs( dx2 ) > 0.20 ) continue; // angle cut

	double xavg2 = (xA + 2*xD)/3; // equidistant

	for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	  double xC = cC->col*0.15 - alignx[2];
	  double yC = cC->row*0.10 - aligny[2];
	  double xmid = xC - 32.4;
	  double ymid = yC -  8.1;
	  double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
	  double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
	  xC = xrot + 32.4;
	  yC = yrot +  8.1;

	  double dx3 = xC - xavg2;

	  if( abs( dx3 ) > 0.15 + 0.02*abs(dx2) ) continue;

	  double xavg3 = 0.5*(xA + xC); // equidistant

	  cout << "\ttri  dx " << dx3 << endl;

	  // quad with B:

	  for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	    double xB = cB->col*0.15; // REF

	    double dx4 = xB - xavg3;

	    if( abs( dx4 ) > 0.15 + 0.02*abs(dx2) ) continue;

	    cout << "\tquad dx " << dx4 << endl;
	    ++n4x;
	    track t;
	    t.xA = xA;
	    t.yA = yA;
	    t.xD = xD;
	    t.yD = yD;
	    xtrks.push_back(t);

	  } // cl B

	} // cl C

      } // cl D

    } // cl A

    cout << "x tracks " << n4x << endl;
    //if( n4x < 5 ) continue; // skip event

    //tracks in y:

    int n4y = 0;
    vector <track> ytrks;

    for( vector<cluster>::iterator cA = cl[0].begin(); cA != cl[0].end(); ++cA ) {

      double xA = cA->col*0.15 - alignx[0];
      double yA = cA->row*0.10 - aligny[0];
      double xmid = xA - 32.4;
      double ymid = yA -  8.1;
      double xrot = xmid     - ymid*fx[0] - tx[0]*xmid;
      double yrot = xmid*fy[0] + ymid     - ty[0]*ymid;
      xA = xrot + 32.4;
      yA = yrot +  8.1;

      for( vector<cluster>::iterator cD = cl[3].begin(); cD != cl[3].end(); ++cD ) {

	double xD = cD->col*0.15 - alignx[3];
	double yD = cD->row*0.10 - aligny[3];
	double xmid = xD - 32.4;
	double ymid = yD -  8.1;
	double xrot = xmid     - ymid*fx[3] - tx[3]*xmid;
	double yrot = xmid*fy[3] + ymid     - ty[3]*ymid;
	xD = xrot + 32.4;
	yD = yrot +  8.1;

	double dy2 = yD - yA;

	double yavg2 = (yA + 2*yD)/3; // C is closer to D: more weight for D

	for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	  double xC = cC->col*0.15 - alignx[2];
	  double yC = cC->row*0.10 - aligny[2];
	  double xmid = xC - 32.4;
	  double ymid = yC -  8.1;
	  double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
	  double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
	  xC = xrot + 32.4;
	  yC = yrot +  8.1;

	  double dy3 = yC - yavg2;

	  if( abs( dy3 ) > 0.15 + 0.02*abs(dy2) ) continue;

	  double yavg3 = 0.5*(yA + yC);

	  // quad with B:

	  for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	    double yB = cB->row*0.10;

	    double dy4 = yB - yavg3;

	    if( abs( dy4 ) > 0.15 + 0.02*abs(dy2) ) continue;

	    ++n4y;
	    track t;
	    t.xA = xA;
	    t.yA = yA;
	    t.xD = xD;
	    t.yD = yD;
	    ytrks.push_back(t);

	  } // cl B

	} // cl C

      } // cl D

    } // cl A

    cout << "y tracks " << n4y << endl;
    //if( n4y < 5 ) continue; // skip event

    //tracks in x and y:

    int n4 = 0;
    vector <track> trks;

    for( vector<cluster>::iterator cA = cl[0].begin(); cA != cl[0].end(); ++cA ) {

      double xA = cA->col*0.15 - alignx[0];
      double yA = cA->row*0.10 - aligny[0];
      double xmid = xA - 32.4;
      double ymid = yA -  8.1;
      double xrot = xmid     - ymid*fx[0] - tx[0]*xmid;
      double yrot = xmid*fy[0] + ymid     - ty[0]*ymid;
      xA = xrot + 32.4;
      yA = yrot +  8.1;

      for( vector<cluster>::iterator cD = cl[3].begin(); cD != cl[3].end(); ++cD ) {

	double xD = cD->col*0.15 - alignx[3];
	double yD = cD->row*0.10 - aligny[3];
	double xmid = xD - 32.4;
	double ymid = yD -  8.1;
	double xrot = xmid     - ymid*fx[3] - tx[3]*xmid;
	double yrot = xmid*fy[3] + ymid     - ty[3]*ymid;
	xD = xrot + 32.4;
	yD = yrot +  8.1;

	double dx2 = xD - xA;
	double dy2 = yD - yA;

	//if( abs( dx2 ) > 0.20 ) continue; // angle cut
	//if( abs( dy2 ) > 0.20 ) continue;

	double xavg2 = (xA + 2*xD)/3; // equidistant
	double yavg2 = (yA + 2*yD)/3; // C is closer to D: more weight for D

	for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	  double xC = cC->col*0.15 - alignx[2];
	  double yC = cC->row*0.10 - aligny[2];
	  double xmid = xC - 32.4;
	  double ymid = yC -  8.1;
	  double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
	  double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
	  xC = xrot + 32.4;
	  yC = yrot +  8.1;

	  double dx3 = xC - xavg2;
	  double dy3 = yC - yavg2;

	  if( abs( dx3 ) > 0.15 + 0.02*abs(dx2) ) continue;
	  if( abs( dy3 ) > 0.15 + 0.02*abs(dy2) ) continue;

	  double xavg3 = 0.5*(xA + xC); // equidistant
	  double yavg3 = 0.5*(yA + yC);

	  cout << "\t\ttri  dx " << dx3 << "\tdy " << dy3 << endl;

	  // quad with B:

	  for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	    double xB = cB->col*0.15; // REF
	    double yB = cB->row*0.10;

	    double dx4 = xB - xavg3;
	    double dy4 = yB - yavg3;

	    if( abs( dx4 ) > 0.15 + 0.02*abs(dx2) ) continue;
	    if( abs( dy4 ) > 0.15 + 0.02*abs(dy2) ) continue;

	    cout << "\t\tquad dx " << dx4 << "\tdy " << dy4 << endl;
	    ++n4;
	    track t;
	    t.xA = xA;
	    t.yA = yA;
	    t.xD = xD;
	    t.yD = yD;
	    trks.push_back(t);

	  } // cl B

	} // cl C

      } // cl D

    } // cl A

    cout << "xy tracks " << n4 << endl;
    if( n4 < 3 ) continue; // skip event

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // book event histos:

    ++kev;

    cout << "display " << kev << endl;

    TH2D xview( Form( "ev%ix", kev ),
		Form( "display %i trigger %i x;x [mm];plane;ADC", kev, (int) event_nr ),
		8*54, 0, 8*8.1, 20, 0.5, 4.5 );

    TH2D yview( Form( "ev%iy", kev ),
		Form( "display %i trigger %i y;y [mm];plane;ADC", kev, (int) event_nr ),
		2*81, 0, 2*8.1, 20, 0.5, 4.5 );

    for( vector<cluster>::iterator c = cl[0].begin(); c != cl[0].end(); ++c ) {

      double x = c->col*0.15 - alignx[0];
      double y = c->row*0.10 - aligny[0];
      double xmid = x - 32.4;
      double ymid = y -  8.1;
      double xrot = xmid     - ymid*fx[0] - tx[0]*xmid;
      double yrot = xmid*fy[0] + ymid     - ty[0]*ymid;
      x = xrot + 32.4;
      y = yrot +  8.1;

      xview.Fill( x, 1.0, c->sumA );
      yview.Fill( y, 1.0, c->sumA );

    }

    for( vector<cluster>::iterator c = cl[1].begin(); c != cl[1].end(); ++c ) {

      double x = c->col*0.15; // REF
      double y = c->row*0.10;

      xview.Fill( x, 2.0, c->sumA );
      yview.Fill( y, 2.0, c->sumA );

    }

    for( vector<cluster>::iterator c = cl[2].begin(); c != cl[2].end(); ++c ) {

      double x = c->col*0.15 - alignx[2];
      double y = c->row*0.10 - aligny[2];
      double xmid = x - 32.4;
      double ymid = y -  8.1;
      double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
      double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
      x = xrot + 32.4;
      y = yrot +  8.1;

      xview.Fill( x, 3.0, c->sumA );
      yview.Fill( y, 3.0, c->sumA );

    }

    for( vector<cluster>::iterator c = cl[3].begin(); c != cl[3].end(); ++c ) {

      double x = c->col*0.15 - alignx[3];
      double y = c->row*0.10 - aligny[3];
      double xmid = x - 32.4;
      double ymid = y -  8.1;
      double xrot = xmid     - ymid*fx[3] - tx[3]*xmid;
      double yrot = xmid*fy[3] + ymid     - ty[3]*ymid;
      x = xrot + 32.4;
      y = yrot +  8.1;

      xview.Fill( x, 4.0, c->sumA );
      yview.Fill( y, 4.0, c->sumA );

    }

    xview.Write();
    xview.GetXaxis()->SetNdivisions(-408);
    xview.GetYaxis()->SetNdivisions(-4);
    xview.GetYaxis()->SetRangeUser( xview.GetYaxis()->GetXmin(),
				    xview.GetYaxis()->GetXmax() );
    //xview.Draw("colz");
    xview.Draw("");

    vector <TLine> xlines;
    for( size_t ii = 0; ii < xtrks.size(); ++ii ) {
      TLine l( xtrks[ii].xA, 1, xtrks[ii].xD, 4 );
      l.SetLineColor(1);
      xlines.push_back(l);
    }
    for( size_t ii = 0; ii < trks.size(); ++ii ) {
      TLine l( trks[ii].xA, 1, trks[ii].xD, 4 );
      l.SetLineColor(2);
      xlines.push_back(l);
    }
    for( size_t ii = 0; ii < xlines.size(); ++ii ) {
      xlines[ii].Draw("same");
    }

    c1->Update();

    string input;
    cout << "hit enter for y...";
    getline( cin, input, '\n' );

    yview.Write();
    yview.GetXaxis()->SetNdivisions(-402);
    yview.GetYaxis()->SetNdivisions(-4);
    yview.GetYaxis()->SetRangeUser( yview.GetYaxis()->GetXmin(),
				    yview.GetYaxis()->GetXmax() );
    //yview.Draw("colz");
    yview.Draw("");

    vector <TLine> ylines;
    for( size_t ii = 0; ii < ytrks.size(); ++ii ) {
      TLine l( ytrks[ii].yA, 1, ytrks[ii].yD, 4 );
      l.SetLineColor(1);
      ylines.push_back(l);
    }
    for( size_t ii = 0; ii < trks.size(); ++ii ) {
      TLine l( trks[ii].yA, 1, trks[ii].yD, 4 );
      l.SetLineColor(4);
      ylines.push_back(l);
    }
    for( size_t ii = 0; ii < ylines.size(); ++ii ) {
      ylines[ii].Draw("same");
    }

    c1->Update();

    cout << "hit enter for next...";
    getline( cin, input, '\n' );
      
  } while( reader.NextEvent() && kev < lev );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  histoFile->Write();
  histoFile->Close();
  cout << "write " << histoFile->GetName() << endl;
  delete histoFile;

  delete myMF;

  return 0;
}
