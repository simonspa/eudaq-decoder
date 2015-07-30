#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>

using namespace eudaq;
using namespace std;

struct pixel {
  int col;
  int row;
  int adc;
  double cal;
  bool big;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int sumA; // DP
  double charge;
  double col,row;
  bool big;
  int ncol, nrow;
};

// globals:

pixel pb[66560]; // global declaration: vector of pixels with hit
int fNHit; // global

//------------------------------------------------------------------------------
// inverse decorrelated Weibull PH -> large Vcal DAC
double PHtoVcal( double ph, double a0, double a1, double a2, double a3, double a4, double a5, int mod )
{
  // modph2ps decorrelated: ph = a4 - a3*exp(-t^a2), t = a0 + q/a1

  //return ph; // test !!

  double Ared = ph - a4; // a4 is asymptotic maximum

  if( Ared >= 0 ) {
    Ared = -0.1; // avoid overflow
  }

  // large Vcal = ( (-ln(-(A-a4)/a3))^1/a2 - a0 )*a1

  if( a3 < 1E-9 ) {
    cout << "PHtoVcal zero a3  " << a3 << endl;
    return ph;
  }
  else if( -Ared > a3 ) {
    cout << "PHtoVcal small a3  " << a3 << "  " << -Ared << " mod " << mod << endl;
    return ph;
  }

  double vc =
    a1 * ( pow( -log( -Ared / a3 ), 1/a2 ) - a0 );

  if( vc > 999 )
    cout << "overflow " << vc << ", Ared " << Ared << ", a3 " << a3 << endl;

  if( vc != vc ) {

    cout << "PHtoVcal NaN at "
	 << "  " << a0
	 << "  " << a1
	 << "  " << a2
	 << "  " << a3
	 << "  " << a4
	 << "  " << a5
	 << endl;

    return ph;

  }

  return vc * a5; // small Vcal
  //return vc; // large Vcal
}

//------------------------------------------------------------------------------
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

    // added all I could. determine position and append it to the list of clusters:

    c.sumA = 0;
    c.charge = 0;
    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    c.big = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      c.sumA += p->adc; // Aout
      double Qpix = p->cal; // calibrated [Vcal]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      if( p->big ) c.big = 1;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

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

    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

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

  // alignments:

  double alignx[4];
  double aligny[4];
  double fx[4];
  double fy[4];
  double tx[4];
  double ty[4];

  if( run >= 70 ) { // Sun 26.7.2015

    alignx[0] = -1.162; // [mm] same sign as dx
    aligny[0] = -0.448; // [mm] same sign as dy
    fx[0] = -0.004; // [rad] rot, same     sign dxvsy
    fy[0] = -0.007; // [rad] rot, opposite sign dyvsx
    tx[0] = 0.0000; // from dxvsx
    ty[0] = 0.0012; // from dyvsy

    alignx[2] = -0.413; // [mm] same sign as dx
    aligny[2] =  0.292; // [mm] same sign as dy
    fx[2] = -0.0055; // [rad] rot, same     sign dxvsy
    fy[2] =  0.0048; // [rad] rot, opposite sign dyvsx
    tx[2] = 0.0000; // from dxvsx
    ty[2] =-0.0050; // from dyvsy, same sign

    alignx[3] =  0.288; // [mm] same sign as dx
    aligny[3] =  0.588; // [mm] same sign as dy
    fx[3] = -0.0015; // [rad] rot, same     sign dxvsy
    fy[3] =  0.0056; // [rad] rot, opposite sign dyvsx
    tx[3] = 0.0004; // from dxvsx
    ty[3] =-0.0024; // from dyvsy, same sign

  }

  string gainFileName[4];
  double ke[4];

  if( run > 70 ) {

    gainFileName[0] = "D4028/D4028-tb24-gaincal.dat"; // Jul 2015
    ke[0] = 0.046; // small Vcal -> ke

    gainFileName[1] = "D4035/D4035-tb24-gaincal.dat"; // Jul 2015
    ke[1] = 0.050; // small Vcal -> ke

    gainFileName[2] = "D4030/D4030-tb24-gaincal.dat"; // Jul 2015
    ke[2] = 0.049; // small Vcal -> ke

    gainFileName[3] = "D4023/D4023-tb24-gaincal.dat"; // Jul 2015
    ke[3] = 0.050; // small Vcal -> ke

  }

  double ****p0 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p0[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p0[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p0[i][j][k] = new double[80];

  double ****p1 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p1[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p1[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p1[i][j][k] = new double[80];

  double ****p2 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p2[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p2[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p2[i][j][k] = new double[80];

  double ****p3 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p3[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p3[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p3[i][j][k] = new double[80];

  double ****p4 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p4[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p4[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p4[i][j][k] = new double[80];

  double ****p5 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p5[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p5[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p5[i][j][k] = new double[80];

  bool haveGain[4] = {0};

  for( int mod = 0; mod < 4; ++mod ) {

    if( gainFileName[mod].length(  ) > 0 ) {

      ifstream gainFile( gainFileName[mod].c_str() );

      if( gainFile ) {

	haveGain[mod] = 1;
	cout << "gain " << mod << ": " << gainFileName[mod] << endl;

	int roc;
	int col;
	int row;
	double a0, a1, a2, a3, a4, a5;

	while( gainFile >> roc ) {
	  gainFile >> col;
	  gainFile >> row;
	  gainFile >> a0;
	  gainFile >> a1;
	  gainFile >> a2;
	  gainFile >> a3;
	  gainFile >> a4;
	  gainFile >> a5;
	  p0[mod][roc][col][row] = a0;
	  p1[mod][roc][col][row] = a1;
	  p2[mod][roc][col][row] = a2;
	  p3[mod][roc][col][row] = a3;
	  p4[mod][roc][col][row] = a4;
	  p5[mod][roc][col][row] = a5;
	}

      } // gainFile open

    } // gainFileName

  } // mod


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile* histoFile = new TFile( "main.root", "RECREATE" );

  // book histos:
  /*
  TH1D hcolA( "colA", "A col;col;A pixels", 416, -0.5, 415.5 );
  TH1D hrowA( "rowA", "A row;row;A pixels", 160, -0.5, 159.5 );
  TH1D hpxqA( "pxqA", "A pixel charge;pixel q [ke];A pixels",
	      100, 0, 25 );
  TH2D * hmapA;
  hmapA = new TH2D( "pxmapA", "A pixel map;column;row;A pixels",
		    416, -0.5, 415.5, 160, -0.5, 159.5 );
  TH1D hnpxA( "npxA", "A pixel per event;pixels;A events",
	      51, -0.5, 50.5 );
  TH1D hsizA( "clszA", "A cluster size;pixels/cluster;A clusters",
	      51, -0.5, 50.5 );
  TH1D hclqA( "clqA", "A cluster charge;cluster charge [ke];A clusters",
	      100, 0, 100 );
  TH1D hncolA( "ncolA", "A cluster size;columns/cluster;A clusters",
	       21, -0.5, 20.5 );
  TH1D hnrowA( "nrowA", "A cluster size;rows/cluster;A clusters",
	       21, -0.5, 20.5 );
  */

  TH1D hncl[4];
  for( int mod = 0; mod < 4; ++mod ) {

    hncl[mod] = TH1D( Form( "ncl%i", mod ),
		      Form( "plane %i cluster per event;cluster;plane %i events", mod, mod ),
		      51, -0.5, 50.5 );

  } // module planes

  TH2D hxxAB( "xxAB", "A vs B;col B;col A;clusters",
	      432, 0, 64.8, 432, 0, 64.8 );
  TH2D hyyAB( "yyAB", "A vs B;row B;row A;clusters",
	      162, 0, 16.2, 162, 0, 16.2 );

  TH1D hdxAB( "dxAB", "Ax-Bx;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyAB( "dyAB", "Ay-By;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcAB( "dxcAB", "Ax-Bx;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycAB( "dycAB", "Ay-By;y-y [mm];cluster pairs", 200, -1, 1 );

  TProfile dxvsxAB( "dxvsxAB", "A-B dx vs x;x [mm];A-B <dx>",
		    216, 0, 64.8, -1, 1 );
  TProfile dxvsyAB( "dxvsyAB", "A-B dx vs y;y [mm];A-B <dx>",
		    81, 0, 16.2, -1, 1 );
  TProfile dyvsxAB( "dyvsxAB", "A-B dy vs x;x [mm];A-B <dy>",
		    216, 0, 64.8, -1, 1 );
  TProfile dyvsyAB( "dyvsyAB", "A-B dy vs y;y [mm];A-B <dy>",
		    81, 0, 16.2, -1, 1 );

  TH2D hxxCB( "xxCB", "C vs B;col B;col C;clusters",
	      432, 0, 64.8, 432, 0, 64.8 );
  TH2D hyyCB( "yyCB", "C vs B;row B;row C;clusters",
	      162, 0, 16.2, 162, 0, 16.2 );

  TH1D hdxCB( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyCB( "dyCB", "Cy-By;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcCB( "dxcCB", "Cx-Bx;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycCB( "dycCB", "Cy-By;y-y [mm];cluster pairs", 200, -1, 1 );

  TProfile dxvsxCB( "dxvsxCB", "C-B dx vs x;x [mm];C-B <dx>",
		    216, 0, 64.8, -1, 1 );
  TProfile dxvsyCB( "dxvsyCB", "C-B dx vs y;y [mm];C-B <dx>",
		    81, 0, 16.2, -1, 1 );
  TProfile dyvsxCB( "dyvsxCB", "C-B dy vs x;x [mm];C-B <dy>",
		    216, 0, 64.8, -1, 1 );
  TProfile dyvsyCB( "dyvsyCB", "C-B dy vs y;y [mm];C-B <dy>",
		    81, 0, 16.2, -1, 1 );

  TH2D hxxDC( "xxDC", "D vs C;col C;col D;clusters",
	      432, 0, 64.8, 432, 0, 64.8 );
  TH2D hyyDC( "yyDC", "D vs C;row C;row D;clusters",
	      162, 0, 16.2, 162, 0, 16.2 );

  TH1D hdxDC( "dxDC", "Dx-Cx;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyDC( "dyDC", "Dy-Cy;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcDC( "dxcDC", "Dx-Cx;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycDC( "dycDC", "Dy-Cy;y-y [mm];cluster pairs", 200, -1, 1 );

  TProfile dxvsxDC( "dxvsxDC", "D-C dx vs x;x [mm];D-C <dx>",
		    216, 0, 64.8, -1, 1 );
  TProfile dxvsyDC( "dxvsyDC", "D-C dx vs y;y [mm];D-C <dx>",
		    81, 0, 16.2, -1, 1 );
  TProfile dyvsxDC( "dyvsxDC", "D-C dy vs x;x [mm];D-C <dy>",
		    216, 0, 64.8, -1, 1 );
  TProfile dyvsyDC( "dyvsyDC", "D-C dy vs y;y [mm];D-C <dy>",
		    81, 0, 16.2, -1, 1 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  size_t event_nr = 0;
  
  do {
    // Get next event:
    DetectorEvent evt = reader.GetDetectorEvent();

    if (evt.IsBORE()) { eudaq::PluginManager::Initialize(evt); }

    bool ldbg = 0;

    if (event_nr%1000==0)
      ldbg = 1;

    if( ldbg ) std::cout<<"Processing event "<< event_nr<<std::endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    int xm = 0;
    int ym = 0;
    int adc = 0;
    double cal = 0;

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

	xm = plane.GetX(ipix); // global column 0..415
	ym = plane.GetY(ipix); // global row 0..159
	adc = plane.GetPixel(ipix); // ADC 0..255

	// to do: gaincal

	cal = adc;

	// leave space for big pixels:

	int roc = xm / 52; // 0..7
	int col = xm % 52; // 0..51
	int row = ym;
	int x = 1 + xm + 2*roc; // 1..52 per ROC with big pix
	int y = ym;
	if( ym > 79 ) y += 2;

	// flip for upper ROCs into local addresses:

	if( ym > 79 ) {
	  roc = 15 - roc; // 15..8
	  col = 51 - col; // 51..0
	  row = 159 - ym; // 79..0
	}

	if( xm < 0 || xm > 415 || ym < 0 || ym > 159 || adc < 0 || adc > 255 )
	  cout << "invalid pixel at event " << event_nr << endl;
	else if( haveGain[mod] ) {
	  double a0 = p0[mod][roc][col][row];
	  double a1 = p1[mod][roc][col][row];
	  double a2 = p2[mod][roc][col][row];
	  double a3 = p3[mod][roc][col][row];
	  double a4 = p4[mod][roc][col][row];
	  double a5 = p5[mod][roc][col][row];
	  cal = PHtoVcal( adc, a0, a1, a2, a3, a4, a5, mod ); // [Vcal]
	}

	pb[npx].col = x;
	pb[npx].row = y;
	pb[npx].adc = adc;
	pb[npx].cal = cal;
	pb[npx].big = 0;
	++npx;

	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	col = xm % 52; // 0..51

	if( col == 0 ) {
	  pb[npx].col = x-1; // double
	  pb[npx].row = y;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

	if( col == 51 ) {
	  pb[npx].col = x+1; // double
	  pb[npx].row = y;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

	if( ym == 79 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 80;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

	if( ym == 80 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 81;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

      } // pix

      if( ldbg ) std::cout << std::endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[mod] = getClus();

      if( ldbg ) cout << "A clusters " << cl[mod].size() << endl;

      hncl[mod].Fill( cl[mod].size() );

    } // planes = mod

    event_nr++;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-B cluster correlations:

    for( vector<cluster>::iterator cA = cl[0].begin(); cA != cl[0].end(); ++cA ) {

      double xA = cA->col*0.15 - alignx[0];
      double yA = cA->row*0.10 - aligny[0];
      double xmid = xA - 32.4;
      double ymid = yA -  8.1;
      double xrot = xmid     - ymid*fx[0] - tx[0]*xmid;
      double yrot = xmid*fy[0] + ymid     - ty[0]*ymid;
      xA = xrot + 32.4;
      yA = yrot +  8.1;

      for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	double xB = cB->col*0.15; // REF
	double yB = cB->row*0.10;

	hxxAB.Fill( xB, xA );
	hyyAB.Fill( yB, yA );

	double dx = xA - xB;
	double dy = yA - yB;
	hdxAB.Fill( dx );
	hdyAB.Fill( dy );
	if( abs( dy ) < 0.15 && cA->big == 0 && cB->big == 0 ) {
	  hdxcAB.Fill( dx );
	  dxvsxAB.Fill( xB, dx );
	  dxvsyAB.Fill( yB, dx );
	}
	if( abs( dx ) < 0.20 && cA->big == 0 && cB->big == 0 ) {
	  hdycAB.Fill( dy );
	  dyvsxAB.Fill( xA, dy );
	  dyvsyAB.Fill( yA, dy );
	}

      } // clusters

    } // cl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // C-B cluster correlations:

    for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

      double xC = cC->col*0.15 - alignx[2];
      double yC = cC->row*0.10 - aligny[2];
      double xmid = xC - 32.4;
      double ymid = yC -  8.1;
      double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
      double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
      xC = xrot + 32.4;
      yC = yrot +  8.1;

      for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	double xB = cB->col*0.15; // REF
	double yB = cB->row*0.10;

	hxxCB.Fill( xB, xC );
	hyyCB.Fill( yB, yC );

	double dx = xC - xB;
	double dy = yC - yB;
	hdxCB.Fill( dx );
	hdyCB.Fill( dy );
	if( abs( dy ) < 0.15 && cB->big == 0 && cC->big == 0 ) {
	  hdxcCB.Fill( dx );
	  dxvsxCB.Fill( xC, dx );
	  dxvsyCB.Fill( yC, dx );
	}
	if( abs( dx ) < 0.20 && cB->big == 0 && cC->big == 0 ) {
	  hdycCB.Fill( dy );
	  dyvsxCB.Fill( xC, dy );
	  dyvsyCB.Fill( yC, dy );
	}

      } // clusters

    } // cl

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // D-C correlations:

    for( vector<cluster>::iterator cD = cl[3].begin(); cD != cl[3].end(); ++cD ) {

      double xD = cD->col*0.15 - alignx[3];
      double yD = cD->row*0.10 - aligny[3];
      double xmid = xD - 32.4;
      double ymid = yD -  8.1;
      double xrot = xmid     - ymid*fx[3] - tx[3]*xmid;
      double yrot = xmid*fy[3] + ymid     - ty[3]*ymid;
      xD = xrot + 32.4;
      yD = yrot +  8.1;

      for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	double xC = cC->col*0.15 - alignx[2];
	double yC = cC->row*0.10 - aligny[2];
	double xmid = xC - 32.4;
	double ymid = yC -  8.1;
	double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
	double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
	xC = xrot + 32.4;
	yC = yrot +  8.1;

	hxxDC.Fill( xC, xD );
	hyyDC.Fill( yC, yD );

	double dx = xD - xC;
	double dy = yD - yC;
	hdxDC.Fill( dx );
	hdyDC.Fill( dy );
	if( abs( dy ) < 0.15 && cD->big == 0 && cC->big == 0 ) {
	  hdxcDC.Fill( dx );
	  dxvsxDC.Fill( xD, dx );
	  dxvsyDC.Fill( yD, dx );
	}
	if( abs( dx ) < 0.20 && cD->big == 0 && cC->big == 0 ) {
	  hdycDC.Fill( dy );
	  dyvsxDC.Fill( xD, dy );
	  dyvsyDC.Fill( yD, dy );
	}

      } // cl

    } // cl

    // tracking


  } while (reader.NextEvent());

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

  return 0;
}
