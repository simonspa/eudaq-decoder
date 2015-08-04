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

  // Landau peak cuts: Mon 27.7.2015

  double qL = 20; // [ke]
  double qR = 36; // [ke]

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
 
  TH1D hcol[4];
  TH1D hrow[4];
  TH1D hpxq[4];
  TH2D * hmap[4];
  TH1D hnpx[4];
  TH1D hsiz[4];
  TH1D hclq[4];
  TH1D hncol[4];
  TH1D hnrow[4];

  TH1D hncl[4];
  for( int mod = 0; mod < 4; ++mod ) {
    char modtos;
    switch(mod){
    case 0: modtos = 'A'; break;
    case 1: modtos =  'B'; break;
    case 2: modtos = 'C'; break;  
    case 3: modtos = 'D'; break;  
    }
 
    hncl[mod] = TH1D( Form( "ncl%c", modtos ),
		      Form( "plane %c cluster per event;cluster;plane %c events", modtos, modtos ),
		      51, -0.5, 50.5 );
    hcol[mod] = TH1D( Form("col%c", modtos),
		      Form("%c col;col;%c pixels", modtos, modtos), 
		      416, -0.5, 415.5 );
    hrow[mod] = TH1D( Form("row%c",modtos),
		      Form("%c row;row;%c pixels",modtos,modtos),
		      160, -0.5, 159.5 );
    hpxq[mod] = TH1D( Form("pxq%c",modtos),
		      Form("%c pixel charge;pixel q [ke];%c pixels",modtos,modtos),
		      100, 0, 25 );
    hmap[mod] = new  TH2D( Form("pxmap%c",modtos),
		      Form("%c pixel map;column;row;%c pixels",modtos,modtos),
		      416, -0.5, 415.5, 160, -0.5, 159.5 );
    hnpx[mod] = TH1D( Form("npx%c",modtos),
		      Form("%c pixel per event;pixels;%c events",modtos,modtos),
		      51, -0.5, 50.5 );
    hsiz[mod] = TH1D( Form("clsz%c",modtos),
		      Form("%c cluster size;pixels/cluster;%c clusters",modtos,modtos),
		      51, -0.5, 50.5 );
    hclq[mod] = TH1D( Form("clq%c",modtos),
		      Form("%c cluster charge;cluster charge [ke];%c clusters",modtos,modtos),
		      100, 0, 100 );
    hncol[mod]= TH1D( Form("ncol%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hnrow[mod]= TH1D( Form("nrow%c",modtos),
		      Form("%c cluster size;rows/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );

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
  TH2D hxxDA( "xxDA", "D vs A;col A;col D;clusters",
	      432, 0, 64.8, 432, 0, 64.8 );
  TH2D hyyDA( "yyDA", "D vs A;row A;row D;clusters",
	      162, 0, 16.2, 162, 0, 16.2 );

  TH1D hdxDA( "dxDA", "Dx-Ax;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyDA( "dyDA", "Dy-Ay;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcDA( "dxcDA", "Dx-Ax;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycDA( "dycDA", "Dy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );

  TProfile dxvsxDA( "dxvsxDA", "D-A dx vs x;x [mm];D-A <dx>",
		    216, 0, 64.8, -1, 1 );
  TProfile dxvsyDA( "dxvsyDA", "D-A dx vs y;y [mm];D-A <dx>",
		    81, 0, 16.2, -1, 1 );
  TProfile dyvsxDA( "dyvsxDA", "D-A dy vs x;x [mm];D-A <dy>",
		    216, 0, 64.8, -1, 1 );
  TProfile dyvsyDA( "dyvsyDA", "D-A dy vs y;y [mm];D-A <dy>",
		    81, 0, 16.2, -1, 1 );

  // triplets ADC:

  TH1D hdxADC( "dxADC", "ADC dx;x-x [mm];ADCplets", 200, -1, 1 );
  TH1D hdyADC( "dyADC", "ADC dy;y-y [mm];ADCplets", 200, -1, 1 );
  TH1D hdxcADC( "dxcADC", "ADC dx;x-x [um];ADCplets", 200, -200, 200 );
  TH1D hdycADC( "dycADC", "ADC dy;y-y [um];ADCplets", 200, -200, 200 );
  TH1D hdxciADC( "dxciADC", "ADC dx;x-x [um];isolated ADCplets",
		 200, -200, 200 );
  TH1D hdyciADC( "dyciADC", "ADC dy;y-y [um];isolated ADCplets",
		 200, -200, 200 );

  TProfile dxvsxADC( "dxvsxADC", "ADCplet dx vs x;x [mm];ADCplet <dx>",
		     216, 0, 64.8, -1, 1 );
  TProfile dxvsyADC( "dxvsyADC", "ADCplet dx vs y;y [mm];ADCplet <dx>",
		     81, 0, 16.2, -1, 1 );
  TProfile dyvsxADC( "dyvsxADC", "ADCplet dy vs x;x [mm];ADCplet <dy>",
		     216, 0, 64.8, -1, 1 );
  TProfile dyvsyADC( "dyvsyADC", "ADCplet dy vs y;y [mm];ADCplet <dy>",
		     81, 0, 16.2, -1, 1 );
  TH1D hxADC( "xADC", "ADCplets;col;ADCplets", 216, 0, 64.8 );
  TH1D hyADC( "yADC", "ADCplets;row;ADCplets",  81, 0, 16.2 );
  TH2D * hmapADC;
  hmapADC = new TH2D( "mapADC", "ADCplet map;ADCplet col;ADCplet row;ADCplets",
		      8*54, 0, 8*8.1, 2*81, 0, 2*8.1 );
  TH1D htxADC( "txADC", "ADCplet angle x;ADCplet angle x;ADCplets", 100, -1, 1 );
  TH1D htyADC( "tyADC", "ADCplet angle y;ADCplet angle y;ADCplets", 100, -1, 1 );

  // triplets ABC:

  TH1D hdxACB( "dxACB", "ACB dx;x-x [mm];ACBs", 200, -1, 1 );
  TH1D hdyACB( "dyACB", "ACB dy;y-y [mm];ACBs", 200, -1, 1 );
  TH1D hdxcACB( "dxcACB", "ACB dx;x-x [um];ACBs", 200, -200, 200 );
  TH1D hdycACB( "dycACB", "ACB dy;y-y [um];ACBs", 200, -200, 200 );
  TH1D hdxciACB( "dxciACB", "ACB dx;x-x [um];isolated ACBs",
		 200, -200, 200 );
  TH1D hdyciACB( "dyciACB", "ACB dy;y-y [um];isolated ACBs",
		 200, -200, 200 );
  TH1D hdycfACB( "dycfACB", "ACB dy;y-y [um];inner ACBs",
		 200, -200, 200 );
  TH1D hdycfqACB( "dycfqACB", "ACB dy;y-y [um];Landau peak inner ACBs",
		  200, -200, 200 );

  TProfile dxvsxACB( "dxvsxACB",
		     "ACB dx vs x;x [mm];ACB <dx>",
		     216, 0, 64.8, -1, 1 );
  TProfile dxvsyACB( "dxvsyACB",
		     "ACB dx vs y;y [mm];ACB <dx>",
		     81, 0, 16.2, -1, 1 );
  TProfile dyvsxACB( "dyvsxACB",
		     "ACB dy vs x;x [mm];ACB <dy>",
		     216, 0, 64.8, -1, 1 );
  TProfile dyvsyACB( "dyvsyACB",
		     "ACB dy vs y;y [mm];ACB <dy>",
		     81, 0, 16.2, -1, 1 );
  TProfile madyvsyACB( "madyvsyACB",
		       "ACB mady vs y;y [mm];ACB MAD(y) [um]",
		       81, 0, 16.2, 0, 100 );
  TProfile madyvsymACB( "madyvsymACB",
			"ACB mady vs ymod;y mod 200 [um];ACB MAD(y) [um]",
			40, 0, 200, 0, 100 );
  TH2D * hmapACB;
  hmapACB = new TH2D( "mapACB",
		      "ACBplet map;ACBplet col;ACBplet row;ACBplets",
		      8*54, 0, 8*8.1, 2*81, 0, 2*8.1 );

  // Bvs AC:

  TH1D hsizB4( "clszB4", "B cluster size;pixels/cluster;B4 clusters",
	       51, -0.5, 50.5 );
  TH1D hclqB4( "clqB4", "B cluster charge;cluster charge [ke];B4 clusters",
	       100, 0, 100 );
  TH1D hncolB4( "ncolB4", "B cluster size;columns/cluster;B4 clusters",
		21, -0.5, 20.5 );
  TH1D hnrowB4( "nrowB4", "B cluster size;rows/cluster;B4 clusters",
		21, -0.5, 20.5 );
  TProfile ncolvsxmB( "ncolvsxmB",
		      "B cols vs xmod;x mod 200 [um];<cols> B4 ",
		      60, 0, 300, 0, 10 );
  TProfile nrowvsymB( "nrowvsymB",
		      "B rows vs ymod;y mod 200 [um];<rows> B4 ",
		      40, 0, 200, 0, 10 );
  // B vs ADC:

  TProfile effBvsx0( "effBvsx0", "effB vs lower x;lower ADCplet x;eff B",
		     216, 0, 64.8, -1, 2 );
  TProfile effBvsx1( "effBvsx1", "effB vs upper x;upper ADCplet x;eff B",
		     216, 0, 64.8, -1, 2 );
  TProfile effBvsy( "effBvsy", "effB vs y;ADCplet y;eff B",  81, 0, 16.2, -1, 2 );
  TProfile effBvst1( "effBvst1", "effB vs time;trigger;eff B",
		     500, 0, 1E6, -1, 2 );
  TProfile effBvst5( "effBvst5", "effB vs time;trigger;eff B",
		     500, 0, 5E6, -1, 2 );
  TProfile effBvst40( "effBvst40", "effB vs time;trigger;eff B",
		      1000, 0, 40E6, -1, 2 );
  TProfile effBvsti1( "effBvsti1", "effB vs time;trigger;iso eff B",
		      500, 0, 1E6, -1, 2 );
  TProfile effBvsti8( "effBvsti8", "effB vs time;trigger;iso eff B",
		      500, 0, 8E6, -1, 2 );
  TProfile effBvsti2( "effBvst21", "effB vs time;trigger;iso eff B",
		     400, 0, 2E6, -1, 2 );
  TProfile effBvsti10( "effBvsti10", "effB vs time;trigger;iso eff B",
		     500, 0, 10E6, -1, 2 );

  TProfile effBvswi( "effBvswi", "effB vs window;link window [mm];iso eff B",
		     20, 0.025, 1.025, -1, 2 );
  TProfile2D * effBmap1;
  effBmap1 = new TProfile2D( "effBmap1",
			     "B efficiency map;col;row;eff B",
			     8*54, 0, 8*8.1, 2*81, 0, 2*8.1, -1, 2 );
  TProfile2D * effBmap4;
  effBmap4 = new TProfile2D( "effBmap4",
			     "B efficiency map;col;row;eff B",
			     4*54, 0, 8*8.1, 1*81, 0, 2*8.1, -1, 2 );
  // triplets ADB:

  TH1D hdxADB( "dxADB", "ADB dx;x-x [mm];ADBplets", 200, -1, 1 );
  TH1D hdyADB( "dyADB", "ADB dy;y-y [mm];ADBplets", 200, -1, 1 );
  TH1D hdxcADB( "dxcADB", "ADB dx;x-x [um];ADBplets", 200, -200, 200 );
  TH1D hdycADB( "dycADB", "ADB dy;y-y [um];ADBplets", 200, -200, 200 );
  TH1D hdxciADB( "dxciADB", "ADB dx;x-x [um];isolated ADBplets",
		 200, -200, 200 );
  TH1D hdyciADB( "dyciADB", "ADB dy;y-y [um];isolated ADBplets",
		 200, -200, 200 );

  TProfile dxvsxADB( "dxvsxADB", "ADBplet dx vs x;x [mm];ADBplet <dx>",
		     216, 0, 64.8, -1, 1 );
  TProfile dxvsyADB( "dxvsyADB", "ADBplet dx vs y;y [mm];ADBplet <dx>",
		     81, 0, 16.2, -1, 1 );
  TProfile dyvsxADB( "dyvsxADB", "ADBplet dy vs x;x [mm];ADBplet <dy>",
		     216, 0, 64.8, -1, 1 );
  TProfile dyvsyADB( "dyvsyADB", "ADBplet dy vs y;y [mm];ADBplet <dy>",
		     81, 0, 16.2, -1, 1 );
  TH1D hxADB( "xADB", "ADBplets;col;ADBplets", 216, 0, 64.8 );
  TH1D hyADB( "yADB", "ADBplets;row;ADBplets",  81, 0, 16.2 );
  TH2D * hmapADB;
  hmapADB = new TH2D( "mapADB", "ADBplet map;ADBplet col;ADBplet row;ADBplets",
		      8*54, 0, 8*8.1, 2*81, 0, 2*8.1 );
  TH1D htxADB( "txADB", "ADBplet angle x;ADBplet angle x;ADBplets", 100, -1, 1 );
  TH1D htyADB( "tyADB", "ADBplet angle y;ADBplet angle y;ADBplets", 100, -1, 1 );

  // triplets BDC:

  TH1D hdxBDC( "dxBDC", "BDC dx;x-x [mm];BDCs", 200, -1, 1 );
  TH1D hdyBDC( "dyBDC", "BDC dy;y-y [mm];BDCs", 200, -1, 1 );
  TH1D hdxcBDC( "dxcBDC", "BDC dx;x-x [um];BDCs", 200, -200, 200 );
  TH1D hdycBDC( "dycBDC", "BDC dy;y-y [um];BDCs", 200, -200, 200 );
  TH1D hdxciBDC( "dxciBDC", "BDC dx;x-x [um];isolated BDCs",
		 200, -200, 200 );
  TH1D hdyciBDC( "dyciBDC", "BDC dy;y-y [um];isolated BDCs",
		 200, -200, 200 );
  TH1D hdycfBDC( "dycfBDC", "BDC dy;y-y [um];inner BDCs",
		 200, -200, 200 );
  TH1D hdycfqBDC( "dycfqBDC", "BDC dy;y-y [um];Landau peak inner BDCs",
		  200, -200, 200 );

  TProfile dxvsxBDC( "dxvsxBDC",
		     "BDC dx vs x;x [mm];BDC <dx>",
		     216, 0, 64.8, -1, 1 );
  TProfile dxvsyBDC( "dxvsyBDC",
		     "BDC dx vs y;y [mm];BDC <dx>",
		     81, 0, 16.2, -1, 1 );
  TProfile dyvsxBDC( "dyvsxBDC",
		     "BDC dy vs x;x [mm];BDC <dy>",
		     216, 0, 64.8, -1, 1 );
  TProfile dyvsyBDC( "dyvsyBDC",
		     "BDC dy vs y;y [mm];BDC <dy>",
		     81, 0, 16.2, -1, 1 );
  TProfile madyvsyBDC( "madyvsyBDC",
		       "BDC mady vs y;y [mm];BDC MAD(y) [um]",
		       81, 0, 16.2, 0, 100 );
  TProfile madyvsxBDC( "madyvsxBDC",
		       "BDC mady vs x;x [mm];BDC MAD(y) [um]",
		       216, 0, 64.8, 0, 100 );
  TProfile madyvsymBDC( "madyvsymBDC",
			"BDC mady vs ymod;y mod 200 [um];BDC MAD(y) [um]",
			40, 0, 200, 0, 100 );
  TH2D * hmapBDC;
  hmapBDC = new TH2D( "mapBDC",
		      "BDCplet map;BDCplet col;BDCplet row;BDCplets",
		      8*54, 0, 8*8.1, 2*81, 0, 2*8.1 );

  // C vs BD:

  TH1D hsizC4( "clszC4", "C cluster size;pixels/cluster;C4 clusters",
	       51, -0.5, 50.5 );
  TH1D hclqC4( "clqC4", "C cluster charge;cluster charge [ke];C4 clusters",
	       100, 0, 100 );
  TH1D hncolC4( "ncolC4", "C cluster size;columns/cluster;C4 clusters",
		21, -0.5, 20.5 );
  TH1D hnrowC4( "nrowC4", "C cluster size;rows/cluster;C4 clusters",
		21, -0.5, 20.5 );
  TProfile ncolvsxmC( "ncolvsxmC",
		      "C cols vs xmod;x mod 200 [um];<cols> C4 ",
		      60, 0, 300, 0, 10 );
  TProfile nrowvsymC( "nrowvsymC",
		      "C rows vs ymod;y mod 200 [um];<rows> C4 ",
		      40, 0, 200, 0, 10 );
  TH1D hminxC4( "minxC4", "C first pixel;first pixel mod 2;C4 clusters",
		2, -0.5, 1.5 );
  TH1D hmaxxC4( "maxxC4", "C last pixel;last pixel mod 2;C4 clusters",
		2, -0.5, 1.5 );

  // C vs ADB:

  TProfile effCvsx0( "effCvsx0", "effC vs lower x;lower ADCplet x;eff C",
		     216, 0, 64.8, -1, 2 );
  TProfile effCvsx1( "effCvsx1", "effC vs upper x;upper ADCplet x;eff C",
		     216, 0, 64.8, -1, 2 );
  TProfile effCvsy( "effCvsy", "effC vs y;ADCplet y;eff C",  81, 0, 16.2, -1, 2 );
  TProfile effCvst1( "effCvst1", "effC vs time;trigger;eff C",
		     500, 0, 1E6, -1, 2 );
  TProfile effCvst5( "effCvst5", "effC vs time;trigger;eff C",
		     500, 0, 5E6, -1, 2 );
  TProfile effCvst40( "effCvst40", "effC vs time;trigger;eff C",
		      1000, 0, 40E6, -1, 2 );
  TProfile effCvst( "effCvst", "effC vs time;trigger;eff C",
		    500, 8.3E6, 8.4E6, -1, 2 );
  TProfile effCvsti1( "effCvsti1", "effC vs time;trigger;iso eff C",
		      500, 0, 1E6, -1, 2 );
  TProfile effCvsti2( "effCvsti2", "effC vs time;trigger;iso eff C",
		      400, 0, 2E6, -1, 2 );
  TProfile effCvsti10( "effCvsti10", "effC vs time;trigger;iso eff C",
		       500, 0, 10E6, -1, 2 );
  TProfile effCvswi( "effCvswi", "effC vs window;link window [mm];iso eff C",
		     20, 0.025, 1.025, -1, 2 );
  TProfile2D * effCmap1;
  effCmap1 = new TProfile2D( "effCmap1",
			     "C efficiency map;col;row;eff C",
			     8*54, 0, 8*8.1, 2*81, 0, 2*8.1, -1, 2 );
  TProfile2D * effCmap4;
  effCmap4 = new TProfile2D( "effCmap4",
			     "C efficiency map;col;row;eff C",
			     4*54, 0, 8*8.1, 1*81, 0, 2*8.1, -1, 2 );

  TH1D hnADC( "nADC", "ADCplets;ADCplets;events", 21, -0.5, 20.5 );
  TH1D hnADB( "nADB", "ADBplets;ADBplets;events", 21, -0.5, 20.5 );
  TH1D hn4ev( "n4ev", "4plets;4plets;events", 21, -0.5, 20.5 );

 
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
	
	hcol[mod].Fill( xm );
	hrow[mod].Fill( ym );
	hmap[mod]->Fill( xm, ym );

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

	cal = adc;
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
	
	hpxq[mod].Fill( cal*ke[mod] );

	// fill pixel block for clustering
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
      
      hnpx[mod].Fill(npx);

      if( ldbg ) std::cout << std::endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[mod] = getClus();

      if( ldbg ) cout << "A clusters " << cl[mod].size() << endl;

      hncl[mod].Fill( cl[mod].size() );

      for( vector<cluster>::iterator cA = cl[mod].begin(); cA != cl[mod].end(); ++cA ) {

      hsiz[mod].Fill( cA->size );
      hclq[mod].Fill( cA->charge*ke[mod] );
      hncol[mod].Fill( cA->ncol );
      hnrow[mod].Fill( cA->nrow );

      }
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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // B vs D-C-A triplets:

    int nADC = 0;
    int nADB = 0;
    bool iso = cl[0].size() == 1 && cl[2].size() == 1 && cl[3].size() == 1;

    for( vector<cluster>::iterator cA = cl[0].begin(); cA != cl[0].end(); ++cA ) {

      double xA = cA->col*0.15 - alignx[0];
      double yA = cA->row*0.10 - aligny[0];
      double xmid = xA - 32.4;
      double ymid = yA -  8.1;
      double xrot = xmid     - ymid*fx[0] - tx[0]*xmid;
      double yrot = xmid*fy[0] + ymid     - ty[0]*ymid;
      xA = xrot + 32.4;
      yA = yrot +  8.1;
      double qA = cA->charge*ke[0];
      bool lqA = 1;
      if(      qA < 17 ) lqA = 0;
      else if( qA > 30 ) lqA = 0;

      for( vector<cluster>::iterator cD = cl[3].begin(); cD != cl[3].end(); ++cD ) {

	double xD = cD->col*0.15 - alignx[3];
	double yD = cD->row*0.10 - aligny[3];
	double xmid = xD - 32.4;
	double ymid = yD -  8.1;
	double xrot = xmid     - ymid*fx[3] - tx[3]*xmid;
	double yrot = xmid*fy[3] + ymid     - ty[3]*ymid;
	xD = xrot + 32.4;
	yD = yrot +  8.1;
	double qD = cD->charge*ke[3];
	bool lqD = 1;
	if(      qD < 17 ) lqD = 0;
	else if( qD > 30 ) lqD = 0;

	hxxDA.Fill( xA, xD );
	hyyDA.Fill( yA, yD );

	double dx2 = xD - xA;
	double dy2 = yD - yA;
	hdxDA.Fill( dx2 );
	hdyDA.Fill( dy2 );
	if( abs( dy2 ) < 0.15 && cA->big == 0 && cD->big == 0 ) {
	  hdxcDA.Fill( dx2 );
	  dxvsxDA.Fill( xD, dx2 );
	  dxvsyDA.Fill( yD, dx2 );
	}
	if( abs( dx2 ) < 0.20 && cA->big == 0 && cD->big == 0 ) {
	  hdycDA.Fill( dy2 );
	  dyvsxDA.Fill( xD, dy2 );
	  dyvsyDA.Fill( yD, dy2 );
	}

	if( abs( dx2 ) > 0.20 ) continue;
	if( abs( dy2 ) > 0.20 ) continue;

	double xavg2 = (xA + 2*xD)/3; // equidistant
	double yavg2 = (yA + 2*yD)/3; // C is closer to D: more weight for D

	double slpx = xD - xA; // angle
	double slpy = yD - yA; // angle

	for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	  double xC = cC->col*0.15 - alignx[2];
	  double yC = cC->row*0.10 - aligny[2];
	  double xmid = xC - 32.4;
	  double ymid = yC -  8.1;
	  double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
	  double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
	  xC = xrot + 32.4;
	  yC = yrot +  8.1;
	  double qC = cC->charge*ke[2];
	  bool lqC = 1;
	  if(      qC < 17 ) lqC = 0;
	  else if( qC > 30 ) lqC = 0;

	  double dx3 = xC - xavg2;
	  double dy3 = yC - yavg2;
	  hdxADC.Fill( dx3 );
	  hdyADC.Fill( dy3 );
	  if( abs( dy3 ) < 0.15 && cD->big == 0 && cA->big == 0 && cC->big == 0 ) {
	    hdxcADC.Fill( dx3*1E3 );
	    if( iso ) hdxciADC.Fill( dx3*1E3 );
	    dxvsxADC.Fill( xC, dx3 );
	    dxvsyADC.Fill( yC, dx3 );
	  }
	  if( abs( dx3 ) < 0.20 && cD->big == 0 && cA->big == 0 && cC->big == 0 ) {
	    hdycADC.Fill( dy3*1E3 );
	    if( iso ) hdyciADC.Fill( dy3*1E3 );
	    dyvsxADC.Fill( xC, dy3 );
	    dyvsyADC.Fill( yC, dy3 );
	  }
	  if( abs( dx3 ) < 0.2 && abs( dy3 ) < 0.15 ) {
	    hxADC.Fill( xC );
	    hyADC.Fill( yC );
	    hmapADC->Fill( xC, yC ); // D-C-A
	  }

	  if( abs( dx3 ) > 0.10 ) continue; // tight tri
	  if( abs( dy3 ) > 0.10 ) continue;

	  ++nADC;

	  htxADC.Fill( slpx );
	  htyADC.Fill( slpy );

	  double xavg3 = 0.5*(xA + xC); // equidistant
	  double yavg3 = 0.5*(yA + yC);
	  double xmod = fmod( xavg3+60, 0.300 )*1E3;
	  double ymod = fmod( yavg3+20, 0.200 )*1E3;

	  // efficiency of B:

	  int nm[21] = {0};

	  for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	    double xB = cB->col*0.15; // REF
	    double yB = cB->row*0.10;
	    double qB = cB->charge*ke[1];
	    bool lqB = 1;
	    if(      qB < 17 ) lqB = 0;
	    else if( qB > 30 ) lqB = 0;

	    // tri ACB:

	    double dx4 = xB - xavg3;
	    double dy4 = yB - yavg3;
	    hdxACB.Fill( dx4 );
	    hdyACB.Fill( dy4 );
	    if( abs( dy4 ) < 0.10
		&& cB->big == 0 && cD->big == 0 && cA->big == 0 && cC->big == 0 ) {
	      hdxcACB.Fill( dx4*1E3 );
	      if( iso ) hdxciACB.Fill( dx4*1E3 );
	      dxvsxACB.Fill( xB, dx4 );
	      dxvsyACB.Fill( yB, dx4 );
	    }
	    if( abs( dx4 ) < 0.10
		&& cB->big == 0 && cD->big == 0 && cA->big == 0 && cC->big == 0 ) {
	      hdycACB.Fill( dy4*1E3 );
	      if( iso ) hdyciACB.Fill( dy4*1E3 );
	      dyvsxACB.Fill( xB, dy4 );
	      dyvsyACB.Fill( yB, dy4 );
	      madyvsyACB.Fill( yB, fabs(dy4)*1E3 );
	      if( yB > 2 && yB < 13 ) { // module handle cutout
		madyvsymACB.Fill( ymod, fabs(dy4)*1E3 );
		hdycfACB.Fill( dy4*1E3 );
		if( lqA && lqB && lqC )
		  hdycfqACB.Fill( dy4*1E3 );
	      }
	    }

	    if( abs( dx4 ) < 0.20 && abs( dy4 ) < 0.20 ) {
	      hmapACB->Fill( xB, yB );
	      hsizB4.Fill( cB->size );
	      hclqB4.Fill( cB->charge*ke[1] );
	      hncolB4.Fill( cB->ncol );
	      hnrowB4.Fill( cB->nrow );
	      ncolvsxmB.Fill( xmod, cB->ncol );
	      nrowvsymB.Fill( ymod, cB->nrow );
	    }
	    
	    for( int iw = 1; iw < 21; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) { // for eff
		nm[iw] = 1;
	      }

	  } // cl B

	  effBvst1.Fill( event_nr, nm[14] );
	  effBvst5.Fill( event_nr, nm[14] );
	  effBvst40.Fill( event_nr, nm[14] );
	  if( iso ) {
	    if( yavg3 < 8.1 )
	      effBvsx0.Fill( xavg3, nm[14] );
	    else
	      effBvsx1.Fill( xavg3, nm[14] );
	    effBvsy.Fill( yavg3, nm[14] );
	    effBmap1->Fill( xavg3, yavg3, nm[14] );
	    effBmap4->Fill( xavg3, yavg3, nm[14] );
	    effBvsti1.Fill( event_nr, nm[14] );
	    effBvsti8.Fill( event_nr, nm[14] );
	    for( int iw = 1; iw < 21; ++ iw )
	      effBvswi.Fill( iw*0.050+0.005, nm[iw] );
	  }

	  effBvst1.Fill( event_nr, nm[14] );
	  effBvst5.Fill( event_nr, nm[14] );
	  effBvst40.Fill( event_nr, nm[14] );
	  if( iso ) {
	    effBvsti1.Fill( event_nr, nm[14] );
	    effBvsti2.Fill( event_nr, nm[14] );
	    effBvsti10.Fill( event_nr, nm[14] );
	  }
	  if( iso && event_nr > 999 ) {
	    if( yavg3 < 8.1 )
	      effBvsx0.Fill( xavg3, nm[14] );
	    else
	      effBvsx1.Fill( xavg3, nm[14] );
	    effBvsy.Fill( yavg3, nm[14] );
	    effBmap1->Fill( xavg3, yavg3, nm[14] );
	    effBmap4->Fill( xavg3, yavg3, nm[14] );
	    for( int iw = 1; iw < 21; ++ iw )
	      effBvswi.Fill( iw*0.050+0.005, nm[iw] );
	  }
	} // cl C
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// tri ADB:

	xavg2 = (2*xA + xD)/3; // equidistant
	yavg2 = (2*yA + yD)/3; // B is closer to A: more weight for A

	for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	  double xB = cB->col*0.15;
	  double yB = cB->row*0.10;
	  double qB = cB->charge*ke[1];
	  bool lqB = 1;
	  if(      qB < qL ) lqB = 0;
	  else if( qB > qR ) lqB = 0;

	  double dx3 = xB - xavg2;
	  double dy3 = yB - yavg2;
	  hdxADB.Fill( dx3 );
	  hdyADB.Fill( dy3 );
	  if( abs( dy3 ) < 0.15 && cD->big == 0 && cA->big == 0 && cB->big == 0 ) {
	    hdxcADB.Fill( dx3*1E3 );
	    if( iso ) hdxciADB.Fill( dx3*1E3 );
	    dxvsxADB.Fill( xB, dx3 );
	    dxvsyADB.Fill( yB, dx3 );
	  }
	  if( abs( dx3 ) < 0.20 && cD->big == 0 && cA->big == 0 && cB->big == 0 ) {
	    hdycADB.Fill( dy3*1E3 );
	    if( iso ) hdyciADB.Fill( dy3*1E3 );
	    dyvsxADB.Fill( xB, dy3 );
	    dyvsyADB.Fill( yB, dy3 );
	  }
	  if( abs( dx3 ) < 0.2 && abs( dy3 ) < 0.15 ) {
	    hxADB.Fill( xB );
	    hyADB.Fill( yB );
	    hmapADB->Fill( xB, yB ); // D-B-A
	  }

	  if( abs( dx3 ) > 0.15 ) continue; // tight tri
	  if( abs( dy3 ) > 0.15 ) continue;

	  ++nADB;

	  htxADB.Fill( slpx );
	  htyADB.Fill( slpy );

	  double xavg3 = 0.5*(xB + xD); // equidistant
	  double yavg3 = 0.5*(yB + yD);
	  double xmod = fmod( xavg3+60, 0.300 )*1E3;
	  double ymod = fmod( yavg3+20, 0.200 )*1E3;

	  // efficiency of C vs ADB:

	  int nm[21] = {0};

	  for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	    double xC = cC->col*0.15 - alignx[2];
	    double yC = cC->row*0.10 - aligny[2];
	    double xmid = xC - 32.4;
	    double ymid = yC -  8.1;
	    double xrot = xmid     - ymid*fx[2] - tx[2]*xmid;
	    double yrot = xmid*fy[2] + ymid     - ty[2]*ymid;
	    xC = xrot + 32.4;
	    yC = yrot +  8.1;

	    double qC = cC->charge*ke[2];
	    bool lqC = 1;
	    if(      qC < qL ) lqC = 0;
	    else if( qC > qR ) lqC = 0;

	    int minx = 999;
	    int maxx = 0;
	    int miny = 999;
	    int maxy = 0;

	    for( vector<pixel>::iterator p = cC->vpix.begin(); p != cC->vpix.end(); ++p ) {
	      if( p->col > maxx ) maxx = p->col;
	      if( p->col < minx ) minx = p->col;
	      if( p->row > maxy ) maxy = p->row;
	      if( p->row < miny ) miny = p->row;
	    }

	    double dx4 = xC - xavg3;
	    double dy4 = yC - yavg3;
	    hdxBDC.Fill( dx4 );
	    hdyBDC.Fill( dy4 );
	    if( abs( dy4 ) < 0.10
		&& cC->big == 0 && cD->big == 0 && cA->big == 0 && cC->big == 0 ) {
	      hdxcBDC.Fill( dx4*1E3 );
	      if( iso ) hdxciBDC.Fill( dx4*1E3 );
	      dxvsxBDC.Fill( xC, dx4 );
	      dxvsyBDC.Fill( yC, dx4 );
	    }
	    if( abs( dx4 ) < 0.10
		&& cC->big == 0 && cD->big == 0 && cA->big == 0 && cC->big == 0 ) {
	      hdycBDC.Fill( dy4*1E3 );
	      if( iso ) hdyciBDC.Fill( dy4*1E3 );
	      dyvsxBDC.Fill( xC, dy4 );
	      dyvsyBDC.Fill( yC, dy4 );
	      madyvsyBDC.Fill( yC, fabs(dy4)*1E3 );
	      if( yC > 2 && yC < 13 ) { // module handle cutout
		madyvsxBDC.Fill( xC, fabs(dy4)*1E3 );
		madyvsymBDC.Fill( ymod, fabs(dy4)*1E3 );
		hdycfBDC.Fill( dy4*1E3 );
		if( lqA && lqC && lqC )
		  hdycfqBDC.Fill( dy4*1E3 );
	      }
	    }

	    if( abs( dx4 ) < 0.20 && abs( dy4 ) < 0.20 ) {
	      hmapBDC->Fill( xC, yC );
	      if( cC->big == 0 && cD->big == 0 && cA->big == 0 && cC->big == 0 ) {
		hsizC4.Fill( cC->size );
		hclqC4.Fill( cC->charge*ke[2] );
		hncolC4.Fill( cC->ncol );
		hnrowC4.Fill( cC->nrow );
		ncolvsxmC.Fill( xmod, cC->ncol );
		nrowvsymC.Fill( ymod, cC->nrow );
		hminxC4.Fill( (minx-1)%2 ); 
		hmaxxC4.Fill( (maxx-1)%2 ); 
	      }
	    }

	    for( int iw = 1; iw < 21; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) { // for eff
		nm[iw] = 1;
	      }

	  } // cl C

	  effCvst.Fill( event_nr, nm[14] );
	  effCvst1.Fill( event_nr, nm[14] );
	  effCvst5.Fill( event_nr, nm[14] );
	  effCvst40.Fill( event_nr, nm[14] );
	  if( iso ) {
	    effCvsti1.Fill( event_nr, nm[14] );
	    effCvsti2.Fill( event_nr, nm[14] );
	    effCvsti10.Fill( event_nr, nm[14] );
	  }
	  if( iso && event_nr > 999 ) {
	    if( yavg3 < 8.1 )
	      effCvsx0.Fill( xavg3, nm[14] );
	    else
	      effCvsx1.Fill( xavg3, nm[14] );
	    effCvsy.Fill( yavg3, nm[14] );
	    effCmap1->Fill( xavg3, yavg3, nm[14] );
	    effCmap4->Fill( xavg3, yavg3, nm[14] );
	    for( int iw = 1; iw < 21; ++ iw )
	      effCvswi.Fill( iw*0.050+0.005, nm[iw] );
	  }

	} // cl B

      } // cl D

    } // cl A


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    // tracking


  } while (reader.NextEvent());

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

  return 0;
}
