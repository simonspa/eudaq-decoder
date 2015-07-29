#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include "TH1D.h"

using namespace eudaq;

int main() {

  FileReader reader = FileReader("20202", "run$6R$X");

  TH1D * h1 = new TH1D();
  
  size_t event_nr = 0;
  
  do {
    // Get next event:
    DetectorEvent evt = reader.GetDetectorEvent();

    if (evt.IsBORE()) { eudaq::PluginManager::Initialize(evt); }
    
    if (event_nr%1000==0) {
      std::cout<<"Processing event "<< event_nr<<std::endl;

      StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

      for (size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane) {
	const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);
	std::vector<double> pxl = plane.GetPixels<double>();

	std::cout << "PLANE " << plane.ID() << ": ";
	for (size_t ipix = 0; ipix < pxl.size(); ++ipix) {
	  std::cout << plane.GetX(ipix) << " " << plane.GetY(ipix) << " " << plane.GetPixel(ipix) << " ";
	}
	std::cout << std::endl;
      }
    }
    
    event_nr++;
      
  } while (reader.NextEvent());

  return 0;
}
