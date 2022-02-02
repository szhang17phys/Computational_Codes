#include <cmath>
#include <cstdint>
#include <memory>

#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"

int main(int argc, char *argv[])
{
  // How many points?
  uint32_t N     = 100;
  double   xlow  = 1.0;
  double   delta = 1.0;

  // Decay parameters
  double numberofevents = 1000.0;
  double halflife       = 50.0;

  // Make graph
  auto graph = std::make_unique<TGraph>(N);
  for (uint32_t i = 0; i < N; ++i) {
    double x = xlow + i*delta;
    graph->SetPoint(i, x, numberofevents*std::exp(-x*std::log(2)/halflife));
  }

  // Draw graph on canvas
  auto canvas = std::make_unique<TCanvas>("c", "c", 600, 600);

  graph->GetXaxis()->SetTitle("Time [s]");
  graph->GetYaxis()->SetTitle("Number Events");
  graph->Draw("ape");

  canvas->SaveAs("decay.png");
  
  return 0;
}
