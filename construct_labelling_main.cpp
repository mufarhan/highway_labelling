#include <iostream>
#include "highway_cover_labelling.h"

using namespace std;

int main(int argc, char **argv) {
  int k = atoi(argv[2]);
  HighwayLabelling *hl = new HighwayLabelling(argv[1], k);

  int topk[k];
  hl->SelectLandmarks_HD(topk);
  sort(topk, topk + k);

  // construct labelling
  hl->ConstructHighwayLabelling(topk);
  hl->StoreIndex(argv[3]);

  delete hl;
  
  exit(EXIT_FAILURE);
}
