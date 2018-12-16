#include <iostream>
#include "highway_cover_labelling.h"

using namespace std;

int main(int argc, char **argv) {

  int k = atoi(argv[2]);
  HighwayLabelling *hl = new HighwayLabelling(argv[1], k);
  hl->LoadIndex(argv[3]);

  int topk[k];
  hl->SelectLandmarks_HD(topk);
  hl->RemoveLandmarks(topk);

  // query distance
  for (int s, t; cin >> s >> t; ) {
    cout << (int) hl->QueryDistance(s, t) << endl;
  }

  exit(EXIT_SUCCESS);
}
