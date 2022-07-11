#ifndef HGHWAY_LABELING_H_
#define HGHWAY_LABELING_H_

#include <stdint.h>
#include <iostream>
#include <thread>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <utility>

#include "two_layer_queue.h"

//
// NOTE: Currently only unweighted and undirected graphs are supported.
//

class HighwayLabelling {
 public:
  // Constructs labelling from a graph, given as a list of edges.
  HighwayLabelling(const char *filename, int k);
  HighwayLabelling();
  ~HighwayLabelling();

  void ConstructHighwayLabelling(int topk[]);

  void SelectLandmarks_HD(int topk[]);

  long LabellingSize();
  uint8_t min(uint8_t a, uint8_t b);

  // Returns distance vetween vertices v and w if they are connected.
  // Otherwise, returns 99.
  uint8_t QueryDistanceUB_naive(int s, int t);
  uint8_t QueryDistanceUB_opt(int s, int t);
  uint8_t QueryDistance(int s, int t);

  bool IsLandmark(int arr[], int left, int right, int value);
  void RemoveLandmarks(int landmarks[]);

  // Loads the highway labelling.
  void LoadIndex(const char *filename);

  // Stores the highway labelling.
  void StoreIndex(const char *filename);

 private:
  int V;  // total number of vertices
  long E; // total number of edges
  int K; // total number of landmarks

  uint8_t **vertices;
  uint8_t **distances;
  uint8_t **highway;
  uint8_t *C;

  std::unordered_map<int, std::vector<int> > adj;
};

HighwayLabelling::HighwayLabelling() { }

HighwayLabelling::~HighwayLabelling() {

  for(int i = 0; i < V; i++) {
    delete [] distances[i];
  }
  delete [] distances;

  if(vertices != NULL) {
    for(int i = 0; i < V; i++) {
      delete [] vertices[i];
    }
    delete [] vertices;
  }
  delete [] C;

  for(int i = 0; i < K; i++)
    delete [] highway[i];
  delete [] highway;

}

HighwayLabelling::HighwayLabelling(const char *filename, int k) {
  V = 0; E = 0; K = k;

  std::ifstream ifs(filename);
  if (ifs.is_open()) {

    int v, w; std::string query;
    std::unordered_map<int, int> vertex2id;
    while (getline(ifs, query)){
      std::istringstream iss(query);
      iss >> v >> w;
      
      if (vertex2id.count(v) == 0) vertex2id[v] = V++;
      if (vertex2id.count(w) == 0) vertex2id[w] = V++;
      v = vertex2id[v];
      w = vertex2id[w];
      if (v != w) {
        adj[v].push_back(w);
        adj[w].push_back(v);
      }
    }
    ifs.close();

    for (int v = 0; v < V; v++) {
      std::sort(adj[v].begin(), adj[v].end());
      adj[v].erase(std::unique(adj[v].begin(), adj[v].end()), adj[v].end());
    }

    for(int i = 0; i < V; i++)
      E += adj[i].size();
    std::cout << "V : " << V << " E : " << E << std::endl << std::endl;

  } else
      std::cout << "Unable to open file" << std::endl;
}

void HighwayLabelling::RemoveLandmarks(int landmarks[]) {

  for(int i = 0; i < K; i++) {
    for (std::vector<int>::iterator it=adj[landmarks[i]].begin(); it!=adj[landmarks[i]].end(); ++it)
      adj[*it].erase(std::find(adj[*it].begin(), adj[*it].end(), landmarks[i]));
    adj.erase(landmarks[i]);
  }

}

long HighwayLabelling::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        size++;
    }
  }

  return size;
}

bool HighwayLabelling::IsLandmark(int arr[], int left, int right, int value) {

  while (left <= right) {
    int middle = (left + right) / 2;
    if (arr[middle] == value)
      return true;
    else if (arr[middle] > value)
      right = middle - 1;
    else
      left = middle + 1;
  }
  return false;
}

void HighwayLabelling::ConstructHighwayLabelling(int topk[]) {

  // Initialization
  distances = new uint8_t*[V];
  for(int i = 0; i < V; i++) {
    distances[i] = new uint8_t[K];
    for(int j = 0; j < K; j++)
      distances[i][j] = 111;
  }

  highway = new uint8_t*[K];
  for(int i = 0; i < K; i++)
    highway[i] = new uint8_t[K];

  // Start computing Highway Labelling (HL)
  for (int i = 0; i < K; i++) {
    uint8_t *P = new uint8_t[V];
    for(int j = 0; j < V; j++)
      P[j] = 99;

    std::queue<int> que[2];

    que[0].push(topk[i]); que[0].push(-1);
    distances[topk[i]][i] = 0; P[topk[i]] = 0; int use = 0;
    while (!que[0].empty()) {
      int u = que[use].front();
      que[use].pop();

      if(u == -1) {
        use = 1 - use;
        que[use].push(-1);
        continue;
      }

      for (int w : adj[u]) {
        if (P[w] == 99) {
          P[w] = P[u] + 1;
          if(use == 1 || IsLandmark(topk, 0, K, w))
            que[1].push(w);
          else {
            que[0].push(w);
            distances[w][i] = P[w];
          }
        }
      }
    }

    for(int j = 0; j < K; j++) {
      if(P[topk[j]] != 111) {
        highway[i][j] = P[topk[j]];
        highway[j][i] = P[topk[j]];
      }
    }

    delete [] P;
  }

  std::cout << "Labelling Size : " << LabellingSize() << std::endl;
}

void HighwayLabelling::SelectLandmarks_HD(int topk[]) {
  std::vector<std::pair<int, int> > deg(V );
  for (int v = 0; v < V; v++)
    deg[v] = std::make_pair(adj[v].size(), v);

  std::sort(deg.rbegin(), deg.rend());

  for (int v = 0; v < K; v++)
    topk[v] = deg[v].second;
}

uint8_t HighwayLabelling::min(uint8_t a, uint8_t b) {
  return (a < b) ? a : b;
}

uint8_t HighwayLabelling::QueryDistanceUB_naive(int s, int t) {

  uint8_t m = 99; int i, j;
  for(i = 0; i < C[s]; i++) {
    for (j = 0; j < C[t]; j++)
      m = min(m, distances[s][i] + highway[vertices[s][i]][vertices[t][j]] + distances[t][j]);
  }

  return m;
}

uint8_t HighwayLabelling::QueryDistanceUB_opt(int s, int t) {

  uint8_t m = 99, uni1[C[s]] = {0}, uni2[C[t]] = {0}; int i = 0, j = 0, i1 = 0, j1 = 0;
  while (i < C[s] && j < C[t]) {
    if (vertices[s][i] < vertices[t][j]) {
      uni1[i1] = i; i++; i1++;
    } else if (vertices[t][j] < vertices[s][i]) {
      uni2[j1] = j; j++; j1++;
    } else {
      m = min(m, distances[s][i] + distances[t][j]);
      i++; j++;
    }
  }

  while (i < C[s]) {
    uni1[i1] = i; i++; i1++;
  }

  while (j < C[t]) {
    uni2[j1] = j; j++; j1++;
  }

  i = 0;
  while (i < i1) { j = 0;
    while (j < j1) {
      m = min(m, distances[s][uni1[i]] + highway[vertices[s][uni1[i]]][vertices[t][uni2[j]]] + distances[t][uni2[j]]);
      j++;
    }
    i++;
  }

  return m;
}

uint8_t HighwayLabelling::QueryDistance(int s, int t) {
  if (s == t) {   return 0;  }

  std::vector<TwoLayerQueue> qque;
  std::vector<uint8_t> qdist[2];

  qdist[0].resize(V, 99);
  qdist[1].resize(V, 99);

  qque.push_back(TwoLayerQueue(V));
  qque.push_back(TwoLayerQueue(V));

  for (int dir = 0; dir < 2; dir++){
    int v = dir == 0 ? s : t;
    qque[dir].clear();
    qque[dir].push(v);
    qque[dir].next();
    qdist[dir][v] = 0;
  }

  uint8_t dist_upper = QueryDistanceUB_opt(s, t);
  uint8_t res = dist_upper, dis[2] = {0, 0};
  while (!qque[0].empty() && !qque[1].empty()) {
    int use = 0;
    use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
    dis[use]++;

    if (dis[0] + dis[1] == dist_upper) {
      res = dis[0] + dis[1];
      goto LOOP_END;
    }

    while (!qque[use].empty()) {

      int v = qque[use].front();
      qque[use].pop();

      for (int w : adj[v]) {

        uint8_t &src_d = qdist[use][w];
        uint8_t &dst_d = qdist[1 - use][w];
        if (src_d != 99) continue;
        if (dst_d != 99) {
          res = qdist[use][v] + 1 + dst_d;
          goto LOOP_END;
        } else {
          qque[use].push(w);
          qdist[use][w] = qdist[use][v] + 1;
        }
      }
    }
    qque[use].next();
  }
  LOOP_END:

  for (int dir = 0; dir < 2; dir++) {
    for (int v : qque[dir]) {
      qdist[dir][v] = 99;
    }
    qque[dir].clear();
  }

  return min(res, dist_upper);
}

void HighwayLabelling::LoadIndex(const char *filename) {
  std::ifstream ifs(std::string(filename) + std::string("_index"));
  
  C = new uint8_t[V];
  vertices = new uint8_t*[V];
  distances = new uint8_t*[V];

  long size = 0;
  for(int i = 0; i < V; i++) {
    ifs.read((char*)&C[i], sizeof(C[i]));
    size += C[i];
    vertices[i] = new uint8_t[C[i]];
    distances[i] = new uint8_t[C[i]];
    for(uint8_t j = 0; j < C[i]; j++) {
      ifs.read((char*)&vertices[i][j], sizeof(vertices[i][j]));
      ifs.read((char*)&distances[i][j], sizeof(distances[i][j]));
    }
  }
  std::cout << "Labelling Size : " << size << std::endl;
  ifs.close();

  ifs.open(std::string(filename) + std::string("_highway"));

  highway = new uint8_t*[K];
  for(uint8_t i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    for(uint8_t j = 0; j < K; j++)
      ifs.read((char*)&highway[i][j], sizeof(highway[i][j]));
  }
  ifs.close();
}

void HighwayLabelling::StoreIndex(const char *filename) {
  std::ofstream ofs(std::string(filename) + std::string("_index"));

  for(int i = 0; i < V; i++) {
    uint8_t C = 0;
    for(int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        C++;
    }

    ofs.write((char*)&C, sizeof(C));
    for(uint8_t j = 0; j < K; j++) {
      if(distances[i][j] != 111) {
        ofs.write((char*)&j, sizeof(j));
        ofs.write((char*)&distances[i][j], sizeof(distances[i][j]));
      }
    }
  }
  ofs.close();

  ofs.open(std::string(filename) + std::string("_highway"));
  for(int i = 0; i < K; i++) {
    for(int j = 0; j < K; j++) {
      ofs.write((char*)&highway[i][j], sizeof(highway[i][j]));
    }
  }
  ofs.close();
}

#endif  // HGHWAY_LABELING_H_
