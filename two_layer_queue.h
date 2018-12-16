#ifndef TWO_LAYER_QUEUE_H
#define TWO_LAYER_QUEUE_H

#include <assert.h>
#include <iostream>
#include <vector>

class TwoLayerQueue {
  std::vector<int> data;
  int V;
  int curr_;
  int next_;
  int end_;

public:
  TwoLayerQueue(int V) : data(V), V(V), curr_(0), next_(0), end_(0){ }
  bool empty() const { return curr_ == next_; }
  bool full() const { return end_ == V; }
  int &front() { return data[curr_];}
  int size() const { return end_; }
  void pop() { ++curr_; assert(curr_ <= end_);}
  void push(const int &val){ data[end_++] = val; assert(end_ <= V);}
  void next() { assert(curr_ == next_); next_ = end_; }
  void clear() { curr_ = next_ = end_ = 0; }

  std::vector<int>::iterator begin() { return data.begin();}
  std::vector<int>::iterator end() { return data.begin() + end_;}
};

#endif /* TWO_LAYER_QUEUE_H */
