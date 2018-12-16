# highway_labelling (HL)

This is an implementation of a Highway Labelling (HL) algorithm, which answers shortest-path distance queries over real-world massive complex networks in the order of milliseconds(ms).

## Usage
Given a graph, it first constructs a labelling. Then, using the labelling and a querying framework, it can quickly answer distance between two vertices.

### CUI Interface
We can follow the following instructions in order to use our CUI:

    $ make
    $ bin/construct_index graph_file_name num_of_landmarks index_file_name
    $ bin/query_distance graph_file_name num_of_landmarks index_file_name <<< "1 2"


### From Our Implementation

We implement the following methods to answer distance queries:

* Call `ConstructHighwayLabelling` to construct a highway cover labelling from a graph (a sample graph file shown in sample folder).
* Call `QueryDistance` to query distance between two vertices.
* Call `LoadIndex` and `StoreIndex` to load and store an index.

For further information, please see `highway_cover_labelling.h`.

## References

* Muhammad Farhan, Qing Wang, Yu Lin, and Brendan Mckay, **[A Highly Scalable Labelling Approach for Exact Distance
Queries in Complex Networks](https://arxiv.org/pdf/1812.02363.pdf)**.
