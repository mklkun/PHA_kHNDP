#ifndef GRAPH_H
#define GRAPH_H

#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>


using namespace std;


enum colors { BLACK, GREY, RED, ORANGE, YELLOW, GREEN, BLUE, WHITE};


 using namespace lemon;


    typedef lemon::ListGraph lemonGraph;
    typedef lemon::ListDigraph lemonDigraph;

    typedef lemonGraph::EdgeIt lemonEdgeIt;
    typedef lemonGraph::NodeIt lemonNodeIt;
    typedef lemonDigraph::ArcIt lemonArcIt;
    typedef lemonDigraph::NodeIt lemonDiNodeIt;

    typedef lemonGraph::Node lemonNode;
    typedef lemonGraph::Edge lemonEdge;
    typedef lemonDigraph::Node lemonDiNode;
    typedef lemonDigraph::Arc lemonArc;

    // une Map sert à donner des valeurs aux composants (noeuds, arêtes, ...) d'un graphe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    typedef lemonGraph::NodeMap<int> lemonNodeIndexMap;
    typedef lemonGraph::NodeMap<int> lemonNodeDegreeMap;
    typedef lemonGraph::NodeMap<int> lemonNodeValMap;

    typedef lemonGraph::EdgeMap<int> lemonEdgeIndexMap;
    typedef lemonGraph::EdgeMap<int> lemonEdgeCapacityMap;
    typedef lemonGraph::EdgeMap<double> lemonEdgeWeightMap;
    typedef lemonGraph::EdgeMap<int> lemonEdgeValMap;

    typedef lemonGraph::EdgeMap<float> lemonEdgeLambdaMap;


    typedef lemonDigraph::NodeMap<int> lemonDiNodeIndexMap;
    typedef lemonDigraph::NodeMap<int> lemonDiNodeDegreeMap;
    typedef lemonDigraph::NodeMap<int> lemonDiNodeValMap;

    typedef lemonDigraph::ArcMap<int> lemonArcIndexMap;
    typedef lemonDigraph::ArcMap<int> lemonArcCapacityMap;
    typedef lemonDigraph::ArcMap<double> lemonArcWeightMap;
    typedef lemonDigraph::ArcMap<int> lemonArcValMap;

    typedef lemonDigraph::ArcMap<float> lemonArcLambdaMap;


    using lemon::INVALID;



struct graph_str
{
    string name;

    vector<lemonNode> &LNodes;
    vector<lemonEdge> &LEdges;

    lemonNodeIndexMap &NIM;
    lemonNodeDegreeMap &NDM;
    lemonNodeValMap &NVM;

    lemonEdgeIndexMap &EIM;
    lemonEdgeCapacityMap &ECM;
    lemonEdgeWeightMap &EWM;
    lemonEdgeValMap &EVM;

    map <int, int> linking;
}; typedef struct graph_str GRAPH_STR;



struct digraph_str
{
    vector<lemonDiNode> &LDNodes1;
    vector<lemonDiNode> &LDNodes2;
    vector<lemonArc> &LArcs;
    lemonDiNodeIndexMap &DNIM;
    lemonDiNodeDegreeMap &DNDM;
    lemonDiNodeValMap &DNVM;

    lemonArcIndexMap &AIM;
    lemonArcCapacityMap &ACM;
    lemonArcWeightMap &AWM;
    lemonArcValMap &AVM;

    map <int, int> linking;

    lemonDiNode s;
    lemonDiNode t;

    unsigned n;
    unsigned m;
}; typedef struct digraph_str DIGRAPH_STR;



    template <typename GRAPH, typename GRAPH_STRUCT>
    void write_graph_viz(GRAPH const &graph, GRAPH_STRUCT const grapht_struct, unsigned const n, unsigned const m, string name);

    template <typename GRAPH, typename GRAPH_STRUCT>
    void write_graph_viz_sol(GRAPH const &graph, GRAPH_STRUCT const grapht_struct, unsigned const n, unsigned const m, string name, vector < vector <int> > const demands);

    template <typename GRAPH, typename GRAPH_STRUCT>
    void write_graph(GRAPH const &graph, GRAPH_STRUCT const &grapht_struct);

    template <typename GRAPH, typename GRAPH_STRUCT>
    void write_graph(GRAPH const &graph, GRAPH_STRUCT const &grapht_struct,string s);



class graph
{
    public:
        graph();
        virtual ~graph();
    protected:
    private:
};

#endif // GRAPH_H
