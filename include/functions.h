#ifndef FUNCTIONS_INCLUDED
#define FUNCTIONS_INCLUDED

#include <vector>
#include <deque>

#include "../include/graph.h"


using namespace std;


//template <typename GRAPH, typename DIGRAPHS, typename G_STR, typename DG_STRS>
//void Lagrangien_Thread(GRAPH const &G, DIGRAPHS const &DGS, G_STR const &g_str, DG_STRS const &dg_strs, vector < vector <int> > const &demands, int const k);

//template <typename GRAPH, typename DIGRAPHS, typename G_STR, typename DG_STRS>
//void Successive_Thread(GRAPH const &G, DIGRAPHS const &DGS, G_STR const &g_str, DG_STRS const &dg_strs, vector < vector <int> > const &demands, int const k);

template <typename GRAPH>
bool open_graph(GRAPH &graph, string name);

template <typename GRAPH, typename G_STR>
bool open_bundled_graph(GRAPH &graph, vector < vector <int> > &demands, unsigned &n, unsigned &m, G_STR &g_str, string name);

template <typename GRAPH, typename G_STR>
bool open_TSP_graph(GRAPH &graph, vector < vector <int> > &demands, unsigned &n, unsigned &m, G_STR &g_str, unsigned &k, string name);

template <typename GRAPH, typename G_STR>
bool open_SOL_graph(GRAPH &graph, G_STR &g_str, string name);

template <typename GRAPH>
bool open_graph(GRAPH &graph, char *name, int k);

template <typename GRAPH, typename POPULATION, typename POPRANKING>
void write_population(GRAPH const &G, POPULATION const &Population, POPRANKING const &PopRanking, string const name);

template <typename OLDGRAPH, typename NEWGRAPH, typename G_STR, typename DG_STR>
void transform_graph(OLDGRAPH const &OG, NEWGRAPH *NG, G_STR g_str, DG_STR* dg_str, int const s_i, int const t_i);

template <typename GRAPH, typename DIGRAPHS, typename G_STR, typename DG_STRS, typename VTAVM, typename TEVM>
bool is_feasable(GRAPH const &G, DIGRAPHS const &DGS, G_STR const &g_str, DG_STRS const &dg_strs, VTAVM const &VARVM, TEVM const &EVM);

template <typename DIGRAPH_STRUCT>
void MyShortPaths(lemonDigraph const &graph, DIGRAPH_STRUCT *dg_str);

template <typename DIGRAPH_STRUCT>
void MyMinCostFlow(lemonDigraph const &graph, DIGRAPH_STRUCT *dg_str);

void test_HA1(vector < vector <int> > const &demands, string const name);

void Selection(unsigned const l,unsigned &p1, unsigned &p2, unsigned &alchimy, vector <string> &pairs);

vector <pair <unsigned, unsigned> > Crossover(unsigned const l , unsigned const alchimy);

vector <pair <unsigned, unsigned> > CrossoverFull(unsigned const l);


#endif // FUNCTIONS_INCLUDED
