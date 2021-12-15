#include "../include/graph.h"

#include <lemon/lgf_writer.h>


template <>
void write_graph_viz(lemonGraph const &g, GRAPH_STR const g_str, unsigned const n, unsigned const m, string name)
{
    ofstream fichier(name.c_str(), ios::out | ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    if(fichier)
    {
        fichier << "graph G {" << endl;

        for(unsigned i(0); i < n; i++)
            fichier << i << " ;" << endl;                                                                                  /// ADDING NODE PROPERTIES

        for(unsigned i(0); i < m; i++)
            fichier << g.id(g.u(g.edgeFromId(i))) << "--" << g.id(g.v(g.edgeFromId(i))) << " ;" << endl;                   /// ADDING EDGE PROPERTIES

        fichier << "}" << endl;
        fichier.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;
    //for ()
    //cout << " [index=\"" << g_str.NIM[v] << "\", degree=\"" << g_str.NDM[e] << "\"]";

    //cout << " [index=\"" << g_str.EIM[e] << "\", weight=\"" << g_str.EWM[e] << "\", capacity=\"" << g_str.ECM[e] << "\"]";
}


template <>
void write_graph_viz(lemonDigraph const &g, GRAPH_STR const g_str, unsigned const n, unsigned const m, string name)
{
    ofstream fichier(name.c_str(), ios::out | ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    if(fichier)
    {
        fichier << "graph G {" << endl;

        for(unsigned i(0); i < n; i++)
            fichier << i << " ;" << endl;                                                                                  /// ADDING NODE PROPERTIES

        for(unsigned i(0); i < m; i++)
            fichier << g.id(g.source(g.arcFromId(i))) << "->" << g.id(g.target(g.arcFromId(i))) << " ;" << endl;                   /// ADDING EDGE PROPERTIES

        fichier << "}" << endl;
        fichier.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;
    //for ()
    //cout << " [index=\"" << g_str.NIM[v] << "\", degree=\"" << g_str.NDM[e] << "\"]";

    //cout << " [index=\"" << g_str.EIM[e] << "\", weight=\"" << g_str.EWM[e] << "\", capacity=\"" << g_str.ECM[e] << "\"]";
}


bool exists_in_position(int const indice, vector < vector <int> > const demands, int const pos)
{
    for(unsigned i(0); i < demands.size(); i++)
        if (demands[i][pos] == indice)
            return(true);
    return(false);
}


template <>
void write_graph_viz_sol(lemonGraph const &g, GRAPH_STR const g_str, unsigned const n, unsigned const m, string name, vector < vector <int> > const demands)
{
    ofstream fichier(name.c_str(), ios::out | ios::trunc);  // ouverture en écriture avec effacement du fichier ouvert

    if(fichier)
    {
        fichier << "graph G {" << endl;

        for(unsigned i(0); i < n; i++)
        {
            if (exists_in_position(i,demands,0))
                fichier << i << " [color=blue];" << endl;
            else if (exists_in_position(i,demands,1))
                fichier << i << " [color=red];" << endl;
            else
                fichier << i << " ;" << endl;
        }

        for(unsigned i(0); i < m; i++)
            if(g_str.EVM[g.edgeFromId(i)] == 1)
                fichier << g.id(g.u(g.edgeFromId(i))) << "--" << g.id(g.v(g.edgeFromId(i))) << " ;" << endl;

        fichier << "}" << endl;
        fichier.close();
    }
    else
        cerr << "Impossible d'ouvrir le fichier !" << endl;
}


template <>
void write_graph(lemonGraph const &g, GRAPH_STR const &g_str)
{
    graphWriter(g, cout)
    .nodeMap("Index", g_str.NIM)
    .nodeMap("Degree", g_str.NDM)
    .edgeMap("Index", g_str.EIM)
    .edgeMap("Capacity", g_str.ECM)
    .edgeMap("Weight", g_str.EWM)
    .edgeMap("Value", g_str.EVM)
    .run();
}


template <>
void write_graph(lemonDigraph const &g, DIGRAPH_STR const &dg_str)
{
    digraphWriter(g, cout)
    .nodeMap("Index", dg_str.DNIM)
    .nodeMap("Degree", dg_str.DNDM)
    .nodeMap("Value", dg_str.DNVM)
    .arcMap("Index", dg_str.AIM)
    .arcMap("Capacity", dg_str.ACM)
    .arcMap("Weight", dg_str.AWM)
    .arcMap("Value", dg_str.AVM)
    //.arcMap("lambda", dg_str.ALM)
    .run();
}


template <>
void write_graph(lemonGraph const &g, GRAPH_STR const &g_str, string s)
{
    ofstream filename(s.c_str());

    graphWriter(g, filename)
    .nodeMap("Index", g_str.NIM)
    .nodeMap("Degree", g_str.NDM)
    .edgeMap("Index", g_str.EIM)
    .edgeMap("Capacity", g_str.ECM)
    .edgeMap("Weight", g_str.EWM)
    .run();
}


template <>
void write_graph(lemonDigraph const &g, DIGRAPH_STR const &dg_str, string s)
{
    ofstream filename(s.c_str());

    digraphWriter(g, filename)
    .nodeMap("Index", dg_str.DNIM)
    .nodeMap("Degree", dg_str.DNDM)
    .nodeMap("Value", dg_str.DNVM)
    .arcMap("Index", dg_str.AIM)
    .arcMap("Capacity", dg_str.ACM)
    .arcMap("Weight", dg_str.AWM)
    .arcMap("Value", dg_str.AVM)
    .run();
}


graph::graph()
{
    //ctor
}

graph::~graph()
{
    //dtor
}

