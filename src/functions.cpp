#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

#include "../include/functions.h"
#include "../include/graph.h"

#include <lemon/lgf_writer.h>



using namespace std;



bool extract_number(string &str, int &number)
{
    string temp;

    unsigned i=0;

    while ( i < str.size() )
    {
        if(isdigit(str[i]))
        {
            while (i < str.size() && (isdigit(str[i]) || str[i] == '.'))
            {
                temp += str[i];
                i++;
            }
/*            else if (str[i] == '.')
            {
                str.erase(str.begin() + i);
                while(isdigit(str[i]) && str[i] == '0')
                    str.erase(str.begin() + i);
            }*/
            istringstream stream(temp);
            stream >> number;

            if(i == str.size())
                str.clear();
            else
                str.erase( 0, i);

            return(true);
        }

        i++;
    }

    str.clear();
    return(false);
}



template <>
bool open_TSP_graph(lemonGraph &g, vector < vector <int> > &demands, unsigned &n, unsigned &m, GRAPH_STR &g_str, unsigned &k, string name)
{
    string s = "../../Instances/New/" + name + ".tsp";
    ifstream fichier(s.c_str(), ios::in);


    if(fichier)  // si l'ouverture a fonctionné
    {
        int number;
        vector <int> temp;
        vector < vector <int> > v;
        string ligne;  // déclaration d'une chaîne qui contiendra la ligne lue
        unsigned flag(0);

        while (getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
        {
            if(flag == 1)
            {
                while (!ligne.empty())
                    if (extract_number(ligne, number))
                        temp.push_back(number);

                if (temp.size() == 3)
                    v.push_back(vector <int> (temp));
                else if(temp.size() == 0)
                    flag = 0;
            }
            else if (flag == 2)
            {
                while (!ligne.empty())
                    if (extract_number(ligne, number))
                        temp.push_back(number);

                if (temp.size() == 3 && temp[0]!=temp[1])
                {
                    k = temp[2];
                    temp.erase(temp.begin() + 2);
                    demands.push_back(vector <int> (temp));
                }
            }
            else
            {
                if(!ligne.compare("NODE_COORD_SECTION"))
                    flag = 1;
                else if (!ligne.compare("DEMAND_SECTION"))
                    flag = 2;
            }

            temp.clear();
        }

        fichier.close();

        for (unsigned i=0; i < v.size(); i++)
        {
            g_str.LNodes.push_back(g.addNode());

            g_str.linking[v[i][1]] = i;
            g_str.NIM[g_str.LNodes[g_str.LNodes.size()-1]] = i;
            g_str.NDM[g_str.LNodes[g_str.LNodes.size()-1]] = v[i][1];
            g_str.NVM[g_str.LNodes[g_str.LNodes.size()-1]] = v[i][2];

            n++;
        }

        for (unsigned i=0; i < v.size(); i++)
        {
            for (unsigned j=i+1; j < v.size(); j++)
            {
                g_str.LEdges.push_back(g.addEdge( g_str.LNodes[i] , g_str.LNodes[j]));

                g_str.EIM[g_str.LEdges[g_str.LEdges.size()-1]] = m;
                g_str.ECM[g_str.LEdges[g_str.LEdges.size()-1]] = 1;

                double dist = (g_str.NDM[g_str.LNodes[i]] - g_str.NDM[g_str.LNodes[j]]) * (g_str.NDM[g_str.LNodes[i]] - g_str.NDM[g_str.LNodes[j]]);
                dist += (g_str.NVM[g_str.LNodes[i]] - g_str.NVM[g_str.LNodes[j]]) * (g_str.NVM[g_str.LNodes[i]] - g_str.NVM[g_str.LNodes[j]]);
                dist = sqrt(dist);

                g_str.EWM[g_str.LEdges[g_str.LEdges.size()-1]] = dist;

                m++;
            }
        }

        for(unsigned i=0; i < demands.size(); i++)
        {
            demands[i][0]--;
            demands[i][1]--;
        }
    }
    else
        return(false);
    return(true);
}



template <>
bool open_bundled_graph(lemonGraph &g, vector < vector <int> > &demands, unsigned &n, unsigned &m, GRAPH_STR &g_str, string name)
{
    string s = "input/weighted_" + name + ".dot";
    ifstream fichier(s.c_str(), ios::in);

    if(fichier)  // si l'ouverture a fonctionné
    {
        int number;
        vector <int> temp;
        vector < vector <int> > e,v;
        string ligne;  // déclaration d'une chaîne qui contiendra la ligne lue
        while (getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
        {
            while (!ligne.empty())
                if (extract_number(ligne, number))
                    temp.push_back(number);

            if (temp.size() == 3)
                v.push_back(vector <int> (temp));
            else if (temp.size() == 5 && temp[0]!=temp[1])
                e.push_back(vector <int> (temp));
            else if (temp.size() == 2 && temp[0]!=temp[1])
                demands.push_back(vector <int> (temp));

            temp.clear();
        }

        fichier.close();

        for (unsigned i=0; i < v.size(); i++)
        {
            g_str.LNodes.push_back(g.addNode());

            g_str.linking[v[i][1]] = i;
            g_str.NIM[g_str.LNodes[g_str.LNodes.size()-1]] = v[i][1];
            g_str.NDM[g_str.LNodes[g_str.LNodes.size()-1]] = v[i][2];

            n++;
        }

        for (unsigned i=0; i < e.size(); i++)
        {
            g_str.LEdges.push_back(g.addEdge( g.nodeFromId(e[i][0]) , g.nodeFromId(e[i][1])));

            g_str.EIM[g_str.LEdges[g_str.LEdges.size()-1]] = e[i][2];
            g_str.ECM[g_str.LEdges[g_str.LEdges.size()-1]] = e[i][4];
            g_str.EWM[g_str.LEdges[g_str.LEdges.size()-1]] = e[i][3];

            m++;
        }
    }
    else
        return(false);
    return(true);
}



template <>
bool open_SOL_graph(lemonGraph &g, GRAPH_STR &g_str, string name)
{
    string s = "Solutions/" + name + "_solution.dot";
    ifstream fichier(s.c_str(), ios::in);

    if(fichier)  // si l'ouverture a fonctionné
    {
        int number;
        vector <int> temp;
        vector < vector <int> > e,v;
        string ligne;  // déclaration d'une chaîne qui contiendra la ligne lue
        while (getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
        {
            while (!ligne.empty())
                if (extract_number(ligne, number))
                    temp.push_back(number);

            if (temp.size() == 1)
                v.push_back(vector <int> (temp));
            else if (temp.size() == 4 && temp[0]!=temp[1])
                e.push_back(vector <int> (temp));


            temp.clear();
        }

        fichier.close();

        for (unsigned i=0; i < v.size(); i++)
        {
            g_str.LNodes.push_back(g.addNode());

            g_str.linking[v[i][1]] = i;
            g_str.NIM[g_str.LNodes[g_str.LNodes.size()-1]] = v[i][1];
            g_str.NDM[g_str.LNodes[g_str.LNodes.size()-1]] = v[i][2];

            //n++;
        }

        for (unsigned i=0; i < e.size(); i++)
        {
            g_str.LEdges.push_back(g.addEdge( g.nodeFromId(e[i][0]) , g.nodeFromId(e[i][1])));

            g_str.EIM[g_str.LEdges[g_str.LEdges.size()-1]] = e[i][2];
            g_str.ECM[g_str.LEdges[g_str.LEdges.size()-1]] = e[i][4];
            g_str.EWM[g_str.LEdges[g_str.LEdges.size()-1]] = e[i][3];

            //m++;
        }
    }
    else
        return(false);
    return(true);
}




template <>
void write_population(lemonGraph const &G,vector < lemonEdgeValMap* > const &Population, vector < int > const &PopRanking, string const name)
{
    unsigned m(countEdges(G));
    ofstream fichier(name.c_str(), ios::out | ios::trunc);
    for (unsigned i(0); i<Population.size(); i++)
    {
        fichier << endl << PopRanking[i] << endl;
        for (unsigned j(0); j<m; j++)
            fichier << (*Population[i])[G.edgeFromId(j)] << ";";
    }
    fichier.close();
}



template <typename VERTEXLIST, typename GRAPH_STRUCT>
unsigned nodeFromMyId(VERTEXLIST vertices, GRAPH_STRUCT g_struct, int index)
{
    for(unsigned i=0; i < vertices.size(); i++)
    {
        if(g_struct.NIM[vertices[i]]==index)
            return(i);
    }
    return(-1);
}


template <>
void transform_graph(lemonGraph const &OG, lemonDigraph *NG, GRAPH_STR g_str, DIGRAPH_STR *dg_str, int const s_i, int const t_i)
{

    dg_str->s = NG->addNode();
    dg_str->DNIM[dg_str->s] = g_str.NIM[g_str.LNodes[s_i]];
    dg_str->linking[g_str.NIM[g_str.LNodes[s_i]]] = 0;

    dg_str->t = NG->addNode();
    dg_str->DNIM[dg_str->t] = g_str.NIM[g_str.LNodes[t_i]];
    dg_str->linking[g_str.NIM[g_str.LNodes[t_i]]] = 1;

    vector<lemonNode> LNodes_st(g_str.LNodes);

    //cout << "erasing the " << find_vertex(LNodes_st, g_str, s_i) << " vertex" << endl;
    LNodes_st.erase(LNodes_st.begin() + nodeFromMyId(LNodes_st, g_str, s_i));
    //cout << "erasing the " << find_vertex(LNodes_st, g_str, t_i) << " vertex" << endl;
    LNodes_st.erase(LNodes_st.begin() + nodeFromMyId(LNodes_st, g_str, t_i));


    for (unsigned i = 0; i < LNodes_st.size(); i++)
    {
        dg_str->LDNodes1.push_back(NG->addNode());
        dg_str->DNIM[dg_str->LDNodes1[dg_str->LDNodes1.size()-1]] = g_str.NIM[LNodes_st[i]];
        dg_str->linking[g_str.NIM[LNodes_st[i]]] = (dg_str->LDNodes1.size()+dg_str->LDNodes2.size()+1);

    }

    for (unsigned i = 0; i < LNodes_st.size(); i++)
    {
        dg_str->LDNodes2.push_back(NG->addNode());
        dg_str->DNIM[dg_str->LDNodes2[dg_str->LDNodes2.size()-1]] = g_str.NIM[LNodes_st[i]] + 1000000;
        dg_str->linking[g_str.NIM[LNodes_st[i]] + 1000000] = dg_str->LDNodes1.size() + dg_str->LDNodes2.size() +1;

    }

    //for(map<int,int>::const_iterator it = dg_str->linking.begin(); it != dg_str->linking.end(); ++it)
        //cout << it->first << " " << it->second << endl;

    for (unsigned i = 0; i < g_str.LEdges.size(); i++)
    {
        int source_i, target_i;
        source_i = g_str.NIM[OG.u(g_str.LEdges[i])];
        target_i = g_str.NIM[OG.v(g_str.LEdges[i])];

        if(source_i == s_i) // case 1: arc sX
        {
            if(target_i == t_i) // case 1.1: arc st
            {
                int new_pos_x = dg_str->linking[s_i];
                int new_pos_y = dg_str->linking[t_i];
                dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x), NG->nodeFromId(new_pos_y)));

                dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
                dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
                dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];
            }
            else // case 1.2: arc su
            {
                int new_pos_x = dg_str->linking[s_i];
                int new_pos_y = dg_str->linking[target_i];
                dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

                dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
                dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
                dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];
            }
        }
        else if(target_i == s_i) // case 2: arc Xs
        {
            if(source_i == t_i) // case 2.1: arc ts
            {
                int new_pos_x = dg_str->linking[s_i];
                int new_pos_y = dg_str->linking[t_i];
                dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

                dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
                dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
                dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];
            }
            else // case 2.2: arc us
            {
                int new_pos_x = dg_str->linking[s_i];
                int new_pos_y = dg_str->linking[source_i];
                dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

                dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
                dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
                dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];
            }
        }
        else if(target_i == t_i) // case 3: arc ut
        {
            int new_pos_x = dg_str->linking[source_i + 1000000];
            int new_pos_y = dg_str->linking[t_i];
            dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

            dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
            dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
            dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];
        }
        else if(source_i == t_i) // case 4: arc tu
        {
            int new_pos_x = dg_str->linking[1000000 + target_i];
            int new_pos_y = dg_str->linking[t_i];
            dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

            dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
            dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
            dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];
        }
        else // case 5: arc uv
        {
            int new_pos_x = dg_str->linking[source_i];
            int new_pos_y = dg_str->linking[1000000 + target_i];
            dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

            dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
            dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
            dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];


            new_pos_x = dg_str->linking[target_i];
            new_pos_y = dg_str->linking[1000000 + source_i];
            dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

            dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EIM[g_str.LEdges[i]];
            dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.ECM[g_str.LEdges[i]];
            dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = g_str.EWM[g_str.LEdges[i]];
        }
    }

        /// THIS CAN BE DONE UP!!!! (TO VERIFY CAPACITY AND WEIGHT IN ARCS?!)
        for (unsigned i = 0; i < LNodes_st.size(); i++)
        {
            int new_pos_x = dg_str->linking[g_str.NIM[LNodes_st[i]]];
            int new_pos_y = dg_str->linking[g_str.NIM[LNodes_st[i]]+1000000];
            dg_str->LArcs.push_back(NG->addArc(NG->nodeFromId(new_pos_x),NG->nodeFromId(new_pos_y)));

            dg_str->AIM[dg_str->LArcs[dg_str->LArcs.size()-1]] = 276447231;
            dg_str->ACM[dg_str->LArcs[dg_str->LArcs.size()-1]] = 1;
            dg_str->AWM[dg_str->LArcs[dg_str->LArcs.size()-1]] = 0;
        }

        dg_str->n = countNodes(*NG);
        dg_str->m = countArcs(*NG);

        //write_graph(*NG, *dg_str);
}



template <>
bool is_feasable(lemonGraph const &G, vector < lemonDigraph* > const &DGS, GRAPH_STR const &g_str, vector < DIGRAPH_STR* > const &dg_strs, vector<lemonArcValMap*> const &VRAVM, lemonEdgeValMap const &REVM)
{
    for(unsigned i=0; i < dg_strs.size(); i++)
    {
        for(unsigned j=0; j < dg_strs[i]->m; j++)
        {
            lemonArc a = DGS[i]->arcFromId(j);
            if (dg_strs[i]->AIM[a] < 276447231)
            {
                lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);
                if ((*VRAVM[i])[a] > REVM[e])
                    return (false);
            }
        }
    }
    return(true);
}



template <>
void MyShortPaths(lemonDigraph const &DG, DIGRAPH_STR *dg_str)
{
    dg_str->DNVM[dg_str->s] = 0;
    dg_str->DNVM[dg_str->t] = dg_str->AWM[findArc(DG,dg_str->s,dg_str->t)];

    for (unsigned j = 0; j < dg_str->LDNodes1.size(); j++)
        dg_str->DNVM[dg_str->LDNodes1[j]] = dg_str->AWM[findArc(DG, dg_str->s, dg_str->LDNodes1[j])];

    for (unsigned j = 0; j < dg_str->LDNodes2.size(); j++)
    {
        dg_str->DNVM[dg_str->LDNodes2[j]] = 9999999;
        for (ListDigraph::InArcIt e(DG, dg_str->LDNodes2[j]); e != INVALID; ++e)
            if(dg_str->DNVM[DG.source(e)]+dg_str->AWM[e] < dg_str->DNVM[dg_str->LDNodes2[j]])
                dg_str->DNVM[dg_str->LDNodes2[j]] =  dg_str->DNVM[DG.source(e)] + dg_str->AWM[e];
    }

    // Adjusting the shortest path value of t
    for (ListDigraph::InArcIt e(DG, dg_str->t); e != INVALID; ++e)
        if(dg_str->DNVM[DG.source(e)]+dg_str->AWM[e] < dg_str->DNVM[dg_str->t])
            dg_str->DNVM[dg_str->t] =  dg_str->DNVM[DG.source(e)] + dg_str->AWM[e];
}



template <>
void MyMinCostFlow(lemonDigraph const &DG, DIGRAPH_STR *dg_str)
{
    /*
    cout << "Calculating the shortest paths from s : Done ... OK!" << endl;

    for (unsigned j = 0; j < dg_str.LDNodes1.size(); j++)
    {
       lemonArc su = findArc(G, dg_str.s, dg_str.LDNodes1[j]);
       for (ListDigraph::OutArcIt e(G, dg_str.LDNodes1[j]); e != INVALID; ++e)
           dg_str.AVM[e] =
    }
    */
}


void test_HA1(vector < vector <int> > const &demands, string const name)
{
    lemonGraph G;

    vector<lemonNode> LNodes;
    vector<lemonEdge> LEdges;

    lemonNodeIndexMap NIM(G,0);	//Pour numéroter les noeuds
    lemonNodeDegreeMap NDM(G,0);
    lemonNodeValMap NVM(G,0);

    lemonEdgeIndexMap EIM(G,0);	//Pour numéroter les arêtes
    lemonEdgeCapacityMap ECM(G,0); //Pour donner une capacité pour chaque arête
    lemonEdgeWeightMap EWM(G,0);
    lemonEdgeValMap EVM(G,0);

    GRAPH_STR g_str_sol = {name,LNodes,LEdges,NIM,NDM,NVM,EIM,ECM,EWM,EVM};

    open_SOL_graph(G,g_str_sol,name);
}



void Selection(unsigned const l,unsigned &p1, unsigned &p2, unsigned &alchimy, vector <string> &pairs)
{
    string temp1,temp2;

    do
    {
        if(l > 5)
        {
            unsigned l1,l2,l3,l4;

            l1 = l / 5;
            l2 = ((double)l / 5)*2;
            l3 = ((double)l / 5)*3;
            l4 = ((double)l / 5)*4;

            alchimy = 0;

            /* generate secret number between 0 and 100: */
            unsigned val = rand() % 100;

            if (val <= 67)          //  67%
            {
                //p1 = rand() % l1 + 0;
                p1 = 0 + rand() / (RAND_MAX / (l1 - 0 + 1) + 1);
                alchimy += 50;
            }
            else if (val <= 86)     //  19%
            {
                //p1 = rand() % l2 + l1;
                p1 = l1 + rand() / (RAND_MAX / (l2 - l1 + 1) + 1);
                alchimy += 40;

            }
            else if (val <= 96)     //  10%
            {
                //p1 = rand() % l3 + l2;
                p1 = l2 + rand() / (RAND_MAX / (l3 - l2 + 1) + 1);
                alchimy += 30;
            }
            else if (val <= 99)     //  3%
            {
                //p1 = rand() % l4 + l3;
                p1 = l3 + rand() / (RAND_MAX / (l4 - l3 + 1) + 1);
                alchimy += 20;
            }
            else                    //  1%
            {
                //p1 = rand() % l + l4;
                p1 = l4 + rand() / (RAND_MAX / (l - l4 + 1) + 1);
                alchimy += 10;
            }

            val = rand() % 100;

            if (val <= 67)          //  67%
            {
                do
                    //p2 = rand() % l1 + 0;
                    p2 = 0 + rand() / (RAND_MAX / (l1 - 0 + 1) + 1);
                while (p2==p1);
                alchimy += 50;
            }
            else if (val <= 86)     //  19%
            {
                do
                    //p2 = rand() % l2 + l1;
                    p2 = l1 + rand() / (RAND_MAX / (l2 - l1 + 1) + 1);
                while (p2==p1);
                alchimy += 40;
            }
            else if (val <= 96)     //  10%
            {
                do
                    //p2 = rand() % l3 + l2;
                    p2 = l2 + rand() / (RAND_MAX / (l3 - l2 + 1) + 1);
                while (p2==p1);
                alchimy += 30;
            }
            else if (val <= 99)     //  3%
            {
                do
                    //p2 = rand() % l4 + l3;
                    p2 = l3 + rand() / (RAND_MAX / (l4 - l3 + 1) + 1);
                while (p2==p1);
                alchimy += 20;
            }
            else                    //  1%
            {
                do
                    //p2 = rand() % l + l4;
                    p2 = l4 + rand() / (RAND_MAX / (l - l4 + 1) + 1);
                while (p2==p1);
                alchimy += 10;
            }
        }
        else
        {
            //p1 = rand() % l + 0;
            p1 = 0 + rand() / (RAND_MAX / (l - 0 + 1) + 1);
            do
                //p2 = rand() % l + 0;
                p2 = 0 + rand() / (RAND_MAX / (l - 0 + 1) + 1);
            while (p2==p1);
            alchimy = 100;
        }


        temp1 =  to_string(p1) + to_string(p2);
        temp2 =  to_string(p2) + to_string(p1);
    }
    while ((find(pairs.begin(), pairs.end(), temp1) != pairs.end() || find(pairs.begin(), pairs.end(), temp2) != pairs.end()) && l > 4);

    pairs.push_back(temp1);
}



vector <pair <unsigned, unsigned> > Crossover(unsigned const l, unsigned const alchimy)
{
    vector <pair <unsigned, unsigned> > v;
    vector <string> pairs;
    unsigned val1, val2;
    string temp1;

    unsigned nb = ((double)l/300)*alchimy;

    //cout << " nb = " << nb << " l = " << l << endl;

    do
    {
        do
        {
            val1 = 0 + rand() / (RAND_MAX / (l - 1 - 0 + 1) + 1);
            val2 = val1 + rand() / (RAND_MAX / (l - val1 + 1) + 1);

            temp1 =  to_string(val1) + to_string(val2);
        }
        while (find(pairs.begin(), pairs.end(), temp1) != pairs.end());

        //cout << "val1 = " << val1 << "   val 2 = " << val2 << endl;

        pair <unsigned, unsigned> p(val1,val2);
        v.push_back(p);
        pairs.push_back(temp1);
    }
    while (v.size() < nb);

    return(v);
}

vector <pair <unsigned, unsigned> > CrossoverFull(unsigned const l)
{
    vector <pair <unsigned, unsigned> > v;

    //cout << " nb = " << nb << " l = " << l << endl;
    for (unsigned i=0; i <l; i++)
    {
        for (unsigned j=i+2; j <l; j++)
        {
            pair <unsigned, unsigned> p(i,j);
            v.push_back(p);
        }
    }

    return(v);
}
