#include <sys/time.h>
#include <iostream>

#include <thread>
#include <mutex>

#include <omp.h>

#include "include/functions.h"
#include "include/graph.h"

#include <lemon/lgf_writer.h>
#include <lemon/preflow.h>
#include <lemon/network_simplex.h>





using namespace std;
using namespace lemon;





//door handle
mutex printLock, VELock, ZLock;

bool Lagflag (true), Succflag(true);

vector < lemonEdgeValMap* >  Population;
vector < int > PopRanking;
vector < vector < lemonArcValMap* > >  PopOrgans;

unsigned LB;
unsigned UB(0);




void fusionerPop(const int debut1,const unsigned fin1,const unsigned fin2, const char c)
{
    vector < lemonEdgeValMap* > P;
    vector < int > PR;
    vector < vector < lemonArcValMap* > >  POs;

    unsigned debut2 = fin1+1;
    unsigned compteur1 = debut1;
    unsigned compteur2 = debut2;


    // copie des éléments du début de tableau
    for(unsigned i = debut1; i<=fin1; i++)
    {
        P.push_back(Population[i]);
        PR.push_back(PopRanking[i]);
    }

    for(unsigned j = 0; j<PopOrgans.size(); j++)
    {
        POs.push_back(vector < lemonArcValMap* > ());
        for(unsigned i = debut1; i<=fin1; i++)
            POs[j].push_back(PopOrgans[j][i]);
    }


    // fusion des deux tableaux
    for(unsigned i=debut1; i<=fin2; i++)
    {
        if(compteur1==debut2) // éléments du 1er tableau tous utilisés
            break; // éléments tous classés

        else if(compteur2==(fin2+1)) // éléments du 2nd tableau tous utilisés
        {
            // copie en fin de tableau des éléments du 1er sous tableau
            Population[i] = P[compteur1-debut1];
            PopRanking[i] = PR[compteur1-debut1];
            for(unsigned j = 0; j<PopOrgans.size(); j++)
                PopOrgans[j][i] = POs[j][compteur1-debut1];

            compteur1++;
        }

        else if( PR[compteur1-debut1] >= PopRanking[compteur2] && c == '>')                                  //tableau2[compteur1-debut1]<tableau[compteur2])
        {
            // ajout d'1 élément du 1er sous tableau
            Population[i] = P[compteur1-debut1];
            PopRanking[i] = PR[compteur1-debut1];
            for(unsigned j = 0; j<PopOrgans.size(); j++)
                PopOrgans[j][i] = POs[j][compteur1-debut1];

            compteur1++; // on avance ds le 1er sous tableau
        }

        else if( PR[compteur1-debut1] <= PopRanking[compteur2] && c == '<')                                  //tableau2[compteur1-debut1]<tableau[compteur2])
        {
            // ajout d'1 élément du 1er sous tableau
            Population[i] = P[compteur1-debut1];
            PopRanking[i] = PR[compteur1-debut1];
            for(unsigned j = 0; j<PopOrgans.size(); j++)
                PopOrgans[j][i] = POs[j][compteur1-debut1];

            compteur1++; // on avance ds le 1er sous tableau
        }

        else
        {
            // copie de l'élément à la suite du tableau
            Population[i] = Population[compteur2];
            PopRanking[i] = PopRanking[compteur2];
            for(unsigned j = 0; j<PopOrgans.size(); j++)
                PopOrgans[j][i] = PopOrgans[j][compteur2];

            compteur2++; // on avance ds le 2nd sous tableau
        }
    }
    //free(tableau2);
}


void triFusionAuxPop(const unsigned debut, const unsigned fin, const char c)
{
    if(debut!=fin) // condition d'arrêt
    {
        int milieu = (debut+fin)/2;

        #pragma omp parallel sections
        {
            #pragma omp section
            triFusionAuxPop(debut, milieu, c); // trie partie1
            #pragma omp section
            triFusionAuxPop(milieu+1, fin, c); // trie partie2
        }
        fusionerPop(debut, milieu, fin, c); // fusion des deux parties
    }
}


void triFusionPop(const char c)
{
    if(Population.size() > 0)
        triFusionAuxPop(0, Population.size()-1, c);
}





void Lagrangien_Thread(lemonGraph const &G, vector < lemonDigraph* > const &DGS, GRAPH_STR const &g_str, vector < DIGRAPH_STR* > const &dg_strs, vector < vector <int> > const &demands, int const k)
{
    unsigned m(countEdges(G));

    vector < vector <double> > lambda;


    vector<lemonArcValMap*> VRAVM;
    vector<lemonArcWeightMap*> VRAWM;

    lemonEdgeValMap *REVM;


    for (unsigned i=0; i < demands.size(); i++)
    {

        VRAVM.push_back(new lemonArcValMap(*DGS[i],0.0));
        VRAWM.push_back(new lemonArcWeightMap(*DGS[i],0));

        lambda.push_back(vector <double> (dg_strs[i]->m,0.0));
    }


    #pragma omp parallel for
    for (unsigned i=0; i < demands.size(); i++)
    {
        //Calculating minimum cost max flow using NetworkSimplex from LEMON
        NetworkSimplex<lemonDigraph, int> myNetworkSimplex(*DGS[i]);
        myNetworkSimplex.upperMap(dg_strs[i]->ACM).costMap(dg_strs[i]->AWM);
        myNetworkSimplex.stSupply(dg_strs[i]->s, dg_strs[i]->t, k);
        myNetworkSimplex.run();


        myNetworkSimplex.flowMap(*VRAVM[i]);
    }

    //cout << endl << "Starting the sub-gradient phase : START ..." << endl;

    //1- donner une valeur initiale à lambda
    //2- calculer la fonction à minimiser
    //3- résoudre la minimisation pour x(e)
    //4- réitérer en calculant la nouvelle valeur de lambda en fonction de l ancienne et en voyant si on améliore x(e) ou pas

    lemonEdgeLambdaMap ELM(G,0.0);

    //HERE CALCULATING THE OBJECTIVE FORMULA
    double obj(0.0), objOpt(0.0), sum_L2(0.0), norme_2(0.0), alpha(2.0);
    vector < vector <int> > derive;
    unsigned Z(0);
    unsigned old_Z;


    vector<lemonEdge> LEdges;

    for(unsigned i=0; i < m; i++)
        LEdges.push_back(G.edgeFromId(i));

    REVM = new lemonEdgeValMap(G,0);

    for(unsigned i=0; i < dg_strs.size(); i++)
    {
        for(unsigned j=0; j < dg_strs[i]->m; j++)
        {
            lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated

            if(dg_strs[i]->AIM[a] < 276447231)
            {
                lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);
                lambda[i][j] = (double)g_str.EWM[e]/demands.size();      /// Important!

                ELM[e] += lambda[i][j];

                if((*REVM)[e] < 1 && (*VRAVM[i])[a] == 1)
                    (*REVM)[e] = 1;
            }
            else
                lambda[i][j] = 0.0;

            sum_L2 += lambda[i][j] * (*VRAVM[i])[a];       // value of arcs of the form uu' are equal to zero
        }
    }

    for (unsigned i=0; i < m; i++)
    {
        // Calculating feasable solution and Zopt
        if ((*REVM)[G.edgeFromId(i)] == 1)
            Z += g_str.EWM[G.edgeFromId(i)];
    }

    bool feasable = is_feasable(G,DGS,g_str,dg_strs,VRAVM,*REVM);

    if (feasable)
    {
        VELock.lock();
        Population.push_back(REVM);
        PopRanking.push_back(Z);
        for (unsigned i=0; i < demands.size(); i++)
            PopOrgans[i].push_back(VRAVM[i]);
        VELock.unlock();

        REVM = new lemonEdgeValMap(G,0);
    }

    obj += sum_L2;

    for(unsigned i=0; i < m; i++)
    {
        float sum_L1_temp = g_str.EWM[LEdges[i]] - ELM[LEdges[i]];
        if(sum_L1_temp <= 0.0)
        {
            //cout << "updating xe because the sum_L1 is equal to = " << sum_L1_temp << endl;
            (*REVM)[G.edgeFromId(i)] = 1;
            obj += sum_L1_temp;     // That means that x(e) = 1
        }
        else
            (*REVM)[G.edgeFromId(i)] = 0;
    }


    for(unsigned i=0; i < dg_strs.size(); i++)
    {
        derive.push_back(vector <int>());
        for(unsigned j=0; j < dg_strs[i]->m; j++)
        {
            lemonArc a = DGS[i]->arcFromId(j);
            if(dg_strs[i]->AIM[a] < 276447231)
            {
                lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);
                derive[i].push_back((*REVM)[e] - (*VRAVM[i])[a]);   //  = x(e) - sum fi(a)
            }
            else
                derive[i].push_back(0);

            norme_2 += derive[i][j] * derive[i][j];
        }
    }

    //write_graph(G,g_str);

    ZLock.lock();
    if (Z < UB || UB==0)
        UB = Z;
    ZLock.unlock();

    /*printLock.lock();
    if(feasable)
        cout << "Z = " << Z << " and the solution is feasable!" << endl;
    else
        cout << "Z = " << Z << " but the solution is not feasable! :/" << endl;
    //cout << "Z = " << Z << endl;
    cout << "L(lambda,x) = " << obj << endl;
    cout << "Zsup = " << Zsup << endl;
    cout << "Z* = " << Zopt << endl << endl;
    printLock.unlock();*/


    unsigned stabilized(0), stabilized2(0), nb_iteration(2000);

    for(unsigned it(1); it <= nb_iteration && stabilized <= 20; it++)
    {
        if (stabilized2 == 20)
        {
            alpha /= 2;
            stabilized2 = 0;

            /* for(unsigned i=0; i < dg_strs.size(); i++)
                for(unsigned j=0; j < dg_strs[i]->m; j++)
                    lambda[i][j] = 0; */

        }

        double st;
        if (norme_2)
            st = (float) alpha * (UB - obj) / norme_2;
        else
            st = 0;

        // Updating lambda's value
        for(unsigned i=0; i < dg_strs.size(); i++)
        {
            for(unsigned j=0; j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j];   // Arcs have changed from the transformation and LArcs is not dynamically updated

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lambda[i][j] -= st * derive[i][j];
                    if (lambda[i][j] < 0)
                        lambda[i][j] = 0;
                }
            }
        }


        for (unsigned i=0; i < dg_strs.size(); i++)
        {
            for (unsigned j=0; j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);
                (*VRAWM[i])[a] = lambda[i][j];
            }

            VRAVM[i] = new lemonArcValMap(*DGS[i]);
        }


        #pragma omp parallel for
        for (unsigned i=0; i < dg_strs.size(); i++)
        {

            //Calculating minimum cost max flow using NetworkSimplex from LEMON
            NetworkSimplex<lemonDigraph, int> myNetworkSimplex(*DGS[i]);
            myNetworkSimplex.upperMap(dg_strs[i]->ACM).costMap(*VRAWM[i]);
            myNetworkSimplex.stSupply(dg_strs[i]->s, dg_strs[i]->t, k);
            myNetworkSimplex.run();


            myNetworkSimplex.flowMap(*VRAVM[i]);
        }


        // Reinitializing structures
        old_Z = Z;
        Z = 0;
        sum_L2 = obj = norme_2 = 0.0;
        for(unsigned i=0; i < m; i++)
            ELM[G.edgeFromId(i)] = 0;

        // Recalculating the second part of the objective function
        for(unsigned i=0; i < dg_strs.size(); i++)
        {
            for(unsigned j=0; j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j];   // Arcs have changed from the transformation and LArcs is not dynamically updated

                sum_L2 += lambda[i][j] * (*VRAVM[i])[a];

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);
                    ELM[e] += lambda[i][j];
                }
            }
        }

        obj += sum_L2;

        // Recalculating the first part of the objective function and xt(e)
        for(unsigned i=0; i < m; i++)
        {
            float sum_L1_temp = g_str.EWM[G.edgeFromId(i)] - ELM[G.edgeFromId(i)];

            if(sum_L1_temp <= 0.0)
            {
                (*REVM)[G.edgeFromId(i)] = 1;
                obj += sum_L1_temp;     // That means that x(e) = 1
                Z += g_str.EWM[G.edgeFromId(i)];
            }
            else
                (*REVM)[G.edgeFromId(i)] = 0;
        }

        // Recalculate x(e) - f(a)
        for(unsigned i=0; i < dg_strs.size(); i++)
        {
            for(unsigned j=0; j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);
                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);
                    derive[i][j] = (*REVM)[e] - (*VRAVM[i])[a];
                    norme_2 += derive[i][j] * derive[i][j];
                }
            }
        }

        if (old_Z == Z)
            stabilized++;
        else
            stabilized = 0;

        if (obj > objOpt)
        {
            objOpt = obj;
            stabilized2 = 0;
        }
        else
            stabilized2++;

        feasable = is_feasable(G,DGS,g_str,dg_strs,VRAVM,(*REVM));

        if (feasable)
        {
            VELock.lock();
            Population.push_back(REVM);
            PopRanking.push_back(Z);
            for (unsigned i=0; i < demands.size(); i++)
                PopOrgans[i].push_back(VRAVM[i]);
            VELock.unlock();

            REVM = new lemonEdgeValMap(G,0);

            ZLock.lock();
            if (Z < UB || UB==0)
                UB = Z;
            ZLock.unlock();
        }
        else
        {
            // Primal Heuristic to make the solution feasable for the initial P
            for(unsigned i=0; i < dg_strs.size(); i++)
            {
                for(unsigned j=0; j < dg_strs[i]->m; j++)
                {
                    lemonArc a = DGS[i]->arcFromId(j);
                    if (dg_strs[i]->AIM[a] < 276447231)
                    {
                        lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);
                        if ((*VRAVM[i])[a] > (*REVM)[e])
                            (*REVM)[e] = 1;
                    }
                }
            }

            bool feasable2(is_feasable(G,DGS,g_str,dg_strs,VRAVM,*REVM));

            unsigned Zf = 0;
            for(unsigned i=0; i < m; i++)
                Zf += g_str.EWM[G.edgeFromId(i)] * (*REVM)[G.edgeFromId(i)];

            if (feasable2)
            {
                VELock.lock();
                Population.push_back(REVM);
                PopRanking.push_back(Zf);
                for (unsigned i=0; i < demands.size(); i++)
                    PopOrgans[i].push_back(VRAVM[i]);
                VELock.unlock();

                REVM = new lemonEdgeValMap(G,0);

                ZLock.lock();
                if (Zf < UB || UB==0)
                    UB = Zf;
                ZLock.unlock();
            }
        }

        if((long)obj == UB)
            stabilized = 100;

        //cout << it+1 << "'th iteration : alpha =" << alpha << endl;
    }

    Lagflag = false;

    LB = objOpt;

    //cout << "L(lambda,x)* = " << objOpt << endl << endl;
}





void Genetic_Thread(lemonGraph const &G, vector < lemonDigraph* > const &DGS, GRAPH_STR const &g_str, vector < DIGRAPH_STR* > const &dg_strs)
{
    // ProcessTime example
    struct timeval startTime;
    struct timeval endTime;
    double tS, tE;
    // get the current time
    // - NULL because we don't care about time zone
    gettimeofday(&startTime, NULL);

    unsigned childn(0), m(countEdges(G));

    vector<lemonArcValMap*> VRAVM;

    lemonEdgeValMap *REVM;

    unsigned Z;

    unsigned generationLimit = 0;

    while (true)
    {
        if ((Lagflag || Succflag) && Population.size() > 1)
        {
            // Evaluation is made here by sorting the population according to the objective function value z
            VELock.lock();
            triFusionPop('<');
            VELock.unlock();

            //cout << " Best = " << PopRanking[0] << endl;
            //cout << " Worst " << Population.size() << " = " << PopRanking[Population.size()] << endl;

            ZLock.lock();
            if (PopRanking[0] < UB || UB==0)
                UB = PopRanking[0];
            ZLock.unlock();

            unsigned ela;
            unsigned nbc;

            if(generationLimit > 100) {
                //cout << " childn = " << childn << endl;
                ela = childn;
                nbc = generationLimit/10;
            }
            else if (generationLimit > 10)
            {
                //cout << " childn = " << childn << endl;
                ela = childn;
                nbc = generationLimit/2;
            }
            else
            {
                ela = 0;
                nbc = generationLimit/2;
            }

            childn = 0;

            VELock.lock();
            Population.erase(Population.end()-ela,Population.end());
            PopRanking.erase(PopRanking.end()-ela,PopRanking.end());
            for (unsigned j=0; j < PopOrgans.size(); j++)
                PopOrgans[j].erase(PopOrgans[j].end()-ela,PopOrgans[j].end());
            VELock.unlock();

            cout << " Taille = " << Population.size() << "    Best = " << UB << endl;

            //cout << " Worst " << Population.size() << " = " << PopRanking[Population.size()] << endl;

            // Selection is made by selecting parents randomly from population but with probability distribution (from ) as follow 67%, 19%, 10%, 3% and 1%

            generationLimit = Population.size()-1;

            vector <string> pairs;

            #pragma omp parallel for private(pairs, Z, REVM, VRAVM)
            for (unsigned itc=1; itc <= nbc; itc++)
            {
                unsigned p1, p2, alchimy;

                Selection(generationLimit, p1, p2, alchimy, pairs);

                //cout << "p1 = " << p1 << " p2 = " << p2 << endl;

                // Crossover
                vector <pair <unsigned, unsigned> > cuttingPoints;
                if(generationLimit > 5)
                    cuttingPoints = Crossover(PopOrgans.size()-1 , alchimy);
                else
                    cuttingPoints = Crossover(PopOrgans.size()-1, alchimy);

                // Reproduction
                for(unsigned k=0; k < cuttingPoints.size(); k++)
                {
                    Z = 0;

                    REVM = new lemonEdgeValMap(G,0);

                    for (unsigned i=0; i < PopOrgans.size(); i++)
                    {
                        VRAVM.push_back(new lemonArcValMap(*DGS[i],0.0));

                        for(unsigned j=0; j < dg_strs[i]->m; j++)
                        {
                            lemonArc a = DGS[i]->arcFromId(j);

                            if (i < cuttingPoints[k].first)
                                (*VRAVM[i])[a] = (*PopOrgans[i][p1])[a];
                            else if ( i < cuttingPoints[k].second || cuttingPoints[k].first == cuttingPoints[k].second)
                                (*VRAVM[i])[a] = (*PopOrgans[i][p2])[a];
                            else
                                (*VRAVM[i])[a] = (*PopOrgans[i][p1])[a];


                            if(dg_strs[i]->AIM[a] < 276447231)
                            {
                                lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                                if((*REVM)[e] < 1 && (*VRAVM[i])[a] == 1)
                                    (*REVM)[e] = 1;
                            }

                        }
                    }

                    for (unsigned j=0; j < m; j++)
                        if ((*REVM)[G.edgeFromId(j)] == 1)
                            Z += g_str.EWM[G.edgeFromId(j)];

                    //cout << "Z = " << Z << endl;

                    VELock.lock();
                    Population.push_back(REVM);
                    PopRanking.push_back(Z);
                    for (unsigned j=0; j < PopOrgans.size(); j++)
                        PopOrgans[j].push_back(VRAVM[j]);
                    VELock.unlock();

                    ZLock.lock();
                    if (Z < UB || UB==1)
                        UB = Z;
                    ZLock.unlock();

                    childn++;
                }
            }
        }
        else if (Population.size()>0 && !Lagflag && !Succflag)
        {
            // Evaluation is made here by sorting the population according to the objective function value z
            VELock.lock();
            triFusionPop('<');
            VELock.unlock();

            //cout << " Best = " << PopRanking[0] << endl;
            //cout << " Worst " << Population.size() << " = " << PopRanking[Population.size()] << endl;

            ZLock.lock();
            if (PopRanking[0] < UB || UB==0)
                UB = PopRanking[0];
            ZLock.unlock();

            unsigned ela;
            unsigned nbc;

            if(generationLimit > 100)
            {
                //cout << " childn = " << childn << endl;
                ela = childn + Population.size()/20;
                nbc = generationLimit/10;
            }
            else
            {
                //cout << " childn = " << childn << endl;
                ela = childn + 1;
                nbc = generationLimit/2;
            }

            cout << "childn = " << childn << endl;

            childn = 0;

            VELock.lock();
            Population.erase(Population.end()-ela,Population.end());
            PopRanking.erase(PopRanking.end()-ela,PopRanking.end());
            for (unsigned j=0; j < PopOrgans.size(); j++)
                PopOrgans[j].erase(PopOrgans[j].end()-ela,PopOrgans[j].end());
            VELock.unlock();

            cout << "ela = " << ela << endl ;

            cout << " 2 - Taille = " << Population.size() << "    Best = " << UB << endl;

            //cout << " Worst " << Population.size() << " = " << PopRanking[Population.size()] << endl;

            if (Population.size()==1 && !Lagflag && !Succflag)
                break;

            // Selection is made by selecting parents randomly from population but with probability distribution (from ) as follow 67%, 19%, 10%, 3% and 1%

            generationLimit = Population.size()-1;

            vector <string> pairs;

            #pragma omp parallel for private(pairs, Z, REVM, VRAVM)
            for (unsigned itc=1; itc <= nbc; itc++)
            {
                unsigned p1, p2, alchimy;

                Selection(generationLimit, p1, p2, alchimy, pairs);

                //cout << "p1 = " << p1 << " p2 = " << p2 << endl;

                // Crossover
                vector <pair <unsigned, unsigned> > cuttingPoints;
                if(generationLimit > 5)
                    cuttingPoints = Crossover(PopOrgans.size()-1 , alchimy);
                else
                    cuttingPoints = Crossover(PopOrgans.size()-1, alchimy);

                // Reproduction
                for(unsigned k=0; k < cuttingPoints.size(); k++)
                {
                    Z = 0;

                    REVM = new lemonEdgeValMap(G,0);

                    for (unsigned i=0; i < PopOrgans.size(); i++)
                    {
                        VRAVM.push_back(new lemonArcValMap(*DGS[i],0.0));

                        for(unsigned j=0; j < dg_strs[i]->m; j++)
                        {
                            lemonArc a = DGS[i]->arcFromId(j);

                            if (i < cuttingPoints[k].first)
                                (*VRAVM[i])[a] = (*PopOrgans[i][p1])[a];
                            else if ( i < cuttingPoints[k].second || cuttingPoints[k].first == cuttingPoints[k].second)
                                (*VRAVM[i])[a] = (*PopOrgans[i][p2])[a];
                            else
                                (*VRAVM[i])[a] = (*PopOrgans[i][p1])[a];


                            if(dg_strs[i]->AIM[a] < 276447231)
                            {
                                lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                                if((*REVM)[e] < 1 && (*VRAVM[i])[a] == 1)
                                    (*REVM)[e] = 1;
                            }

                        }
                    }

                    for (unsigned j=0; j < m; j++)
                        if ((*REVM)[G.edgeFromId(j)] == 1)
                            Z += g_str.EWM[G.edgeFromId(j)];

                    //cout << "Z = " << Z << endl;

                    VELock.lock();
                    Population.push_back(REVM);
                    PopRanking.push_back(Z);
                    for (unsigned j=0; j < PopOrgans.size(); j++)
                        PopOrgans[j].push_back(VRAVM[j]);
                    VELock.unlock();

                    ZLock.lock();
                    if (Z < UB || UB==1)
                        UB = Z;
                    ZLock.unlock();

                    childn++;
                }
            }

            gettimeofday(&endTime, NULL);
            // calculate time in microseconds
            tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
            tE = endTime.tv_sec*1000000  + (endTime.tv_usec);

            if (((tE - tS)/ 1000000.0L) >= 7200)
                break;
        }
    }
}





void Successive_Thread_Random(lemonGraph const &G, vector < lemonDigraph* > const &DGS, GRAPH_STR const &g_str, vector < DIGRAPH_STR* > const &dg_strs, vector < vector <int> > const &demands, int const k)
{
    unsigned m(countEdges(G));

    #pragma omp parallel for
    for (unsigned it=0; it < demands.size(); it++)
    {
        vector<lemonArcValMap*> VRAVM(demands.size());
        vector<lemonArcWeightMap*> VRAWM(demands.size());

        lemonEdgeValMap *REVM;

        vector < vector <int> > CopyD;

        for (unsigned i(0); i < demands.size(); i++)
            CopyD.push_back(vector <int> (demands[i]));
        if(it > 0)
            random_shuffle( CopyD.begin(), CopyD.end() );

        REVM = new lemonEdgeValMap(G,0);

        /// setting in 1 all fixed edges
        for (unsigned i(0); i < CopyD.size(); i++)
            (*REVM)[findEdge(G, G.nodeFromId(CopyD[i][0]), G.nodeFromId(CopyD[i][1]))] = 1;

        for (unsigned iv(0); iv < CopyD.size(); iv++)
        {
            unsigned i = find(demands.begin(), demands.end(), CopyD[iv]) - demands.begin();

            VRAVM[i] = new lemonArcValMap(*DGS[i],0.0);
            VRAWM[i] = new lemonArcWeightMap(*DGS[i],0);

            // Adjusting costs!!!
            for (unsigned j(0); j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated
                (*VRAWM[i])[a] = dg_strs[i]->AWM[a];

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                    if((*REVM)[e] == 1)
                        (*VRAWM[i])[a] = 0;
                }
            }

            //Calculating minimum cost max flow using NetworkSimplex from LEMON
            NetworkSimplex<lemonDigraph, int> myNetworkSimplex(*DGS[i]);
            myNetworkSimplex.upperMap(dg_strs[i]->ACM).costMap(*VRAWM[i]);
            myNetworkSimplex.stSupply(dg_strs[i]->s, dg_strs[i]->t, k);
            myNetworkSimplex.run();


            myNetworkSimplex.flowMap(*VRAVM[i]);


            // updating superposing values in e
            for(unsigned j(0); j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                    if((*REVM)[e] < 1 && (*VRAVM[i])[a] == 1)
                        (*REVM)[e] = 1;
                }
            }
        }


        //HERE CALCULATING THE OBJECTIVE FORMULA
        unsigned Z(0);

        // Calculating realisable solution and Zopt
        for(unsigned i=0; i < m; i++)
            if ((*REVM)[G.edgeFromId(i)] == 1)
                Z += g_str.EWM[G.edgeFromId(i)];

        VELock.lock();
        Population.push_back(REVM);
        PopRanking.push_back(Z);
        for (unsigned i=0; i < CopyD.size(); i++)
            PopOrgans[i].push_back(VRAVM[i]);
        VELock.unlock();

        ZLock.lock();
        if (Z < UB || UB==1)
            UB = Z;
        ZLock.unlock();

        /*printLock.lock();
        cout << endl <<" Successive Heuristic!!!" << endl;
        if(feasable)
            cout << "the solution is feasable!" << endl;
        else
            cout << "but the solution is not feasable! :/" << endl;
        cout << "Z = " << Zopt << endl << endl;
        printLock.unlock();*/
    }

    Succflag = false;
}



void Successive_Thread(lemonGraph const &G, vector < lemonDigraph* > const &DGS, GRAPH_STR const &g_str, vector < DIGRAPH_STR* > const &dg_strs, vector < vector <int> > const &demands, int const k)
{
    Succflag = true;

    unsigned m(countEdges(G));

    vector<lemonArcValMap*> VRAVM(demands.size());
    vector<lemonArcWeightMap*> VRAWM(demands.size());

    lemonEdgeValMap *REVM;

    for (unsigned it(0); it < demands.size(); it++)
    {
        REVM = new lemonEdgeValMap(G,0);

        /// setting in 1 all fixed edges
        for (unsigned i(0); i < demands.size(); i++)
            (*REVM)[findEdge(G, G.nodeFromId(demands[i][0]), G.nodeFromId(demands[i][1]))] = 1;

        for (unsigned i(it); i < demands.size(); i++)
        {
            VRAVM[i] = new lemonArcValMap(*DGS[i],0.0);
            VRAWM[i] = new lemonArcWeightMap(*DGS[i],0);

            // Adjusting costs!!!
            for (unsigned j(0); j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated
                (*VRAWM[i])[a] = dg_strs[i]->AWM[a];

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                    if((*REVM)[e] == 1)
                        (*VRAWM[i])[a] = 0;
                }
            }

            //Calculating minimum cost max flow using NetworkSimplex from LEMON
            NetworkSimplex<lemonDigraph, int> myNetworkSimplex(*DGS[i]);
            myNetworkSimplex.upperMap(dg_strs[i]->ACM).costMap(*VRAWM[i]);
            myNetworkSimplex.stSupply(dg_strs[i]->s, dg_strs[i]->t, k);
            myNetworkSimplex.run();


            myNetworkSimplex.flowMap(*VRAVM[i]);


            // updating superposing values in e
            for(unsigned j(0); j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                    if((*REVM)[e] < 1 && (*VRAVM[i])[a] == 1)
                        (*REVM)[e] = 1;
                }
            }
        }


        for (unsigned i(0); i < it; i++)
        {
            VRAVM[i] = new lemonArcValMap(*DGS[i],0.0);
            VRAWM[i] = new lemonArcWeightMap(*DGS[i],0);

            // Adjusting costs!!!
            for (unsigned j(0); j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated
                (*VRAWM[i])[a] = dg_strs[i]->AWM[a];

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                    if((*REVM)[e] == 1)
                        (*VRAWM[i])[a] = 0;
                }
            }

            //Calculating minimum cost max flow using NetworkSimplex from LEMON
            NetworkSimplex<lemonDigraph, int> myNetworkSimplex(*DGS[i]);
            myNetworkSimplex.upperMap(dg_strs[i]->ACM).costMap(*VRAWM[i]);
            myNetworkSimplex.stSupply(dg_strs[i]->s, dg_strs[i]->t, k);
            myNetworkSimplex.run();


            myNetworkSimplex.flowMap(*VRAVM[i]);


            // updating superposing values in e
            for(unsigned j(0); j < dg_strs[i]->m; j++)
            {
                lemonArc a = DGS[i]->arcFromId(j);   // dg_strs[i]->LArcs[j]; // Arcs have changed from the transformation and LArcs is not dynamically updated

                if(dg_strs[i]->AIM[a] < 276447231)
                {
                    lemonEdge e = G.edgeFromId(dg_strs[i]->AIM[a]);

                    if((*REVM)[e] < 1 && (*VRAVM[i])[a] == 1)
                        (*REVM)[e] = 1;
                }
            }
        }


        //HERE CALCULATING THE OBJECTIVE FORMULA
        unsigned Z(0);

        // Calculating realisable solution and Zopt
        for(unsigned i=0; i < m; i++)
            if ((*REVM)[G.edgeFromId(i)] == 1)
                Z += g_str.EWM[G.edgeFromId(i)];

        VELock.lock();
        Population.push_back(REVM);
        PopRanking.push_back(Z);
        for (unsigned i=0; i < demands.size(); i++)
            PopOrgans[i].push_back(VRAVM[i]);
        VELock.unlock();

        ZLock.lock();
        if (Z < UB || UB==1)
            UB = Z;
        ZLock.unlock();


        /*printLock.lock();
        cout << endl <<" Successive Heuristic!!!" << endl;
        if(feasable)
            cout << "the solution is feasable!" << endl;
        else
            cout << "but the solution is not feasable! :/" << endl;
        cout << "Z = " << Z << endl;
        cout << "Z* = " << Zopt << endl << endl;
        printLock.unlock();*/

    }

    Succflag = false;
}





int main()
{
    cout << " Hello world!" << endl;
    string myList = "list.txt";
    ifstream list_file(myList.c_str(), ios::in);

    if(list_file)  // si l'ouverture a fonctionné
    {
        string name;  // déclaration d'une chaîne qui contiendra la ligne lue

        while (getline(list_file, name))  // tant que l'on peut mettre la ligne dans "contenu"
        {
            cout << endl << endl << endl << endl << endl << "Instance : " << name << endl ;

            lemonGraph G;

            unsigned n(0),m(0);

            vector < vector <int> > demands;

            vector < lemonDigraph* > d_graphs;
            vector < DIGRAPH_STR* > dg_strs;

            unsigned k = 3;


            Lagflag = true;
            Succflag = true;

            Population.clear();
            PopRanking.clear();
            for (unsigned j=0; j < PopOrgans.size(); j++)
                PopOrgans[j].clear();
            PopOrgans.clear();

            UB = 0;

            // ProcessTime example
            struct timeval startTime;
            struct timeval endTime;
            double tS, tE;
            // get the current time
            // - NULL because we don't care about time zone
            gettimeofday(&startTime, NULL);


            vector<lemonNode> LNodes;
            vector<lemonEdge> LEdges;

            lemonNodeIndexMap NIM(G,0);	//Pour numéroter les noeuds
            lemonNodeDegreeMap NDM(G,0);
            lemonNodeValMap NVM(G,0);

            lemonEdgeIndexMap EIM(G,0);	//Pour numéroter les arêtes
            lemonEdgeCapacityMap ECM(G,0); //Pour donner une capacité pour chaque arête
            lemonEdgeWeightMap EWM(G,0);
            lemonEdgeValMap EVM(G,0);

            GRAPH_STR g_str = {name,LNodes,LEdges,NIM,NDM,NVM,EIM,ECM,EWM,EVM};

            if (!open_TSP_graph(G,demands,n,m,g_str,k,name))
                return(10);

            n = countNodes(G);
            m = countEdges(G);

            lemonDiNodeIndexMap *N_NIM;	//Pour numéroter les noeuds
            lemonDiNodeDegreeMap *N_NDM;
            lemonDiNodeValMap *N_NVM;

            lemonArcIndexMap *N_AIM;	//Pour numéroter les arêtes
            lemonArcCapacityMap *N_ACM; //Pour donner une capacité pour chaque arête
            lemonArcWeightMap *N_AWM;
            lemonArcValMap *N_AVM;

            for (unsigned i=0; i < demands.size(); i++)
            {
                d_graphs.push_back(new lemonDigraph);

                PopOrgans.push_back(vector < lemonArcValMap* > ());

                vector<lemonDiNode> N_LNodes1, N_LNodes2;
                vector<lemonArc> N_LArcs;

                N_NIM = new lemonDiNodeIndexMap (*d_graphs[i],0);	//Pour numéroter les noeuds
                N_NDM = new lemonDiNodeDegreeMap (*d_graphs[i],0);
                N_NVM = new lemonDiNodeValMap (*d_graphs[i],9999999);

                N_AIM = new lemonArcIndexMap (*d_graphs[i],0);	//Pour numéroter les arêtes
                N_ACM = new lemonArcCapacityMap (*d_graphs[i],0); //Pour donner une capacité pour chaque arête
                N_AWM = new lemonArcWeightMap (*d_graphs[i],0);
                N_AVM = new lemonArcValMap (*d_graphs[i],0);

                dg_strs.push_back(new digraph_str ({N_LNodes1,N_LNodes2,N_LArcs,*N_NIM,*N_NDM,*N_NVM,*N_AIM,*N_ACM,*N_AWM,*N_AVM}));

                transform_graph(G, d_graphs[i], g_str, dg_strs[i], demands[i][0], demands[i][1]);
            }

            srand (time(NULL));

            /// CORPS DE L'ALGO HYBRID
            thread Lag(Lagrangien_Thread, ref(G), d_graphs, g_str, dg_strs, demands, k);
            thread Suc(Successive_Thread_Random, ref(G), d_graphs, g_str, dg_strs, demands, k);
            thread Gen(Genetic_Thread, ref(G), d_graphs, g_str, dg_strs);

            Suc.join();
            Lag.join();
            Gen.join();

            Population.clear();
            PopRanking.clear();
            for (unsigned j=0; j < PopOrgans.size(); j++)
                PopOrgans[j].clear();
            PopOrgans.clear();


            gettimeofday(&endTime, NULL);
            // calculate time in microseconds
            tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
            tE = endTime.tv_sec*1000000  + (endTime.tv_usec);

            cout << endl << endl << endl << " ********************* RESULTS *********************** "<< endl << endl << " Instance : " << name << endl << endl;

            cout << " *** UB = " << UB << endl;
            cout << " *** LB = " << LB << endl;

            cout << endl << "Le calcul de solution a dure : " << (tE - tS)/ 1000000.0L << endl << endl;
        }
    }
}
