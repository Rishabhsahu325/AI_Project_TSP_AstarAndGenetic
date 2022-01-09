#include<bits/stdc++.h>
#include<iostream>
#include <fstream>
using namespace std;

class Graph
{

public:
    int vertexCount;
    vector<vector<int>> dist;
    Graph(int vertices)//constructor for graph class object
    {   
        vertexCount = vertices;
        dist = vector<vector<int>>(vertices, vector<int>(vertices, INT_MAX));
        for (int i = 0; i < vertices; i++)
            dist[i][i] = 0;
    }

    void addEdge(int city1, int city2, int weight)//add an edge connecting vertices u and v with weight w
    {
        dist[city1][city2] = weight;
        dist[city2][city1] = weight;
    }

    void printGraph()//display all graph weights in adjacency matrix format
    {
        for (int i = 0; i < vertexCount; i++)
        {
            for (int j = 0; j < vertexCount; j++)
            {
                cout << dist[i][j] << "\t";
            }
            cout << endl;
        }
    }

    int getEdgeWt(int city1, int city2)
    {
        return dist[city1][city2];
    }
};

class FringeListNode // fringe list
{
    public:
    int nodeId;
    vector<int> path;
    int fValue, gValue;
    unordered_set<int> visitedNodes;
//constructor with parameters to store evaluation function f ,node cost function g and heuristic cost function h
    FringeListNode(int nodeId = -1, int fValue = INT_MAX, int gValue = INT_MAX, vector<int> path = {}, unordered_set<int> visitesNodes = {})
    {
        this->fValue = fValue;
        this->gValue = gValue;
        this->nodeId = nodeId;
        this->path = path;
        this->visitedNodes = visitedNodes;
    }

    int getFValue()
    {
        return this->fValue;
    }
};

struct CompareFLNode //overloading operator to comparetwo fringe list nodes based on their fvalue
{
    bool operator()(FringeListNode &n1, FringeListNode &n2)
    {
        return n1.getFValue() > n2.getFValue();
    }
};

class AStar
{
    
public:
    int numOfCities;
    vector<vector<int>> dist;
    priority_queue<FringeListNode, vector<FringeListNode>, CompareFLNode> fringeList;
    int startCity;
    FringeListNode solution;
    int cost=0;
    AStar(Graph g, int startNode)// constructor for A star class object created from a graph object and starting node
    {
        this->dist = g.dist;
        this->startCity = startNode;
        this->numOfCities = g.vertexCount;
        solution = FringeListNode();
        cost = 0;
    }

    void calculateCost()//save the cost of solution path to cost variable
    {   
        vector<int> path = solution.path;
        int n = path.size();
        for (int i = 1; i < path.size(); i++)
        {
            this->cost += dist[path[i - 1]][path[i]];
        }
        this->cost += dist[path[n - 1]][startCity];
    }

    vector<int> getAdjacentNodes(FringeListNode curr)//return all the unvisited nodes that are connected to Node curr
    {
        vector<int> successorNodes;
        unordered_set<int> visited = curr.visitedNodes;
        for (int i = 0; i < numOfCities; i++)
        {
            if (i != curr.nodeId && visited.find(i) == visited.end())
            {
                successorNodes.push_back(i);
            }
        }
        return successorNodes;
    }

    int calculateHeuristic(FringeListNode curr)//calculate MST heuristic values based on Prim's Algorithm
    {
        // We will consider only unvisited nodes (except current) for calcuating MST
        // h = min cost edge from curr node to any node in thr MST + MST + min cost edge from any node in the MST to startNode

        
        if (curr.path.size() == numOfCities - 1)
        {
            return 0;
        }
        unordered_set<int> unvisitedNodes; //nodes that are yet to be considered for MST
        
        for (int i = 0; i < numOfCities; i++)
        {
            if (i != curr.nodeId && curr.visitedNodes.find(i) == curr.visitedNodes.end())
            {
                unvisitedNodes.insert(i);
            }
        }
       

        int MSTCost = 0;
        if (unvisitedNodes.size() > 0)
        {
            MSTCost = getMSTCost(unvisitedNodes);
        }

        int currToMST = INT_MAX;
        int MSTToStart = INT_MAX;
        for (int node : unvisitedNodes)
        {
            currToMST = min(currToMST, dist[curr.nodeId][node]);
            MSTToStart = min(MSTToStart, dist[node][startCity]);
        }
        return currToMST + MSTCost + MSTToStart;
    }

    bool solve()
    {
        //initial state i.e startNode
        fringeList.push( FringeListNode( startCity, 0, 0, {}, {} ) );
        

        // run while we have states in the fringelist or until we find the optimal solution
        while (!fringeList.empty())
        {

            // select node with least f cost
            FringeListNode curr = fringeList.top();
            fringeList.pop();
            
            //Goal test
            if (isGoalState(curr))
            {
                solution = curr;
                solution.path.push_back(curr.nodeId);
                calculateCost();
                
                return true;
            }
            

            //get all adjacent nodes of current node i.e all - {already visited, itself}
            vector<int> successorsNodes = getAdjacentNodes(curr);
            
    
            // for each successor node insert corresponding state in the fringeList
            for (int successor : successorsNodes)
            {
                
                FringeListNode successorState;
                successorState.nodeId = successor;
                successorState.gValue = curr.gValue + dist[curr.nodeId][successor];
                successorState.path = curr.path;
                successorState.path.push_back(curr.nodeId); //add current node id to the path of successor
                successorState.visitedNodes = curr.visitedNodes;
                successorState.visitedNodes.insert(curr.nodeId);   //similarly add current node id to the visited set
                int h = calculateHeuristic(successorState);        // MST heuristic
                successorState.fValue = successorState.gValue + h; // f = g + h
                fringeList.push(successorState);
                
            }
        }
        return false;
    }

    bool isGoalState(FringeListNode curr)
    {
        
        if (curr.path.size() == numOfCities - 1)
            return true;
        return false;
    }

    int getMSTCost(unordered_set<int> nodeList)
    {
        
        if (nodeList.size() == 1)
        {
            return 0;
        }
        
        int key[numOfCities];
        int visited[numOfCities];
        for (int i = 0; i < numOfCities; i++)
        {
            key[i] = INT_MAX;
            visited[i] = 0;
        }
        int src = *nodeList.begin();
        key[src] = 0;
        int cost = 0;
        for (int i = 0; i < nodeList.size(); i++)
        {
            int u = minKey(key, visited, nodeList);
            cost += key[u];
            visited[u] = 1;
            for (int j = 0; j < numOfCities; j++)
            {
                if (nodeList.find(j) != nodeList.end() && !visited[j] && dist[u][j] < INT_MAX)
                {
                    key[j] = min(key[j], dist[u][j]);
                }
            }
        }
        
        return cost;
    }

    int minKey(int key[], int visited[], unordered_set<int> nodeList)
    {

        int mn = INT_MAX;
        int min_index;
        
        for (int i = 0; i < numOfCities; i++)
        {
            if (nodeList.count(i) && visited[i] == 0 && key[i] < mn)
            {
                mn = key[i];
                min_index = i;
            }
        }
        return min_index;
    }
};


void addEdge(vector<int> weights[], int city1, int city2)
{
	weights[city1].push_back(city2);
	weights[city2].push_back(city1);
}

// A utility function to print the adjacency list
// representation of graph
void printGraph(vector<int> weights[], int V)
{
	for (int v = 0; v < V; ++v)
	{
		cout << "\n Adjacency list of vertex "<< v << "\n head ";
		for (auto x : weights[v])
		cout << "-> " << x;
		printf("\n");
	}
}


Graph readMap(string filename)//add city data from city coordinates to Graph object
{
    fstream myfile(filename, std::ios_base::in);
    string check;
    int citiesCount;

    myfile >> check;
    while(strcmp(check.c_str(),"DIMENSION")!=0 )
    {
        myfile>>check;
    }
    char skipCharacter;
    myfile>>skipCharacter;
    myfile>>citiesCount;//read number of cities
    
    //skip few lines not containing coordinates
    myfile>>check;
    while(strcmp(check.c_str(),"NODE_COORD_SECTION") !=0 )
    {
        myfile>>check;
    }
    myfile>>skipCharacter;
    
    Graph g(citiesCount);
    tuple<int,double,double>entry;//tuple to store coordinate values of a single city
    vector<tuple<int,double,double>> coordinates;
    for(int city=0;city<citiesCount;city++)//added all coordinates
    {
        myfile>>get<0>(entry) >> get<1>(entry) >> get<2>(entry);
        coordinates.push_back(entry);
    }
    int distance=0;
    //calculating intercity distances and adding them as weights
    for(int c1=0;c1<citiesCount-1;c1++)
    {
        g.addEdge(c1,c1,0);
        for(int c2=c1+1;c2<citiesCount;c2++)
        {
            double diff1=get<1>(coordinates[c1])-get<1>(coordinates[c2]);
            double diff2=get<2>(coordinates[c1])-get<2>(coordinates[c2]);
            distance=(int)ceil(sqrt(( diff1*diff1)+(diff2*diff2)));
            g.addEdge(c1,c2,distance);
            g.addEdge(c2,c1,distance);
        }
    }
    g.addEdge(citiesCount-1,citiesCount-1,0);
    return g;
}

int main(int argc ,char* argv[])
{
    string fname=argv[1];
    Graph map=readMap("./"+fname);
    //cout<<"\n"<<argv[0]<<" "<<argv[1]<<"\n";
    int startCity = 0;
    float startTicks=(float)clock();//to only account for time spent in solving the problem without storing input values
    AStar solver(map, startCity);
    bool found = solver.solve();
    if (found)
    {
        vector<int> path = solver.solution.path;
        int cost = solver.cost;
        cout << "Cost = " << cost << endl;
        cout << "Path:\n";
        for (int city : path)
        {
            cout << city + 1 << " -> ";//as city indexes start from 1 
        }
        cout << endl;
    }
    else
    {
        cout << "!!! PROBLEM !!!" << endl;
    }
    cout << "time taken : " << ((float)clock() -startTicks)/ CLOCKS_PER_SEC << " secs " << endl;//display time taken to find optimal path
    return 0;
}

