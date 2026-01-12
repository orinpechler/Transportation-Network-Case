import java.util.*;

public class SteinerTreeProblem {
    
    private ArrayList<Integer> terminals;
    private double[][] graph;

    public SteinerTreeProblem(ArrayList<Integer> terminals, double[][] graph) {
        this.terminals = terminals;
        this.graph = graph;
    }


    /**
     * This method implements the Shortest Path-based approximation algorithm for the Steiner Tree problem. It returns the Steiner Tree
     * as a binary matrix, where index is equal to 1 if an edge exists for the corresponding row and column. (Note: if (i,j) = 1, then 
     * (j,i) = 0, to simplify calculations)
     * @return Integer Matrix
    **/
    public int[][] findSteinerTree() {

        // First the algorithm starts with a subtree T consisting of one given terminal vertex: in our case we always start with Eindhoven, which is 64 (in java terms)
        // Note that in the algorithm the first terminal in the terminals list will be Eindhoven
        ArrayList<Integer> T = new ArrayList<Integer>();        // T represents all vertices that are included in the steiner tree      
        T.add(this.terminals.get(0));
        int[][] steinerTree = new int[this.graph.length][this.graph.length];    // This matrix represents the graph of the steiner tree
        
        ArrayList<Integer> terminalsUsed = new ArrayList<Integer>();        // This stores the terminals already considered by the algorithm
        terminalsUsed.add(this.terminals.get(0));                     

        Collections.sort(this.terminals);       

        // We update the steiner tree with shortest paths until all terminals are included in the steiner tree
        while(!this.terminals.equals(terminalsUsed)) {

            ArrayList<Double> shortestPath = new ArrayList<Double>();
            double shortestLength = Integer.MAX_VALUE;
            int terminal = Integer.MAX_VALUE;

            // In this for loop we try to find the vertex already included in the steiner tree, which is closest 
            // to any terminal not yet included. (also this terminal has to be found)
            for(int i = 0; i < this.terminals.size(); i++) {
                int terminalTemp = this.terminals.get(i);
                if(!terminalsUsed.contains(terminalTemp)) {         // only consider terminals not in the steiner tree!
                    for(int j = 0; j < T.size(); j++) {
                        int vertexTemp = T.get(j);
                        ArrayList<Double> spTemp = shortestPath(vertexTemp, terminalTemp);
                        if(spTemp.get(0) < shortestLength) {
                            shortestPath = spTemp;
                            shortestLength = spTemp.get(0);   // length of shortest path is saved at the first index
                            terminal = terminalTemp;
                        }
                    }
                }
            }

            // Add terminal which is found, to the used terminals
            terminalsUsed.add(terminal);

            // Add all vertices included in the path to the steiner tree (and the edges to the graph)
            for(int k = 1; k < shortestPath.size() - 1; k++) {
                int vertexTemp2 = shortestPath.get(k).intValue();
                int vertexTemp3 = shortestPath.get(k+1).intValue();
                if(!T.contains(vertexTemp2)) {
                    T.add(vertexTemp2);
                }
                if(!T.contains(vertexTemp3)) {
                    T.add(vertexTemp3);
                }
                // In the steiner tree graph we, say that there exists an edge between vertexTemp2 and vertexTemp3
                steinerTree[vertexTemp2][vertexTemp3] = 1;
                //steinerTree[vertexTemp3][vertexTemp2] = 1; is not included, to simplify calculations for cost etc.
            }

            Collections.sort(terminalsUsed);            // Sort arraylist to check equivalence
        }

        return steinerTree;
    }


    /**
     * This method implements the Dijkstra algorithm for the shortest path problem. For this a start node and an end node is given,
     * together with a matrix representing the graph. In this matrix the weight corresponds to an entry in the matrix,
     * where if an edge does not exist, then the weight is equal to 0.
     * @return ArrayList of doubles
    **/
    public  ArrayList<Double> shortestPath(int startNode, int endNode) {
        // The method assumes that the startnode and endnote take into account that the first node is numbered by 0

        int n = this.graph.length;               // Number of vertices
        double[] labels = new double[n];
        for(int i = 0; i < n; i++) {
            labels[i] = Integer.MAX_VALUE;
        }
        labels[startNode] = 0;

        // Arraylist that represents all vertices
        ArrayList<Integer> V = new ArrayList<Integer>();    
        for(int i = 0; i < n; i++) {
            V.add(i);
        }

        // Arraylist that represents the vertices that are already considered by the algorithm
        ArrayList<Integer> S = new ArrayList<Integer>();

        while(!S.equals(V)) {

            double minLabel = Integer.MAX_VALUE;        // To find the vertex with the minimum label in V\S
            int vertex = Integer.MAX_VALUE;             // Vertex with the minimum label
            for(int i = 0; i < V.size(); i++) {
                int tempVertex = V.get(i);
                if(!S.contains(tempVertex)){ 
                    if(labels[tempVertex] < minLabel) {
                        minLabel = labels[tempVertex];
                        vertex = tempVertex;
                    }
                }
            }

            S.add(vertex);              // The vertex will now be considered
            Collections.sort(S);        // Sort S such that it can be compared with V      

            // Now we update the labels for all outgoing edges from the vertex
            for(int i = 0; i < n; i++) {
                if(this.graph[vertex][i] > 0) {          // Can only update if the edge exists
                    double edgeWeight = this.graph[vertex][i];
                    if(labels[i] > (labels[vertex] + edgeWeight)) {
                        labels[i] = labels[vertex] + edgeWeight;
                    }
                }
            }
        }

        // Now we find the path of vertices that lead the path from the start node to the end node
        // This is saved in an arraylist, where the first element specifies the length of the shortest path!
        ArrayList<Double> shortestPath = new ArrayList<Double>();
        double pathLength = labels[endNode];
        shortestPath.add(pathLength);           // Saves length of the shortest path
        shortestPath.add( (double) endNode);    // Add the end node

        int node = endNode;
        while(pathLength != 0) {                // We work our way backwards from the endnode
            for(int i = 0; i < n; i++) {
                double edgeWeight = this.graph[i][node];
                if((edgeWeight > 0) && (pathLength == edgeWeight + labels[i])) {
                    pathLength = labels[i];
                    node = i;
                    shortestPath.add(1, (double) i);        // As we work backwards, we add this vertex before all the others
                    break;
                }
            }
        }

        return shortestPath;
    }

}
