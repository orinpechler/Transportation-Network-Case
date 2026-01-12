import java.util.*;
import java.io.*;

public class ProofOfConcept {
// This is the code that implements the algorithm for the Proof of Concept. (The first part of the case)
// The running time is approximately 5 seconds.

    public static void main(String args[]) {

        // Read in the cities contained in the file
        ArrayList<City> citiesList = new ArrayList<City>();
        String filename1 = "cities.txt";
        try{
            citiesList = readCitiesFile(filename1);
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        int numCities = citiesList.size();

        // Read in the connections contained in the file
        ArrayList<Connection> connectionsList = new ArrayList<Connection>();
        String filename2 = "connections.txt";
        try{
            connectionsList = readConnectionsFile(filename2);
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // Read in the targets contained in the file
        ArrayList<Target> targetsList = new ArrayList<Target>();
        String filename3 = "targets.txt";
        try{
            targetsList = readTargetsFile(filename3);
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        double[][] avgGraph = constructGraph(connectionsList, numCities);

        greedyAlgForPOC(avgGraph, targetsList, connectionsList, citiesList);
       
    }


    /**
     * This method implements the greedy algorithm specified for the ProofOfConcept. Every iteration a new target, which has 
        the highest revenue, will be considered. When adding a target, the best network is found by using Steiner Trees.
        This is done until the profit is atleast 50.
     * @param avgGraph The best network will be based upon this graph
    **/
    public static void greedyAlgForPOC(double[][] avgGraph, ArrayList<Target> targetsList, ArrayList<Connection> connectionsList, ArrayList<City> citiesList) {
    
        // This is the set of terminals to be included, where a terminal will correspond to a location of a target.
        ArrayList<Integer> T = new ArrayList<Integer>();    
        T.add(64);      // To this set, we add eindhoven.

        targetsList = sortTargetList(targetsList);      // we sort the targets from most revenue to least

        int count = 0;              // This count is used to activate the failsafe, which stops the while loop after a large number of iterations
        
        double profit = 0;
        int[][][] initialize = new int[avgGraph.length][avgGraph.length][4];             // This matrix only is defined to help initialize the networkSolution

        Network networkSolution = new Network(0, 0, 0, 0, initialize);  // Just initialization

        while(profit < 50) {                    // find solution such that profit is 50
            int newTarget1 = Integer.MAX_VALUE;
            int newTarget2 = Integer.MAX_VALUE;

            // We take the next most profitable target and see whether a network that includes also this target next to 
            // the other targets is possible.
            for(int i = 0; i < targetsList.size(); i++) {
                int originCity = targetsList.get(i).getOriginCity();
                int destinationCity = targetsList.get(i).getDestinationCity();

                // We create a temp of T, which the new cities will be added to 
                if((!T.contains(originCity)) && (!T.contains(destinationCity))) {
                    ArrayList<Integer> Ttemp = new ArrayList<Integer>();
                    for(int j = 0; j < T.size(); j++) {
                        Ttemp.add(T.get(j));
                    }
                    Ttemp.add(originCity);
                    Ttemp.add(destinationCity);

                    // We find the network for the new set of terminals (Ttemp)
                    SteinerTreeProblem stTemp = new SteinerTreeProblem(Ttemp, avgGraph);
                    int[][] tempSolution = stTemp.findSteinerTree();
                    Network networkTemp = createNetwork(tempSolution, connectionsList, targetsList, citiesList);

                    // We see whether it is possible to get a network including the new target (the network has to include all edges
                    // that are also included in the Steiner Tree)
                    if((verifyEdges(networkTemp, tempSolution))) {              // verify edges to make sure that the network contains the necessary edges
                        newTarget1 = originCity;
                        newTarget2 = destinationCity;
                        break;
                    }
                } else {
                    continue;
                }
            }

            // We add the targets found to the set of terminals
            T.add(newTarget1);
            T.add(newTarget2);
            
            // For which we find again the steiner tree and the corresponding network.
            SteinerTreeProblem st = new SteinerTreeProblem(T, avgGraph);
            int[][] stSolution = st.findSteinerTree();
            
            networkSolution = createNetwork(stSolution, connectionsList, targetsList, citiesList);
            profit = networkSolution.getProfit();               // update profit

            // Failsafe
            count++;
            if(count == 1000) {
                System.out.println("Failsafe");
                break;
            }
        }
        // After the while loop the optimal solution is found and we output this to the terminal
        networkSolution.outputNetwork(citiesList);

    }


    /** 
     * This method verifies whether the edges that are in the steiner tree are all included in the resulting network
     * @param network The network corresponding to the Steiner Tree
     * @param stSolution The Steiner Tree used to compute the network
     * @return True if the network includes all edges as specified in the Steiner Tree
     **/ 
    public static boolean verifyEdges(Network network, int[][] stSolution) {

        int[][][] networkGraph = network.getNetwork();

        // We simplify the network to a graph, with has a value 1 if the edge is included
        int[][] edgesIncluded = new int[networkGraph.length][networkGraph.length];
        for(int i = 0; i < networkGraph.length; i++) {
            for(int j = 0; j < networkGraph.length; j++) {
                for(int k = 0; k < 4; k++) {
                    if(networkGraph[i][j][k] > 0) {
                        edgesIncluded[i][j] = 1;
                    }
                }
            }
        }

        // Check if all edges included in the steiner tree and the network are the same
        boolean flag = true;
        for(int i = 0; i < edgesIncluded.length; i++) {
            for(int j = 0; j < edgesIncluded.length; j++) {
                if(edgesIncluded[i][j] != stSolution[i][j]) {
                    flag = false;
                }
            }
        }

        return flag;
    }


    /**
     * This method sorts the target list from a target with highest revenue to a target with lowest revenue
     * @param targetsList It will sort this list
     * @return A sorted target list
    **/
    public static ArrayList<Target> sortTargetList(ArrayList<Target> targetsList) {

        // This will be the sorted list
        ArrayList<Target> sortedTargets = new ArrayList<Target>();

        while(targetsList.size() > 0) {

            int maxRev = Integer.MIN_VALUE;
            Target tempTarget = new Target(0, 0, 0);
            
            // Find target with next highest revenue
            for(int j = 0; j < targetsList.size(); j++) {
                if(targetsList.get(j).getRevenue() > maxRev) {
                    maxRev = targetsList.get(j).getRevenue();
                    tempTarget = targetsList.get(j);
                }
            }

            sortedTargets.add(tempTarget);          // add to sorted list
            targetsList.remove(tempTarget);         // remove from original list so it will not be considered anymore
        }

        return sortedTargets;
    }


    /** 
     * This method creates a network for a given Steiner Tree, it does so greedily by taking for each connection, the modality
     * that uses lowest transportation units.
     * @param stSolution This is the Steiner Tree for which the network is created
     * @param terminals The revenue is calculated based upon which targets are used as terminals
     * @return It returns the network as an object Network
     **/  
    public static Network createNetwork(int[][] stSolution, ArrayList<Connection> connectionsList, ArrayList<Target> targetsList, ArrayList<City> citiesList) {

        // Initializing the budget (in terms of transportation units) for each modality, the budget for each 
        // modality is saved as follows: Truck(0), Sea(1), Air(2), Rail(3)
        int[] budget = new int[4];
        for(int i = 0; i < 4; i++) {
            budget[i] = 50;
        }

        // This array will keep track of the number of times is invested in a new modality, i.e. when 5 more units of 
        // transportation are bought for: Truck(0), Sea(1), Air(2), Rail(3), respectively.
        int[] numberOfInv = new int[4];

        // The networkEdges will save the network as specified for the Network object
        int n = stSolution.length;
        int[][][] networkEdges = new int[n][n][4];

        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(stSolution[i][j] > 0) {              // If greater than 0, than we need to include this connection in the network
                    Connection edge = findConnection(connectionsList, i, j);
                    int[] sortedTrans = edge.getMinTrans();         // Sorts the transportation possibilies, as specified in Connection
                    for(int l = 0; l < sortedTrans.length; l++) {
                        int modality = sortedTrans[l];
                        int costModality;
                        if(modality == 0) {
                            costModality = edge.getCostPerTruck();      // cost if modality is truck
                        } else if((modality == 1)) {
                            costModality = edge.getCostPerSea();        // cost if modality is ship
                        } else if(modality == 2) {
                            costModality = edge.getCostPerAir();        // cost if modality is airplane
                        } else { 
                            costModality = edge.getCostPerRail();       // cost if modality is train
                        }

                        if(costModality <= budget[modality]) {          // if the budget already allows to use the modality, then we use it!
                            networkEdges[i][j][modality] = costModality;
                            budget[modality] -= costModality;           // update budget
                            break;                                      // only use one modality
                        }
                        else {                                          // else we try to invest in order to use the modality
                            boolean invested = false;
                            for(int m = 1; m <= (5 - numberOfInv[modality]); m++) {         // we can invest in total at most 5 times
                                if(costModality <= (m*5 + budget[modality])) {              // every investment counts for 5 units of transportation
                                    networkEdges[i][j][modality] = costModality;
                                    numberOfInv[modality] += m;
                                    budget[modality] = budget[modality] + m*5 - costModality;       // update budget taking into account the investment
                                    invested = true;
                                    break;                                                  // only use one modality
                                }
                            }
                            if(invested) {                          
                                break;
                            } else{                                     // if we also can not invest then we try the next mode of transport
                                continue;
                            }
                        }
                    }
                }
            }
        }
        

        double totalProfit = calculateProfit(targetsList, networkEdges, numberOfInv);
        double totalEmission = calculateEmission(connectionsList, networkEdges);
        int numRegions = findNumRegions(citiesList, networkEdges);                  
        double score = 0*totalProfit - 1*totalEmission + 0*numRegions;   // calculates the score as specified in the assignment

        Network network = new Network(score, totalProfit, totalEmission, numRegions, networkEdges); 
        
        return network;
    }  


    /**
     * This method calculates for a given network and the investments made for the network, the profit. This is equal to the 
     * total revenue obtained from the targets minus the investment costs.
     * @return Profit as a double
    **/
    public static double calculateProfit(ArrayList<Target> targetsList, int[][][] networkEdges, int[] numberOfInv) {

        // This arraylist saves the citynumbers of the cities used in the network
        ArrayList<Integer> citiesInNetwork = new ArrayList<Integer>();

        // Revenue is equal to the summation of the profit gained from including targets
        double revenue = 0;

        for(int i = 0; i < networkEdges.length; i++) {
            for(int j = 0; j < networkEdges.length; j++) {
                for(int k = 0; k < 4; k++) {
                    if(networkEdges[i][j][k] > 0 ) {            // if the two cities are connected 
                        if(!citiesInNetwork.contains(i)) {      // and if the city i is not yet in the list
                            citiesInNetwork.add(i);             // we add the city
                        }
                        if(!citiesInNetwork.contains(j)) {      // same for city j
                            citiesInNetwork.add(j);
                        }
                    }
                }
            }
        }

        for(int i = 0; i < targetsList.size(); i++) {
            int originCity = targetsList.get(i).getOriginCity();
            int destinationCity = targetsList.get(i).getDestinationCity();
            if(citiesInNetwork.contains(originCity) && citiesInNetwork.contains(destinationCity)) {     // We check if a target is contained in the network
                revenue += targetsList.get(i).getRevenue();                         // If so, then we add it's revenue
            }
        }
        
        // The cost from investment follows from the table as specified in the assignment
        // Observe that all the calculations are simplified to formulas for finding the total investment cost per modality
        double totalInv = 0;
        totalInv += 3*numberOfInv[0];            // total cost of investing in capacity for Truck

        // Total cost of investing in capacity for Ship:
        int costInvShip = 0;
        for(int i = 1; i <= (numberOfInv[1]); i++) {
            costInvShip += i;
        }
        totalInv += costInvShip;

        // Total cost of investing in Airplane:
        totalInv += (numberOfInv[2]*numberOfInv[2]);

        // Total cost of investing in Train:
        if(numberOfInv[3] > 0) {
            totalInv += Math.pow(2, (numberOfInv[3] - 1));
        }
        
        // Profit is then revenue minus total investment
        double profit = revenue - totalInv;
        return profit;
    }


    /**
     * This method calculates the total emission for the network that has been found.
     * @return Total Emission as a double
    **/
    public static double calculateEmission(ArrayList<Connection> connectionsList, int[][][] networkEdges) {

        double emission = 0;
        for(int i = 0; i < networkEdges.length; i++) {
            for(int j = 0; j < networkEdges.length; j++) {
                for(int k = 0; k < 4; k++) {
                    if(networkEdges[i][j][k] > 0) {
                        if(k == 0) {
                            emission += networkEdges[i][j][k]*8;        // Emission of connection if modality is Truck
                        } else if(k == 1) {
                            emission += networkEdges[i][j][k]*1;        // Emission of connection if modality is Ship
                        } else if(k == 2) {
                            emission += networkEdges[i][j][k]*25;        // Emission of connection if modality is Plane
                        } else {
                            emission += networkEdges[i][j][k]*3;        // Emission of connection if modality is Train
                        }
                    }
                }
            }
        }

        return emission;
    }


    /**
     * This method finds the number of regions included in a network 
     * @return Number of regions as a integer
     */
    public static int findNumRegions(ArrayList<City> citiesList, int[][][] networkEdges) {

        // We add all regions visited to an arraylist
        ArrayList<Integer> regions = new ArrayList<Integer>();
        for(int i = 0; i < networkEdges.length; i++) {
            for(int j = 0; j < networkEdges.length; j++) {
                for(int k = 0; k < 4; k++) {
                    if(networkEdges[i][j][k] > 0) {
                        int region1 = citiesList.get(i).getRegionNumber();
                        int region2 = citiesList.get(j).getRegionNumber();

                        if(!regions.contains(region1)) {        // Only add it if it was not in there yet
                            regions.add(region1);
                        }
                        if(!regions.contains(region2)) {
                            regions.add(region2);
                        }
                    }
                }
            }
        }

        int numRegions = regions.size();    // Number of regions is the size of the array
        return numRegions;
    }


    /**
     * This method finds the connection of interest given two cities, it returns this is an object Connection
     * @param cityNumber1 The city is either the origin city or the destination city
     * @param cityNumber2 The city is either the origin city or the destination city
     * @return The Connection of the two cities
    **/
    public static Connection findConnection(ArrayList<Connection> connectionsList, int cityNumber1, int cityNumber2) {
        
        Connection edge = new Connection(0,0,0,0,0,0);
        for(int i = 0; i < connectionsList.size(); i++) {
            Connection temp = connectionsList.get(i);
            if(((temp.getOriginCity() == cityNumber1) && (temp.getDestinationCity() == cityNumber2)) || ((temp.getOriginCity() == cityNumber2) && (temp.getDestinationCity() == cityNumber1))) {
                edge = temp;
                break;
            }
        }

        return edge;
    }


    /**
     * This method creates a matrix of the graph of connections. In this graph, we take the average of the transportation units
     * necessary per modality. We multiply the cost for truck and airplane by 2, as a penalty. The reason for this is that the truck and
     * ship modality use the most emission.
     * @return A (double) matrix with an edge between cities if the corresponding index (the weight) is greater than zero
    **/
    public static double[][] constructGraph(ArrayList<Connection> connectionsList, int numCities) {

        // In this matrix, there exists an edge between city i and j, if avgGraph[i][j] > 0.
        double[][] avgGraph= new double[numCities][numCities];

        for(int i = 0; i < connectionsList.size(); i++) {
            Connection connectionTemp = connectionsList.get(i);
            double sum = 0;             // The sum keeps track of the total transportation units if all modalities are used
            double count = 0;           // The count keeps track of how many modalities are possible for the edge.                           

            if(connectionTemp.getCostPerTruck() > 0) {
                sum += (double) (connectionTemp.getCostPerTruck()*2);       // penalty (second highest emission modality)
                count++;
            }
            if(connectionTemp.getCostPerSea() > 0) {
                sum += (double) (connectionTemp.getCostPerSea());          
                count++;
            }
            if(connectionTemp.getCostPerAir() > 0) {
                sum += (double) (connectionTemp.getCostPerAir()*2);          // penalty (highest emission modality)
                count++;
            }
            if(connectionTemp.getCostPerRail() > 0) {
                sum += (double) (connectionTemp.getCostPerRail());
                count++;
            }

            double avg = sum/count;         // Note that count makes sure that we only divide by the number of modalities that are possible.

            // We add the weight for both ways of the connection
            avgGraph[connectionTemp.getOriginCity()][connectionTemp.getDestinationCity()] = avg;
            avgGraph[connectionTemp.getDestinationCity()][connectionTemp.getOriginCity()] = avg;
        }

        return avgGraph;
    }


    /** 
     * This method reads in the data of cities given in the data file. It returns the data in the form of an 
     * arraylist, where each city is stored as an object: City.
    **/ 
    public static ArrayList<City> readCitiesFile(String filename) throws FileNotFoundException {

        ArrayList<City> citiesList = new ArrayList<City>();

        File file = new File(filename);
        Scanner input = new Scanner(file);
        int numCities = input.nextInt();

        for(int i = 0; i < numCities; i++) {
            int cityNumber = input.nextInt() - 1;       // Convert it in a number from 0 to (numCities - 1) = 75 in our case
            String cityName = input.next();
            double latitude = input.nextDouble();
            double longitude = input.nextDouble();
            int regionNumber = input.nextInt();

            City cityTemp = new City(cityNumber, cityName, latitude, longitude, regionNumber);
            citiesList.add(cityTemp);
        }

        input.close();

        return citiesList;
    }


    /**
     * This method reads in the data of connections given in the data file. It returns the data in the form of an 
     * arraylist, where each connection is stored as an object: Connection.
    **/ 
    public static ArrayList<Connection> readConnectionsFile(String filename) throws FileNotFoundException {

        ArrayList<Connection> connectionsList = new ArrayList<Connection>();

        File file = new File(filename);
        Scanner input = new Scanner(file);
        int numConnections = input.nextInt();

        for(int i = 0; i < numConnections; i++) {
            int originCity = input.nextInt() - 1;           // Convert it in a number from 0 to (numCities - 1) = 75 in our case
            int destinationCity = input.nextInt() - 1;
            int costPerTruck = input.nextInt();
            int costPerSea = input.nextInt();
            int costPerAir = input.nextInt();
            int costPerRail = input.nextInt();

            Connection connectionTemp = new Connection(originCity, destinationCity, costPerTruck, costPerSea, costPerAir, costPerRail);
            connectionsList.add(connectionTemp);
        }

        input.close();

        return connectionsList;
    }


    /**
     * This method reads in the data of targets given in the data file. It returns the data in the form of an 
     * arraylist, where each target is stored as an object: Target.
    **/
    public static ArrayList<Target> readTargetsFile(String filename) throws FileNotFoundException {

        ArrayList<Target> targetsList = new ArrayList<Target>();

        File file = new File(filename);
        Scanner input = new Scanner(file);
        int numTargets = input.nextInt();

        for(int i = 0; i < numTargets; i++) {
            int originCity = input.nextInt() - 1;               // Convert it in a number from 0 to (numCities - 1) = 75 in our case
            int destinationCity = input.nextInt() - 1;
            int revenue = input.nextInt();

            Target targetTemp = new Target(originCity, destinationCity, revenue);
            targetsList.add(targetTemp);
        }

        input.close();

        return targetsList;
    }
}
