import java.io.*;
import java.util.*;

public class MinimumViableProduct {
    
    // This code implements an algorithm to find a solution for the Minimum Viable Product (second part of the case)
    // The running time of the code is approximately 1 minute
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

        double[][] graph = constructGraph(connectionsList, citiesList, numCities);

        greedyAlgForMVP(graph, targetsList, connectionsList, citiesList);
        
    }


    /**
     * This method implements the greedy algorithm as specified for the Minimum Viable Product. First it is found which cities should be 
     * included to obtain a reach of 10. Then, this solution will be further optimised by adding most profitable targets.
     * @param graph The best network will be based upon this graph
    **/
    public static void greedyAlgForMVP(double[][] graph, ArrayList<Target> targetsList, ArrayList<Connection> connectionsList, ArrayList<City> citiesList) {

        // This is the set of terminals to be included, where a terminal will correspond to a a location of a city in a particular region.
        ArrayList<Integer> T = new ArrayList<Integer>();    
        T.add(64);      // To this set, we add eindhoven.

        ArrayList<Integer> regions = new ArrayList<Integer>();
        regions.add(10);           // add the region of eindhoven
        int count = 0;              // This count is used to activate the failsafe, which stops the while loop after a large number of iterations
        
        double reach = 1;
        int[][][] initialize = new int[graph.length][graph.length][4];             // This matrix only is defined to help initialize the networkSolution

        Network networkSolution = new Network(0, 0, 0, 0, initialize);

        while(reach < 10) {                     // First we make sure that the set of terminals is guaranteed to give a reach of 10
            int newCity1 = Integer.MAX_VALUE;
            double profit = Integer.MIN_VALUE;

            // We take one city of a target that maximizes the profit and is in a region that is not used yet!
            for(int j = 0; j < targetsList.size(); j++) {                   
                int originCity = targetsList.get(j).getOriginCity();        // We use origin city as this works the best for the algorithm -> it gives highest profit
                City temp = citiesList.get(originCity);
                if(!regions.contains(temp.getRegionNumber())) {             // If the region is not used yet in the network, then we consider what happens when we add the city

                    ArrayList<Integer> Ttemp = new ArrayList<Integer>();
                    for(int l = 0; l < T.size(); l++) {
                        Ttemp.add(T.get(l));
                    }
                    Ttemp.add(temp.getCityNumber());

                    // We find the network for the new set of terminals (Ttemp)
                    SteinerTreeProblem stTemp = new SteinerTreeProblem(Ttemp, graph);
                    int[][] stTempSolution = stTemp.findSteinerTree();

                    Network tempNetwork = createNetwork(stTempSolution, connectionsList, targetsList, citiesList);

                    // We make sure that all the edges are included in the network and we take the city for which profit is maximized
                    if((verifyEdges(tempNetwork, stTempSolution)) && (tempNetwork.getProfit() > profit)) {
                        profit = tempNetwork.getProfit();
                        newCity1 = temp.getCityNumber();
                    }
                }
            }

            // We add the most profitable city that also includes a new region, to our set of terminals.
            T.add(newCity1);

            // We find the network for the new set T 
            SteinerTreeProblem st = new SteinerTreeProblem(T, graph);
            int[][] stSolution = st.findSteinerTree();

            networkSolution = createNetwork(stSolution, connectionsList, targetsList, citiesList);
            reach = networkSolution.getNumRegions();                // update reach
            regions = networkSolution.findRegions(citiesList);      // update the set of regions that is contained in the network

            count++;
            //Failsafe #1
            if(count == 1000) {
                System.out.println("Failsafe #2");
                break;
            }
        }

        boolean nextIteration = true;   
        count = 0;                      // count for new failsafe

        // Now given the set of terminals found, we add most profitable targets to maximize profit.
        while(nextIteration) {
            
            double profit  = networkSolution.getProfit();
            int newTarget1 = Integer.MAX_VALUE;
            int newTarget2 = Integer.MIN_VALUE;

            boolean foundNewTargets = false;                // This boolean gives true if a new target is found which improves profit
            for(int i = 0; i < targetsList.size(); i++) {
                int originCity = targetsList.get(i).getOriginCity();
                int destinationCity = targetsList.get(i).getDestinationCity();

                ArrayList<Integer> Ttemp = new ArrayList<Integer>();
                for(int j = 0; j < T.size(); j++) {
                    Ttemp.add(T.get(j));
                }
                if(!Ttemp.contains(originCity)) {               // only add the target city if it is not yet contained in T
                    Ttemp.add(originCity);
                }
                if(!Ttemp.contains(destinationCity)) {          // only add the target city if it is not yet contained in T
                    Ttemp.add(destinationCity);
                }

                // Again we find a network for the new set of terminals (Ttemp)
                SteinerTreeProblem stTemp = new SteinerTreeProblem(Ttemp, graph);
                int[][] stTempSolution = stTemp.findSteinerTree();

                Network networkTemp = createNetwork(stTempSolution, connectionsList, targetsList, citiesList);

                // If by including the target the profit will be increased, then we update the temporary variables
                if((networkTemp.getProfit() > profit) && (verifyEdges(networkTemp, stTempSolution))) {      // verify edges to make sure that the network contains the necessary edges
                    profit = networkTemp.getProfit();
                    newTarget1 = originCity;
                    newTarget2 = destinationCity;
                    foundNewTargets = true;                     // we found a target in order to update the network, so boolean is true      
                }
            }

            // We add the new found targets to the actual set of terminals (only if T does not contain it yet)
            if(foundNewTargets) {                           
                if(!T.contains(newTarget1)) {
                    T.add(newTarget1);
                }
                if(!T.contains(newTarget2)) {
                    T.add(newTarget2);
                }

                // We create the network for the new set of terminals
                SteinerTreeProblem st = new SteinerTreeProblem(T, graph);
                int[][] stSolution = st.findSteinerTree();

                networkSolution = createNetwork(stSolution, connectionsList, targetsList, citiesList);
            }
            else {
                nextIteration = false;          // If we do not find targets that improve profit, then we stop the while loop
            }

            // Failsafe
            count++;
            if(count == 10000) {
                System.out.println("Failsafe #2");
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
        double score = 1*totalProfit - 0*totalEmission + 0*numRegions;   // calculates the score as specified in the assignment

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
     * necessary per modality. 
     * @return A (double) matrix with an edge between cities if the corresponding index (the weight) is greater than zero
    **/
    public static double[][] constructGraph(ArrayList<Connection> connectionsList, ArrayList<City> citiesList, int numCities) {

          // In this matrix, there exists an edge between city i and j, if avgGraph[i][j] > 0.
          double[][] avgGraph= new double[numCities][numCities];

          for(int i = 0; i < connectionsList.size(); i++) {
              Connection connectionTemp = connectionsList.get(i);
              double sum = 0;             // The sum keeps track of the total transportation units if all modalities are used
              double count = 0;           // The count keeps track of how many modalities are possible for the edge.                           
  
              if(connectionTemp.getCostPerTruck() > 0) {
                  sum += (double) (connectionTemp.getCostPerTruck());
                  count++;
              }
              if(connectionTemp.getCostPerSea() > 0) {
                  sum += (double) (connectionTemp.getCostPerSea());          
                  count++;
              }
              if(connectionTemp.getCostPerAir() > 0) {
                  sum += (double) (connectionTemp.getCostPerAir());          
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

