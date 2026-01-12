import java.util.*;

public class Network {
    
    private double score;           // This will save the score of the network
    private double totalProfit;     // This will save the total profit of the network
    private double totalEmission;   // This will save the total emission of the network
    private int numRegions;         // This will save the number of regions the network passes through
    private int[][][] networkEdges; 
    // This 3d matrix will save the connections in the matrix in the following manner:
    // The row represents the origin city, and the column represents the destination city. The third dimension
    // will save whether this connection is made possible via: Truck (0), Sea(1), Air(2), Rail(3). The value saved 
    // at networkEdges[i][j][k] = the capacity used by the modality. Hence, if a connection is not made or
    // if the modality for the connection is not used, this is equal to 0
    // Observe that i,j = 0,..,n-1 (n number of cities), and k= 0,..,3

    public Network(double score, double totalProfit, double totalEmission, int numRegions, int[][][] networkEdges) {

        this.score = score;
        this.totalProfit = totalProfit;
        this.totalEmission = totalEmission;
        this.numRegions = numRegions;
        this.networkEdges = networkEdges;
    }


    /**
     * Returns the score of the network
     * @return double
    **/
    public double getScore() {
        return this.score;
    }


    /**
     * Returns the total profit of the network, i.e. revenue - total investment
     * @return double
    **/
    public double getProfit() {
        return this.totalProfit;
    }


    /**
     * Returns the total emission of the network
     * @return double
    **/
    public double getEmission() {
        return this.totalEmission;
    }


    /**
     * Returns the number of regions for the network
     * @return double
    **/
    public double getNumRegions() {
        return this.numRegions;
    }


    /**
     * Return the 3d matrix representing the actual network
     * @return 3D Integer Matrix
    **/
    public int[][][] getNetwork() {
        return this.networkEdges;
    }


    /**
     * Returns an arraylist of integers representing the regions that are contained in the network
     * @return ArrayList of Integers
     */
    public ArrayList<Integer> findRegions(ArrayList<City> citiesList) {

        // We first find the cities contained in the network
        ArrayList<City> cities = new ArrayList<City>();
        for(int i = 0; i < this.networkEdges.length; i++) {
            for(int j = 0; j < this.networkEdges.length; j++) {
                for(int k = 0; k < 4; k++) {
                    if(this.networkEdges[i][j][k] > 0) {            // If a connection exists, then both cities are in the network
                        City temp1 = citiesList.get(i);             // The citynumber corresponds to its position in the cities list
                        City temp2 = citiesList.get(j);
                        if(!cities.contains(temp1)) {               // If not yet added to the cities in the network
                            cities.add(temp1);                      // Then add it
                        }
                        if(!cities.contains(temp2)) {               // Same for the other city  
                            cities.add(temp2);
                        }
                    }
                }
            }
        }

        // Now we find the regions that are used in the network
        ArrayList<Integer> regions = new ArrayList<Integer>();
        for(int i = 0; i < cities.size(); i++) {
            int tempRegion = cities.get(i).getRegionNumber();       
            if(!regions.contains(tempRegion)) {                 // If the region is not yet included                 
                regions.add(tempRegion);                        // Then add it
            }
        }

        return regions;
    }


    /**
     * Outputs the network as a solution to the terminal
    **/
    public void outputNetwork(ArrayList<City> citiesList) {

        System.out.println();
        System.out.println("The score of the network is: " + this.score);
        System.out.println("The total profit of the network is: " + this.totalProfit);
        System.out.println("The total emission resulting from the network is: " + this.totalEmission);
        System.out.println("The number of regions in the network is: " + this.numRegions);

        System.out.println();

        double truckCapUsed = 0;
        double shipCapUsed = 0;
        double planeCapUsed = 0;
        double trainCapUsed = 0;

        System.out.println("The following connections are included in the network:");
        System.out.println();

        for(int i = 0; i < this.networkEdges.length; i++) {
            for(int j = 0; j < this.networkEdges.length; j++) {
                for(int k = 0; k < 4; k++) {
                    if(this.networkEdges[i][j][k] > 0) {
                        if(k == 0) {
                            System.out.println(citiesList.get(i).getCityName() + " to " + citiesList.get(j).getCityName() + " by Truck,     using truck capacity of: " + this.networkEdges[i][j][k]);
                            truckCapUsed += this.networkEdges[i][j][k]; 
                        } else if(k == 1) {
                            System.out.println(citiesList.get(i).getCityName() + " to " + citiesList.get(j).getCityName() + " by Ship,     using sea capacity of: " + this.networkEdges[i][j][k]);
                            shipCapUsed += this.networkEdges[i][j][k];
                        } else if(k == 2) {
                            System.out.println(citiesList.get(i).getCityName() + " to " + citiesList.get(j).getCityName() + " by Airplane,     using air capacity of: " + this.networkEdges[i][j][k]);
                            planeCapUsed += this.networkEdges[i][j][k];
                        } else if(k == 3) {
                            System.out.println(citiesList.get(i).getCityName() + " to " + citiesList.get(j).getCityName() + " by Train,     using rail capacity of: " + this.networkEdges[i][j][k]);
                            trainCapUsed += this.networkEdges[i][j][k];
                        }
                        
                        
                    }
                }
            }
        }

        System.out.println();
        System.out.println("The total capacity used for Truck is: " + truckCapUsed);
        System.out.println("The total capacity used for Ship is: " + shipCapUsed);
        System.out.println("The total capacity used for Plane is: " + planeCapUsed);
        System.out.println("The total capacity used for Train is: " + trainCapUsed);
        System.out.println();
    }
}
