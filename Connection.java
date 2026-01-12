public class Connection {
    
    private int originCity;
    private int destinationCity;
    private int costPerTruck;
    private int costPerSea;
    private int costPerAir;
    private int costPerRail;

    public Connection(int originCity, int destinationCity, int costPerTruck, int costPerSea, int costPerAir, int costPerRail) {

        this.originCity = originCity;
        this.destinationCity = destinationCity;
        this.costPerTruck = costPerTruck;
        this.costPerSea = costPerSea;
        this.costPerAir = costPerAir;
        this.costPerRail = costPerRail;
    }

    /**
     * returns city number of the origin city
     * @return int
    **/
    public int getOriginCity() {
        return this.originCity;
    }


    /**
     * returns city number of the destination city
     * @return int
    **/
    public int getDestinationCity() {
        return this.destinationCity;
    }


    /**
     * returns the cost per truck for travel from origin city to destination city
     * @return int
    **/
    public int getCostPerTruck() {
        return this.costPerTruck;
    }

    
    /**
     * returns the cost per sea for travel from origin city to destination city
     * @return int
    **/
    public int getCostPerSea() {
        return this.costPerSea;
    }


    /**
     * returns the cost per air for travel from origin city to destination city
     * @return int
    **/
    public int getCostPerAir() {
        return this.costPerAir;
    }


    /**
     * returns the cost per rail for travel from origin city to destination city
     * @return int
    **/
    public int getCostPerRail() {
        return this.costPerRail;
    }


    /**
     * This method returns a sorted array, where the array is sorted from transport that used the lowest number of transportation units
     * to transport that used the most. Transport is represented as follows: Truck = 0, Sea = 1, Air = 2, Rail = 3
     * If modalitiy does not exist (i.e. costPer.. = 0), then it is not included in the array
     * @return Integer Array
    **/
    public int[] getMinTrans() {

        int[] transCosts = {this.costPerTruck, this.costPerSea, this.costPerAir, this.costPerRail};        // saves the costs in an array with the pre specified order
        int count = 0;                          // count that will save the number of costs that are greater than zero
        for(int i = 0; i < transCosts.length; i++) {
            if(transCosts[i] > 0) {
                count++;
            }
        }
        
        int[] sortTransport = new int[count];       // This array will save the sorted modes of transport as specified before

        // This for loop will sort the modes of transport into the sortTransport
        for(int j = 0; j < count; j++) {
            int minCost = Integer.MAX_VALUE;
            int index = 0;
            for(int l = 0; l < transCosts.length; l++) {
                if((transCosts[l] > 0) && (transCosts[l] < minCost)) {
                    minCost = transCosts[l];
                    index = l;
                }
            }
            sortTransport[j] = index;   // index is what represents the mode of transport
            transCosts[index] = 0;      // set the found modality to cost of 0, so it will not be considered anymore
        }
        
        return sortTransport;
    }
    
}
