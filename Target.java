public class Target {
    
    private int originCity;
    private int destinationCity;
    private int revenue;

    public Target(int originCity, int destinationCity, int revenue) {

        this.originCity = originCity;
        this.destinationCity = destinationCity;
        this.revenue = revenue;
    }


    /**
     * returns the city number of the origin city
     * @return int
    **/
    public int getOriginCity() {
        return this.originCity;
    }


    /**
     * returns the city number of the destination city
     * @return int
    **/
    public int getDestinationCity() {
        return this.destinationCity;
    }


    /**
     * returns the revenue in millions of euros from the connection made between origin city to destination city
     * @return int
    **/
    public int getRevenue() {
        return this.revenue;
    }
}
