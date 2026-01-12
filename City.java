public class City {
    

    private int cityNumber;
    private String cityName;
    private double latitude;
    private double longitude;
    private int regionNumber;

    public City(int cityNumber, String cityName, double latitude, double longitude, int regionNumber) {

        this.cityNumber =  cityNumber;
        this.cityName = cityName;
        this.latitude = latitude;
        this.longitude = longitude;
        this.regionNumber = regionNumber;
    }


    /**
     * returns city number of the city 
     * @return int
    **/ 
    public int getCityNumber() {
        return this.cityNumber;
    }

    /**
     * returns the city name as a string
     * @return String
    **/
    public String getCityName() {
        return this.cityName;
    }


    /**
     * returns the latitude of the city
     * @return double
    **/
    public double getLatitude() {
        return this.latitude;
    }

    /**
     * returns the longitude of the city
     * @return double
    **/
    public double getLongitude() {
        return this.longitude;
    }


    /**
     * returns the region number of the city
     * @return int
    **/
    public int getRegionNumber() {
        return this.regionNumber;
    }
}
