"""
This script generates and processes geospatial data for the city of Graz, Austria.
It performs the following tasks:

1. **Generate Evenly Distributed Coordinates**:
   - The function `generate_evenly_distributed_coordinates` creates a grid of `n_points` latitude-longitude pairs 
     within a defined bounding box for Graz.

2. **Convert Coordinates to Addresses**:
   - The function `coordinates_to_address` uses the Google Maps Geocoding API to retrieve address information 
     for each generated coordinate.

3. **Plot Addresses on an Interactive Map**:
   - The function `plot_addresses_on_map` visualizes the geocoded addresses on a Folium map.
   - Each point is marked with its corresponding address and the map is saved as an interactive HTML file.

4. **Save Address Points as a Shapefile**:
   - The function `save_points_to_shapefile` creates a GeoDataFrame from the geocoded points.
   - The dataset includes address information, latitude, longitude, and geometry.
   - The data is stored in a `.shp` file, which can be used in GIS applications like QGIS or ArcGIS.
"""

import folium
import geopandas as gpd
from shapely.geometry import Point
import numpy as np

def generate_evenly_distributed_coordinates(lat_min, lat_max, lon_min, lon_max, n_points):
    # Generates n_points coordinates evenly distributed within the given bounds.
    # Calculate an approximate grid size based on n_points
    grid_size = int(np.ceil(np.sqrt(n_points)))  # Ensure the grid has at least n_points
    latitudes = np.linspace(lat_min, lat_max, grid_size)
    longitudes = np.linspace(lon_min, lon_max, grid_size)
    
    # Create a grid of coordinates
    coordinates = [
        (lat, lon) for lat in latitudes for lon in longitudes
    ]
    
    # Check if enough points are available in the grid
    if len(coordinates) < n_points:
        raise ValueError("Not enough points in the grid to sample the requested number of points.")

    # Randomly sample n_points from the grid
    selected_indices = np.random.choice(len(coordinates), n_points, replace=False)
    return [coordinates[i] for i in selected_indices]


def coordinates_to_address(gmaps, coordinates):
    # Coordinates to addresses using Google Maps Geocoding API.
    address_points = []
    for lat, lon in coordinates:
        try:
            result = gmaps.reverse_geocode((lat, lon))
            if result:
                address = result[0]["formatted_address"]
                # Check if the address is in Graz
                if "Graz" in address:
                    address_points.append({"address": address, "lat": lat, "lon": lon})
        except Exception as e:
            print(f"Error fetching address for {lat}, {lon}: {e}")
    return address_points

def plot_addresses_on_map(address_points, map_center=(47.0707, 15.4395), zoom_start=13):
    #Plots the addresses on a map using folium.
    # Initialize the map
    map_graz = folium.Map(location=map_center, zoom_start=zoom_start)

    # Add each address as a marker
    for point in address_points:
        folium.Marker(
            location=(point["lat"], point["lon"]),
            popup=point["address"],
        ).add_to(map_graz)

    # Save the map to an HTML file and display it
    map_graz.save("graz_addresses_map.html")
    print("Map saved as 'graz_addresses_map.html'")

def save_points_to_shapefile(address_points, filename="graz_points.shp"):
    #Saves the address points as a Shapefile.
    # Create a GeoDataFrame from the address points
    geometry = [Point(point["lon"], point["lat"]) for point in address_points]
    data = {
        "Address": [point["address"] for point in address_points],
        "Latitude": [point["lat"] for point in address_points],
        "Longitude": [point["lon"] for point in address_points],
    }
    gdf = gpd.GeoDataFrame(data, geometry=geometry, crs="EPSG:4326")

    # Save to a Shapefile
    gdf.to_file(filename)
    print(f"Shapefile saved as '{filename}'")

# if __name__ == "__main__":
#     google_api_key = "ENTER-KEY"  # Set your Google API key here
#     gmaps = googlemaps.Client(key=google_api_key)

#     # Define the bounding box for Graz (approximate values)
#     lat_min, lat_max = 47.020, 47.120  # Latitude range for Graz
#     lon_min, lon_max = 15.380, 15.500  # Longitude range for Graz

#     # Generate evenly distributed coordinates
#     n_points = 1500  # Number of points to generate
#     random_coords = generate_evenly_distributed_coordinates(lat_min, lat_max, lon_min, lon_max, n_points)

#     # Convert coordinates to addresses
#     address_points = coordinates_to_address(gmaps, random_coords)

#     # Plot the addresses on a map
#     plot_addresses_on_map(address_points)

#     # Save the points to a Shapefile
#     save_points_to_shapefile(address_points)
