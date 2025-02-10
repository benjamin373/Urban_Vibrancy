
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup
import folium
import webbrowser
import requests 
import geopandas as gpd
import h3
import branca.colormap as cm
from collections import Counter
import datetime
import osmnx as ox
import matplotlib.pyplot as plt
import json
from collections import defaultdict
from shapely.geometry import Polygon
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib.patches import FancyArrow
import pandas as pd
from matplotlib.lines import Line2D


#user_date_input = input("Please enter a date in the format DD.MM.YYYY: ")
#user_time_input = input("Please enter the time in HH:MM format:")
#user_city_input = input("Please enter a city: ")

user_date_input = "18.11.2024"
user_time_input = "14:24"
user_city_input ="Graz"
# Define H3 resolution
h3_resolution = 8

# Validate city
user_city_input = user_city_input.strip()
if not user_city_input:
    print("Invalid city name.")
    exit()

# Parse and validate dates
try:
    user_date = datetime.datetime.strptime(user_date_input, "%d.%m.%Y").date()
except ValueError:
    print("Invalid date.")
    exit()

# Determine the season
if datetime.date(user_date.year, 6, 21) <= user_date <= datetime.date(user_date.year, 9, 23):
    # Summer
    Tree_multi = 0.5
    Poi_multi = 0.4
    Paths_multi = 0.5
elif datetime.date(user_date.year, 9, 23) < user_date <= datetime.date(user_date.year, 12, 21):
    # Fall
    Tree_multi = 0.4
    Poi_multi = 0.4
    Paths_multi = 0.4
elif datetime.date(user_date.year, 12, 21) <= user_date or user_date < datetime.date(user_date.year, 3, 20):
    # Winter
    Tree_multi = 0.1
    Poi_multi = 0.6
    Paths_multi = 0.3
elif datetime.date(user_date.year, 3, 20) <= user_date < datetime.date(user_date.year, 6, 21):
    # Spring
    Tree_multi = 0.3
    Poi_multi = 0.5
    Paths_multi = 0.4
else:
    print("Invalid date.")
    exit()

# Debugging: Saison- und Multiplikatorwerte ausgeben
print(f"Date: {user_date}, Season: {'Summer' if Tree_multi == 0.5 else 'Fall' if Tree_multi == 0.4 else 'Winter' if Tree_multi == 0.1 else 'Spring'}")
print(f"Tree_multi: {Tree_multi}, Poi_multi: {Poi_multi}")
print(f"Analyzing city: {user_city_input}")



overpass_url = "http://overpass-api.de/api/interpreter"

def get_data(query):
    response = requests.get(overpass_url, params={'data': query})
    if response.status_code == 200:
        data = response.json()
        if not data.get('elements'):
            print(f"No data for {user_city_input} found.")
            exit()
        return data
    else:
        print(f"Error retrieving data: {response.status_code}")
        exit()


def get_city_coordinates(city_name):
    city_query = f"""
    [out:json];
    area["name"="{city_name}"]->.searchArea;
    node(area.searchArea)["place"="city"];
    out center 1;
    """
    city_data = get_data(city_query)
    if 'elements' in city_data and city_data['elements']:
        city_center = city_data['elements'][0]
        return city_center['lat'], city_center['lon']
    else:
        print(f"No coordinates for {city_name} found.")
        exit()
        
# Get city coordinates
city_coordinates = get_city_coordinates(user_city_input)


user_date = datetime.datetime.strptime(user_date_input, "%d.%m.%Y").date()


# Get weather information from Open-Meteo API
def get_weather(city_coordinates, user_date, user_time):
    base_url = "https://api.open-meteo.com/v1/forecast"
    params = {
        "latitude": city_coordinates[0],
        "longitude": city_coordinates[1],
        "hourly": "weathercode",
        "timezone": "auto",
        "start_date": user_date.strftime("%Y-%m-%d"),
        "end_date": user_date.strftime("%Y-%m-%d")
    }

    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        weather_data = response.json()
        if "hourly" in weather_data and "weathercode" in weather_data["hourly"]:
            hourly_weather = weather_data["hourly"]["weathercode"]
            times = weather_data["hourly"]["time"]

            # Find time that is closest to the specified time
            datetime_strs = [datetime.datetime.fromisoformat(t) for t in times]
            user_datetime = datetime.datetime.combine(user_date, datetime.datetime.strptime(user_time, "%H:%M").time())
            closest_time_index = min(range(len(datetime_strs)), key=lambda i: abs(datetime_strs[i] - user_datetime))
            
            return hourly_weather[closest_time_index]
        else:
            print("Error: No weather data available.")
            return None
    else:
        print(f"Error retrieving weather data: {response.status_code}")
        return None


# Determine weather multiplier
def calculate_weather_multiplier(weather_code):
    # Open-Meteo Weathercodes: https://open-meteo.com/en/docs
    if weather_code in [0]:  # Clear sky
        return 1.2
    elif weather_code in [1, 2, 3]:  # Partly or completely cloudy
        return 0.8
    elif weather_code in [61, 63, 65, 80, 81, 82]:  # Rain
        return 0.3
    elif weather_code in [95, 96, 99]:  # Thunderstorm
        return 0.2
    else:
        return 1  # Default value for other weather conditions


weather_code = get_weather(city_coordinates, user_date, user_time_input)


# Weather-Multiplikator berechnen
if weather_code is not None:
    weather_multi = calculate_weather_multiplier(weather_code)
    print(f"Weather-Multiplikator: {weather_multi} (Code: {weather_code})")
else:
    print("Weatherdata could not be retrieved. Default weather multiplier value is used.")
    weather_multi = 1
    
# Determine the day of the week
weekday = user_date.strftime("%A")
print(f"Weekday: {weekday}")


#OSM queries
footway_query = f"""
[out:json];
area["name"="{user_city_input}"]->.searchArea;
(
  way["highway"="footway"](area.searchArea);
);
out center;
"""

tree_query = f"""
[out:json];
area["name"="{user_city_input}"]->.searchArea;
(
  node["natural"="tree"](area.searchArea);
  way["natural"="tree"](area.searchArea);
  relation["natural"="tree"](area.searchArea);
);
out center;
"""

poi_query = f"""
[out:json];
area["name"="{user_city_input}"]->.searchArea;
(
  node["amenity"~"cafe|restaurant|bar|pub|nightclub"](area.searchArea);
  way["amenity"~"cafe|restaurant|bar|pub|nightclub"](area.searchArea);
  relation["amenity"~"cafe|restaurant|bar|pub|nightclub"](area.searchArea);

  node["amenity"~"bus_station|taxi|bus_stop"](area.searchArea);
  way["amenity"~"bus_station|taxi|bus_stop"](area.searchArea);
  relation["amenity"~"bus_station|taxi|bus_stop"](area.searchArea);

  node["public_transport"="stop_position"](area.searchArea);
  way["public_transport"="stop_position"](area.searchArea);
  relation["public_transport"="stop_position"](area.searchArea);

  node["highway"="bus_stop"](area.searchArea);
  way["highway"="bus_stop"](area.searchArea);
  relation["highway"="bus_stop"](area.searchArea);

  node["amenity"~"toilets|waste_basket|bench"](area.searchArea);
  way["amenity"~"toilets|waste_basket|bench"](area.searchArea);
  relation["amenity"~"toilets|waste_basket|bench"](area.searchArea);

  node["amenity"~"school|university|kindergarten"](area.searchArea);
  way["amenity"~"school|university|kindergarten"](area.searchArea);
  relation["amenity"~"school|university|kindergarten"](area.searchArea);

  node["shop"](area.searchArea);
  way["shop"](area.searchArea);
  relation["shop"](area.searchArea);

  node["amenity"~"marketplace|supermarket"](area.searchArea);
  way["amenity"~"marketplace|supermarket"](area.searchArea);
  relation["amenity"~"marketplace|supermarket"](area.searchArea);

  node["amenity"~"fast_food"](area.searchArea);
  way["amenity"~"fast_food"](area.searchArea);
  relation["amenity"~"fast_food"](area.searchArea);

  node["amenity"="bicycle_parking"](area.searchArea);
  way["amenity"="bicycle_parking"](area.searchArea);
  relation["amenity"="bicycle_parking"](area.searchArea);

  way["building:levels"](if:t["building:levels"] > 5)(area.searchArea);
  relation["building:levels"](if:t["building:levels"] > 5)(area.searchArea);
);
out center;
"""


# Retrieve footpath data
footway_data = get_data(footway_query)

# Retrieve tree data
tree_data = get_data(tree_query)

# Retrieve poi data
poi_data = get_data(poi_query)


# Create GeoDataFrame for footpaths
footway_elements = footway_data['elements']
footways = []
for element in footway_elements:
    if 'center' in element:
        footways.append({'latitude': element['center']['lat'], 'longitude': element['center']['lon']})

gdf_footways = gpd.GeoDataFrame(
    footways, 
    geometry=gpd.points_from_xy([f['longitude'] for f in footways], [f['latitude'] for f in footways]), 
    crs="EPSG:4326"
)

# Create GeoDataFrame for trees
tree_elements = tree_data['elements']
trees = []
for element in tree_elements:
    if 'lat' in element and 'lon' in element:
        trees.append({'latitude': element['lat'], 'longitude': element['lon']})

gdf_trees = gpd.GeoDataFrame(
    trees, 
    geometry=gpd.points_from_xy([t['longitude'] for t in trees], [t['latitude'] for t in trees]), 
    crs="EPSG:4326"
)

# Create GeoDataFrame for Poi
poi_elements = poi_data['elements']
poi = []
for element in poi_elements:
    if 'lat' in element and 'lon' in element:
        poi.append({'latitude': element['lat'], 'longitude': element['lon']})

gdf_poi = gpd.GeoDataFrame(
    poi, 
    geometry=gpd.points_from_xy([p['longitude'] for p in poi], [p['latitude'] for p in poi]), 
    crs="EPSG:4326"
)


# Calculate H3 cells for each footpath
gdf_footways['h3_index'] = gdf_footways.geometry.apply(lambda x: h3.latlng_to_cell(x.y, x.x, h3_resolution))

# Calculate H3 cells for each tree
gdf_trees['h3_index'] = gdf_trees.geometry.apply(lambda x: h3.latlng_to_cell(x.y, x.x, h3_resolution))

# Calculate H3 cells for each poi
gdf_poi['h3_index'] = gdf_poi.geometry.apply(lambda x: h3.latlng_to_cell(x.y, x.x, h3_resolution))

# Count number of footpaths per hexagon
footway_hex_counts = Counter(gdf_footways['h3_index'])

# Count number of trees per hexagon
tree_hex_counts = Counter(gdf_trees['h3_index'])

# Count number of poi per hexagon
poi_hex_counts = Counter(gdf_poi['h3_index'])



# Path to GeoJSON file with detected people from Street View
geojson_path = r"D:\Countpeople\Panorama_test\points_1500_geojson.geojson"

people_hex_counts = defaultdict(int)  # For the new people_count assignment


# Load GeoJSON data and extract `people_count`
with open(geojson_path, "r", encoding="utf-8") as f:
    geojson_data = json.load(f)

# Iterate through the features of the GeoJSON file
for feature in geojson_data["features"]:
    # Extract the coordinates and the `people_count`
    coords = feature["geometry"]["coordinates"]  # [lon, lat]
    people_count = feature["properties"].get("people_count", 0)
    
    # Determine the H3 hexagon for this coordinate
    h3_index = h3.latlng_to_cell(coords[1], coords[0],h3_resolution)
    
    # Add the people_count to this hexagon
    people_hex_counts[h3_index] += people_count

# First loop: Calculate and collect combined values
combined_values = {}

for h3_index in set(footway_hex_counts.keys()).union(tree_hex_counts.keys(), poi_hex_counts.keys(), people_hex_counts.keys()):
    footways_count = footway_hex_counts.get(h3_index, 0)  # Number of footpaths in the hexagon
    trees_count = tree_hex_counts.get(h3_index, 0)       # Number of trees in the hexagon
    poi_count = poi_hex_counts.get(h3_index, 0)          # Number of POIs in the hexagon
    person_count = people_hex_counts.get(h3_index, 0)    # Number of people in the hexagon

    # Calculate combined density
    combined_value = (
        (trees_count * Tree_multi) +
        (footways_count * Paths_multi) +
        (poi_count * Poi_multi) +
        (person_count)  
    ) * weather_multi  # Apply weather multiplier

    combined_values[h3_index] = combined_value

# Determine the minimum and maximum values ​​of the combined values
min_combined_value = min(combined_values.values())
max_combined_value = max(combined_values.values())

# Second loop: Calculate and store normalized values
combined_hex_counts = {}

for h3_index, combined_value in combined_values.items():
    if max_combined_value == min_combined_value:
        # If all values ​​are equal, set to 1
        combined_hex_counts[h3_index] = 1
    else:
        # Min-max normalization
        normalized_value = (combined_value - min_combined_value) / (max_combined_value - min_combined_value)
        combined_hex_counts[h3_index] = normalized_value
        

# Create map based on city coordinates
m = folium.Map(location=city_coordinates, zoom_start=14)

# FeatureGroup for the footpaths
footway_layer = folium.FeatureGroup(name="Footways")

# FeatureGroup for the trees
tree_layer = folium.FeatureGroup(name="Trees")

# FeatureGroup for the POIs
poi_layer = folium.FeatureGroup(name="poi")


# FeatureGroup for the combined hexagons
combined_layer = folium.FeatureGroup(name="Urban vibrancy")

# Color scale for footpaths (green to red)
footway_colormap = cm.LinearColormap(
    ['green', 'red'],
    vmin=min(footway_hex_counts.values()),
    vmax=max(footway_hex_counts.values())
)
footway_colormap.caption = 'Footways'

# Color scale for trees (green to red)
tree_colormap = cm.LinearColormap(
    ['green', 'red'],  
    vmin=min(tree_hex_counts.values()),
    vmax=max(tree_hex_counts.values())
)
tree_colormap.caption = 'Trees'


# Color scale for POIs
poi_colormap = cm.LinearColormap(
    ['green', 'red'],
    vmin=min(poi_hex_counts.values()),  
    vmax=max(poi_hex_counts.values())   
)
poi_colormap.caption = 'poi'

# Combined density color scale (green to red)
combined_colormap = cm.LinearColormap(
    ['green', 'red'],  
    vmin=min(combined_hex_counts.values()),
    vmax=max(combined_hex_counts.values())
)
combined_colormap.caption = 'Urban vibrancy (Trees * 0.2 + Footways * 0.8)'



# Draw hexagons for footpaths on the map
for h3_index, count in footway_hex_counts.items():
    boundary = h3.cell_to_boundary(h3_index)
    boundary = [(lat, lng) for lat, lng in boundary]
    tooltip = f'Footways: {count}'
    folium.Polygon(
        locations=boundary,
        color=footway_colormap(count),
        weight=0.5,
        fill=True,
        fill_opacity=0.5,
        tooltip=tooltip
    ).add_to(footway_layer)

# Draw hexagons for trees on the map
for h3_index, count in tree_hex_counts.items():
    boundary = h3.cell_to_boundary(h3_index)
    boundary = [(lat, lng) for lat, lng in boundary]
    tooltip = f'Trees: {count}'
    folium.Polygon(
        locations=boundary,
        color=tree_colormap(count),
        weight=0.5,
        fill=True,
        fill_opacity=0.5,
        tooltip=tooltip
    ).add_to(tree_layer)

# Draw hexagons for POIs on the map
for h3_index, count in poi_hex_counts.items():
    boundary = h3.cell_to_boundary(h3_index)
    boundary = [(lat, lng) for lat, lng in boundary]
    tooltip = f'poi: {count}'
    folium.Polygon(
        locations=boundary,
        color=poi_colormap(count),
        weight=0.5,
        fill=True,
        fill_opacity=0.5,
        tooltip=tooltip
    ).add_to(poi_layer)

# Draw hexagons for the combined density on the map
for h3_index, value in combined_hex_counts.items():
    boundary = h3.cell_to_boundary(h3_index)
    boundary = [(lat, lng) for lat, lng in boundary]
    tooltip = f'Urban vibrancy: {value:.2f}'
    folium.Polygon(
        locations=boundary,
        color=combined_colormap(value),
        weight=0.5,
        fill=True,
        fill_opacity=0.5,
        tooltip=tooltip
    ).add_to(combined_layer)



# Grey Hexagone
# city limits from OSM
graz_boundary = ox.geocode_to_gdf(user_city_input)

# Get the geometry of the city
geo = graz_boundary.geometry[0]

# Convert the geometry to H3 cells
h3_cells = h3.geo_to_cells(geo, res=h3_resolution)


# Iterate through all calculated H3 cells
for cell in h3_cells:
    boundary = h3.cell_to_boundary(cell)

# Collect all used H3 cells (Footways, Trees, POIs, Urban Vibrancy)
used_h3_cells = set(footway_hex_counts.keys()) | set(tree_hex_counts.keys()) | set(poi_hex_counts.keys()) | set(combined_hex_counts.keys())

# Find all unused H3 cells
unused_h3_cells = set(h3_cells) - used_h3_cells

# Draw the unused H3 cells as gray hexagons (no data)
for cell in unused_h3_cells:
    boundary = h3.cell_to_boundary(cell)
    boundary = [(lat, lng) for lat, lng in boundary]
    folium.Polygon(
        locations=boundary,
        color="grey",
        weight=0.5,
        fill=True,
        fill_opacity=0.4,
        tooltip="No Data"
    ).add_to(m)




# old RSS part (must be redone with downloaded data)

from datetime import datetime, timedelta

user_time = datetime.strptime(user_time_input, "%H:%M").time()

def extract_time_from_event_page(event_page_url):
    try:
        # Abrufen der Event-Seite
        event_response = requests.get(event_page_url)
        if event_response.status_code == 200:
            soup_event = BeautifulSoup(event_response.text, 'html.parser')

            # Suche nach dem <div>, das die Uhrzeit enthält
            time_div = soup_event.find('div', {'data-eym': True})
            if time_div:
                # Extrahiere den Text aus dem <div>
                time_text = time_div.get_text(strip=True)

                # Suche nach einem Muster für die Uhrzeit (z. B. "19:00")
                import re
                match = re.search(r'\b\d{1,2}:\d{2}\b', time_text)
                if match:
                    return match.group()  # Gibt die gefundene Uhrzeit zurück, z. B. "19:00"
        else:
            print(f"Fehler beim Abrufen der Event-Seite: {event_page_url}")
    except Exception as e:
        print(f"Fehler beim Extrahieren der Uhrzeit: {e}")
    return None

# RSS-Feed URL
rss_url = "https://kultur.graz.at/ksv2_heute.xml"

# RSS-Feed abrufen
response = requests.get(rss_url)

if response.status_code == 200:
    xml_content = response.text
else:
    print(f"Fehler beim Abrufen des Feeds: {response.status_code}")
    xml_content = ""

if xml_content:
    root = ET.fromstring(xml_content)

    namespaces = {
        'rdf': "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        '': "http://purl.org/rss/1.0/"
    }

    library = []

    for item in root.findall('item', namespaces):
        name = item.find('title', namespaces).text if item.find('title', namespaces) is not None else ""
        link = item.find('link', namespaces).text if item.find('link', namespaces) is not None else ""

        description = item.find('description', namespaces)
        href = ""
        if description is not None and 'href=' in description.text:
            start = description.text.find('href="') + len('href="')
            end = description.text.find('"', start)
            href = description.text[start:end]

        # Uhrzeit aus der Eventseite extrahieren
        event_time_str = extract_time_from_event_page(href)
        event_time = None
        if event_time_str:
            try:
                event_time = datetime.strptime(event_time_str, "%H:%M").time()
            except ValueError:
                pass

        if event_time:
            library.append({
                "name": name.strip() if name else "",
                "link": link.strip() if link else "",
                "href": href.strip() if href else "",
                "time": event_time
            })

    # Funktion zur Geokodierung der Adresse mit OpenStreetMap (Nominatim)
    def geocode_address(address):
        geocode_url = f"https://nominatim.openstreetmap.org/search?format=json&q={address}"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        response = requests.get(geocode_url, headers=headers)
        if response.status_code == 200:
            result = response.json()
            if result:
                lat = float(result[0]['lat'])
                lng = float(result[0]['lon'])
                return lat, lng
        return None, None

    # Funktion zur Erstellung der Karte mit gefilterten Events. Fügt direkt m hinzu
    def create_map_for_filtered_events():
        for entry in library:
            event_time = entry['time']
            if event_time:
                # Zeitdifferenz berechnen
                event_time_delta = timedelta(hours=event_time.hour, minutes=event_time.minute)
                user_time_delta = timedelta(hours=user_time.hour, minutes=user_time.minute)
                time_difference = (event_time_delta - user_time_delta).total_seconds()

                # Event innerhalb von 2 Stunden nach user_time filtern
                if 0 <= time_difference <= 7200:  # Innerhalb von 2 Stunden
                    event_url = entry["href"]
                    if event_url:
                        response = requests.get(event_url)
                        if response.status_code == 200:
                            soup = BeautifulSoup(response.text, 'html.parser')
                            location_div = soup.find('div', id='ort')
                            if location_div:
                                location_link = location_div.find('a')
                                if location_link:
                                    location_name = location_link.text.strip()
                                    location_url = location_link['href']
                                    base_url = "https://kultur.graz.at"
                                    full_location_url = base_url + location_url
                                    location_response = requests.get(full_location_url)
                                    if location_response.status_code == 200:
                                        location_soup = BeautifulSoup(location_response.text, 'html.parser')
                                        address_div = location_soup.find('div', id='adresse')
                                        if address_div:
                                            address = address_div.find('b').text.strip()
                                            lat, lng = geocode_address(address)
                                            if lat and lng:
                                                folium.Marker(
                                                    [lat, lng],
                                                    popup=f"{location_name}\n{entry['time']}"
                                                ).add_to(m)

    

    create_map_for_filtered_events()
else:
    print("Keine Daten im RSS-Feed gefunden.")

# old RSS part (must be redone with downloaded data)


# folium map
# Add layers and color scales to the map
footway_layer.add_to(m)
tree_layer.add_to(m)
poi_layer.add_to(m)
combined_layer.add_to(m)
footway_colormap.add_to(m)
tree_colormap.add_to(m)
combined_colormap.add_to(m)
folium.Map(location=city_coordinates, zoom_start=12).add_to(m)

# Add layer control to map
folium.LayerControl().add_to(m)




# Save and view map
map_file = f"{user_city_input}_combined_map_footways.html"
m.save(map_file)
print(f"Combined map with footpaths saved as: {map_file}")
webbrowser.open(map_file)



# matplotlib map
# convert H3 cells to polygonsln
def h3_to_polygon(h3_index):
    boundary = h3.cell_to_boundary(h3_index)
    return Polygon([(lng, lat) for lat, lng in boundary])

# Create GeoDataFrame for combined density
combined_gdf = gpd.GeoDataFrame(
    {
        "h3_index": list(combined_hex_counts.keys()),
        "density": list(combined_hex_counts.values()),
    },
    geometry=[h3_to_polygon(h3_index) for h3_index in combined_hex_counts.keys()],
    crs="EPSG:4326",
)

# Create GeoDataFrame for gray hexagons (unused cells).
unused_gdf = gpd.GeoDataFrame(
    {
        "h3_index": list(unused_h3_cells),
        "density": [None] * len(unused_h3_cells),
    },
    geometry=[h3_to_polygon(h3_index) for h3_index in unused_h3_cells],
    crs="EPSG:4326",
)

# Combined GeoDataFrame (all cells: data + no data)
all_gdf = gpd.GeoDataFrame(
    pd.concat([combined_gdf, unused_gdf], ignore_index=True),
    crs="EPSG:4326",
)

# Save GeoDataFrames as GeoJSON
combined_gdf.to_file("combined_density_res11.geojson", driver="GeoJSON")
unused_gdf.to_file("unused_cells_res11.geojson", driver="GeoJSON")
all_gdf.to_file("all_cells_res11.geojson", driver="GeoJSON")

# Transform the GeoDataFrame into projected CRS (UTM 33N)
all_gdf = all_gdf.to_crs("EPSG:32633")

# Create Matplotlib plot
fig, ax = plt.subplots(figsize=(12, 12))

# Plot hexagons: "No Data" in gray and data with a gradient
all_gdf.plot(
    ax=ax,
    column='density',               
    cmap='RdYlGn_r',                
    missing_kwds={                  
        "color": "lightgrey",       
        "edgecolor": "none",        
        "label": "No Data"          
    },
    edgecolor='none',               
    legend=True,                    
    legend_kwds={"shrink": 0.5},    
)

# Add additional legend for "No Data".
no_data_patch = Line2D(
    [0], [0], marker='h', color='none', markerfacecolor='lightgrey', markersize=10, 
    label="No Data", linestyle='None', markeredgecolor='none', markeredgewidth=2
)

# Add north arrow
arrow_x, arrow_y = 0.06, 0.05 
arrow_dx, arrow_dy = 0, 0.05
ax.add_patch(
    FancyArrow(
        x=arrow_x, y=arrow_y, dx=arrow_dx, dy=arrow_dy, width=0.01, 
        color="black", transform=ax.transAxes
    )
)
ax.text(arrow_x, arrow_y + 0.06, "N", transform=ax.transAxes,
        fontsize=12, ha='center', va='center', color='black')

# Add scale
scalebar = ScaleBar(1, location='lower left', units='m', scale_loc='bottom', length_fraction=0.25)
ax.add_artist(scalebar)

# Legend for all elements
handles, labels = ax.get_legend_handles_labels()
handles.append(no_data_patch)
labels.append("No Data")
ax.legend(handles=handles, labels=labels, loc='lower right', bbox_to_anchor=(1.1, 0.1))
legend_kwds={"label": "Density Scale", "shrink": 0.5}

source_text = (
    "Sources:\n"
    "Google Street View\n"
    "OpenStreetMap\n"
    "Author: Chassin, Signitzer\n"
    "Projection: EPSG:32633"
)

# position of the text
ax.text(
    1.0, 0.005, source_text, transform=ax.transAxes, fontsize=10,
    ha='left', va='bottom', color='black', bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.5')
)

# Add buffers around the axes
ax.margins(0.1)

# Titel
fig.suptitle(f"Urban Vibrancy Map of {user_city_input}", fontsize=15, y=0.96)

# Subtitle
fig.text(0.5, 0.92, "Based on POIs, Streetview and activity data aggregated to hexagons", 
         ha='center', fontsize=10, color='gray')

ax.set_axis_off()

plt.tight_layout()
plt.show()

