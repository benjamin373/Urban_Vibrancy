"""
This script utilizes Google Street View images and YOLOv11 object detection to analyze urban areas by detecting 
the number of people in specific locations in Graz, Austria. It performs the following tasks:

1. **Load Geographic Data**:
   - Reads a Shapefile (`graz_points.shp`) containing predefined geographic points in Graz.
   - Converts coordinates to `EPSG:4326` for compatibility with mapping services.

2. **Retrieve Street View Panoramas**:
   - The function `find_and_merge_panorama` searches for the nearest available Google Street View panorama 
     at each location.
   - If no panorama is found at the exact coordinates, the function expands the search radius up to 300 meters.
   - If a panorama is available, it downloads four images (headings: `0°, 90°, 180°, 270°`).
   - These images are merged to create a complete 360° panorama.
   - The final panorama is saved as a `.jpg` file.

3. **Detect People Using YOLOv11**:
   - The function `count_people_in_image` loads the pre-trained YOLOv11 model (`yolo11m.pt`) to analyze the panoramas.
   - It detects and counts the number of people in the image.
   - The function also calculates the average confidence score of the detections.

4. **Store and Display Results**:
   - Each analyzed point is assigned the number of detected people.
   - The data is stored in an updated Shapefile (`graz_points_with_people.shp`).
   - The results are visualized using Matplotlib:
     - The locations of analyzed points are plotted on a map.
     - The number of detected people is displayed for each location.
     - Individual panoramas are shown with detection results.
"""
#Alter Name der Datei: points_übergeben_und_analyse
import geopandas as gpd
import requests
from ultralytics import YOLO
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
import os

GOOGLE_API_KEY = "ENTER-KEY"


model = YOLO("yolo11m.pt")

# Function to recognize people using YOLO
def count_people_in_image(image_path):
    results = model(image_path)
    person_count = 0
    confidences = []
    for result in results:
        for detection in result.boxes:
            if detection.cls == 0:  # Class 0 corresponds to "person"
                person_count += 1
                confidences.append(detection.conf)

    average_confidence = sum(confidences) / len(confidences) if confidences else 0
    return person_count, average_confidence

# Function to download a panorama at the nearest available address
def find_and_merge_panorama(lat, lon, output_folder):
    pic_base = 'https://maps.googleapis.com/maps/api/streetview?'
    metadata_base = 'https://maps.googleapis.com/maps/api/streetview/metadata?'

    search_radius = 0  # Start radius
    max_radius = 300  # Maximum search in meters

    while search_radius <= max_radius:
        metadata_params = {
            'key': GOOGLE_API_KEY,
            'location': f"{lat},{lon}",
            'radius': search_radius
        }

        metadata_response = requests.get(metadata_base, params=metadata_params)
        metadata = metadata_response.json()

        if metadata_response.status_code == 200 and metadata.get('status') == 'OK':
            # Available panorama found
            pano_lat, pano_lon = metadata['location']['lat'], metadata['location']['lng']
            print(f"Available panorama found at {pano_lat}, {pano_lon} with radius {search_radius}m.")

            # Define common parameters
            pic_params = {
                'key': GOOGLE_API_KEY,
                'location': f"{pano_lat},{pano_lon}",
                'size': "640x640",  
                'fov': 90,          
                'pitch': 0          
            }

            headings = [0, 90, 180, 270]  # 360° coverage (4 directions)
            images = []

            if not os.path.exists(output_folder):
                os.makedirs(output_folder)

            # Download images for each direction
            for heading in headings:
                pic_params['heading'] = heading
                pic_response = requests.get(pic_base, params=pic_params)

                if pic_response.status_code == 200:
                    img = Image.open(BytesIO(pic_response.content))
                    images.append(img)
                    print(f"Image at viewing direction {heading}° downloaded.")
                else:
                    print(f"Error retrieving image for view direction {heading}: {pic_response.status_code}")

            # Merge images
            if images:
                total_width = sum(img.width for img in images)
                max_height = max(img.height for img in images)

                panorama = Image.new('RGB', (total_width, max_height))

                current_x = 0
                for img in images:
                    panorama.paste(img, (current_x, 0))
                    current_x += img.width

                # Save panorama
                panorama_path = os.path.join(output_folder, f"panorama_{pano_lat}_{pano_lon}.jpg")
                panorama.save(panorama_path)
                print(f"Panorama saved as {panorama_path}")

                return panorama_path

        search_radius += 10  # Increase radius by 10 meters

    print(f"No available panorama for the coordinates ({lat}, {lon}).")
    return None

# Main function
def main():
    # Load the points from the uploaded file
    points_gdf = gpd.read_file("D:/Countpeople/Panorama_test/graz_points.shp")
    points_gdf = points_gdf.to_crs("EPSG:4326")

    person_counts = []
    panorama_paths = []
    output_folder = "/mnt/data/panoramas"

    for i, point in enumerate(points_gdf.geometry):
        lat, lon = point.y, point.x
        panorama_path = find_and_merge_panorama(lat, lon, output_folder)

        if panorama_path:
            panorama_paths.append(panorama_path)
            person_count, _ = count_people_in_image(panorama_path)
            person_counts.append(person_count)

            # View panorama
            img = Image.open(panorama_path)
            plt.figure(figsize=(16, 8))
            plt.imshow(img)
            plt.title(f"Points {i + 1}: {person_count} People detected")
            plt.axis('off')
            plt.show()
        else:
            print(f"No panorama available for point {i + 1}.")
            person_counts.append(0)

#     # Add the number of people to the attributes
#     points_gdf['person_count'] = person_counts
#     points_gdf.to_file("D:/Countpeople/Panorama_test/graz_points_with_people.shp")

#     # Visualize the results
#     fig, ax = plt.subplots(figsize=(12, 12))
#     points_gdf.plot(ax=ax, color="red", markersize=10, label="Punkte", zorder=5)
#     ax.set_title("Points in Graz with number of people")
#     ax.set_xlabel("Longitude")
#     ax.set_ylabel("Latitude")

#     plt.legend()
#     plt.show()

# if __name__ == "__main__":
#     main()
