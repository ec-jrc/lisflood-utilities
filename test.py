class ParentClass:
    def __init__(self, message):
        self.message = message

    def filter(self):
        print("ParentClass filter")
        self.apply_filter()

    def apply_filter(self):
        print("ParentClass apply_filter:", self.message)

class ChildClass(ParentClass):
    def __init__(self, message):
        super().__init__(message)

    def apply_filter(self):
        print("ChildClass apply_filter:", self.message)

class_name = 'ChildClass'
childs = [globals()[class_name]('Hello, I am a child')]
for child in childs:
    child.filter()

import math

def m_to_decimal_degrees(m, latitude):
    # Convert meters to kilometers
    km = m / 1000

    # Convert km to degrees of latitude
    delta_lat = km / 111.32

    # Convert km to degrees of longitude at the given latitude
    delta_long = km / (111.32 * math.cos(math.radians(latitude)))

    return delta_lat, delta_long

latitude = 0  # Assuming the latitude is 0 (equator) for simplicity
delta_lat, delta_long = m_to_decimal_degrees(100000, latitude)
print("Change in latitude:", delta_lat)
print("Change in longitude:", delta_long)

# import math

def haversine(lon1, lat1, lon2, lat2):
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a))
    print(f'dlat {dlat} dlon {dlon} c {c}')
    r = 6371  # Radius of the Earth in kilometers
    distance = c * r
    return distance

def get_radius_in_arcminutes(lon1, lat1, lon2, lat2):
    distance_km = haversine(lon1, lat1, lon2, lat2)
    # Convert distance from kilometers to arcminutes
    distance_arcminutes = distance_km * 60
    radius_arcminutes = distance_arcminutes / 2
    return radius_arcminutes

# Example usage
lon1, lat1 = 39.8972, 47.8125  # Coordinates of point 1 in decimal degrees
lon2, lat2 = 41.1355, 47.47555833  # Coordinates of point 2 in decimal degrees
distance_km = haversine(lon1, lat1, lon2, lat2)
print("The distance between the two points is", distance_km, "kilometers")
radius_arcminutes = get_radius_in_arcminutes(lon1, lat1, lon2, lat2)
print("The radius of the circle containing the two points is approximately", radius_arcminutes, "arcminutes")

def get_radius(lon1, lat1, lon2, lat2):
    # Convert decimal degrees to radians
    # lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    dlon = abs(lon2 - lon1)
    dlat = abs(lat2 - lat1)
    r = math.sqrt(dlat*dlat + dlon*dlon)
    print(f'dlat {dlat} dlon {dlon} r {r}')
    return r
radius = get_radius(lon1, lat1, lon2, lat2)
print("The radius of the circle containing the two points is approximately", radius)
