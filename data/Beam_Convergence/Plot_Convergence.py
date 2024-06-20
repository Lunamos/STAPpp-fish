import numpy as np
import matplotlib.pyplot as plt

# Material Properties
E = 1000
Iz = 1/12
Iy =1/12
J = 0
L = 10

# Load Case
q_y = -10

# Exact Solution
def exact_Mz(x):
    return q_y * (x-L)**2 / 2

def exact_Qy(x):
    return q_y * (L - x)

def exact_deflection(x):
    return q_y * (x**4 - 4 * L * x**3 + 6 * L**2 * x**2) / (24 * E * Iz)

def exact_angle(x):
    return q_y * (x**3 - 3 * L * x**2 + 3 * L**2 * x) / (6 * E * Iz)


# Extract Data
def extract_data(displacement_data, force_data):
    v = displacement_data[:, 2]  # Y-DISPLACEMENT
    theta_z = displacement_data[:, 6]  # Z-ANGLE
    moment_z_1 = force_data[:, 5]  # MOMENT_Z_1
    moment_z_2 = force_data[:, 6]  # MOMENT_Z_2
    shear_y = force_data[:, 7]  # SHEAR_Y
    
    data_dict = {
        'v': v,
        'theta_z': theta_z,
        'moment_z_1': moment_z_1,
        'moment_z_2': moment_z_2,
        'shear_y': shear_y
    }
    
    return data_dict


# Function to calculate x-coordinates for moments(gaussian points)
def calculate_moment_x_coords(mesh_density, num_elements):
    x_coords = []
    delta_x = 10 / mesh_density
    for I in range(1, num_elements + 1):
        x1 = (delta_x * (I - 1)) + (delta_x / 2) - (delta_x / 2) * np.sqrt(3) / 3
        x2 = (delta_x * (I - 1)) + (delta_x / 2) + (delta_x / 2) * np.sqrt(3) / 3
        x_coords.extend([x1, x2])
    return np.array(x_coords)

# Data Set
# Data for Mesh_density = 1
displacement_data_1 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [2, 0.0, -150.0, 0.0, 0.0, 0.0, -20.0]
])


force_data_1 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, -311.004, -22.3291, -50.0, 0.0]
])


# Data for Mesh_density = 2
displacement_data_2 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [2, 0.0, -53.2, 0.0, 0.0, 0.0, -17.53],
    [3, 0.0, -150.3, 0.0, 0.0, 0.0, -20.06]
])


force_data_2 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, -400.42, -183.913, -75.0, 0.0],
    [2, 0.0, 0.0, 0.0, 0.0, -78.2511, -6.08227, -25.0, 0.0]
])


# Data for Mesh_density = 5
displacement_data_5 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [2, 0.0, -10.48, 0.0, 0.0, 0.0, -9.76],
    [3, 0.0, -36.48, 0.0, 0.0, 0.0, -15.68],
    [4, 0.0, -71.28, 0.0, 0.0, 0.0, -18.72],
    [5, 0.0, -110.08, 0.0, 0.0, 0.0, -19.84],
    [6, 0.0, -150.0, 0.0, 0.0, 0.0, -20.0]
])


force_data_5 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, -458.628, -354.705, -90.0, 0.0],
    [2, 0.0, 0.0, 0.0, 0.0, -287.081, -206.252, -70.0, 0.0],
    [3, 0.0, 0.0, 0.0, 0.0, -155.534, -97.7992, -50.0, 0.0],
    [4, 0.0, 0.0, 0.0, 0.0, -63.9872, -29.3462, -30.0, 0.0],
    [5, 0.0, 0.0, 0.0, 0.0, -12.4402, -0.893164, -10.0, 0.0]
])


# Data for Mesh_density = 10
displacement_data_10 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [2, 0.0, -2.805, 0.0, 0.0, 0.0, -5.42],
    [3, 0.0, -10.48, 0.0, 0.0, 0.0, -9.76],
    [4, 0.0, -22.005, 0.0, 0.0, 0.0, -13.14],
    [5, 0.0, -36.48, 0.0, 0.0, 0.0, -15.68],
    [6, 0.0, -53.125, 0.0, 0.0, 0.0, -17.5],
    [7, 0.0, -71.28, 0.0, 0.0, 0.0, -18.72],
    [8, 0.0, -90.405, 0.0, 0.0, 0.0, -19.46],
    [9, 0.0, -110.08, 0.0, 0.0, 0.0, -19.84],
    [10, 0.0, -130.005, 0.0, 0.0, 0.0, -19.98],
    [11, 0.0, -150.0, 0.0, 0.0, 0.0, -20.0]
])


force_data_10 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, -479.091, -424.243, -95.0, 0.0],
    [2, 0.0, 0.0, 0.0, 0.0, -386.204, -337.129, -85.0, 0.0],
    [3, 0.0, 0.0, 0.0, 0.0, -303.317, -260.016, -75.0, 0.0],
    [4, 0.0, 0.0, 0.0, 0.0, -230.431, -192.903, -65.0, 0.0],
    [5, 0.0, 0.0, 0.0, 0.0, -167.544, -135.79, -55.0, 0.0],
    [6, 0.0, 0.0, 0.0, 0.0, -114.657, -88.6763, -45.0, 0.0],
    [7, 0.0, 0.0, 0.0, 0.0, -71.7703, -51.563, -35.0, 0.0],
    [8, 0.0, 0.0, 0.0, 0.0, -38.8835, -24.4498, -25.0, 0.0],
    [9, 0.0, 0.0, 0.0, 0.0, -15.9968, -7.33654, -15.0, 0.0],
    [10, 0.0, 0.0, 0.0, 0.0, -3.11004, -0.223291, -5.0, 0.0]
])


# Data for Mesh_density = 20
displacement_data_20 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [2, 0.0, -0.725313, 0.0, 0.0, 0.0, -2.8525],
    [3, 0.0, -2.805, 0.0, 0.0, 0.0, -5.42],
    [4, 0.0, -6.10031, 0.0, 0.0, 0.0, -7.7175],
    [5, 0.0, -10.48, 0.0, 0.0, 0.0, -9.76],
    [6, 0.0, -15.8203, 0.0, 0.0, 0.0, -11.5625],
    [7, 0.0, -22.005, 0.0, 0.0, 0.0, -13.14],
    [8, 0.0, -28.9253, 0.0, 0.0, 0.0, -14.5075],
    [9, 0.0, -36.48, 0.0, 0.0, 0.0, -15.68],
    [10, 0.0, -44.5753, 0.0, 0.0, 0.0, -16.6725],
    [11, 0.0, -53.125, 0.0, 0.0, 0.0, -17.5],
    [12, 0.0, -62.0503, 0.0, 0.0, 0.0, -18.1775],
    [13, 0.0, -71.28, 0.0, 0.0, 0.0, -18.72],
    [14, 0.0, -80.7503, 0.0, 0.0, 0.0, -19.1425],
    [15, 0.0, -90.405, 0.0, 0.0, 0.0, -19.46],
    [16, 0.0, -100.195, 0.0, 0.0, 0.0, -19.6875],
    [17, 0.0, -110.08, 0.0, 0.0, 0.0, -19.84],
    [18, 0.0, -120.025, 0.0, 0.0, 0.0, -19.9325],
    [19, 0.0, -130.005, 0.0, 0.0, 0.0, -19.98],
    [20, 0.0, -140.0, 0.0, 0.0, 0.0, -19.9975],
    [21, 0.0, -150.0, 0.0, 0.0, 0.0, -20.0]
])


force_data_20 = np.array([
    [1, 0.0, 0.0, 0.0, 0.0, -489.49, -461.344, -97.5, 0.0],
    [2, 0.0, 0.0, 0.0, 0.0, -441.268, -414.565, -92.5, 0.0],
    [3, 0.0, 0.0, 0.0, 0.0, -395.546, -370.287, -87.5, 0.0],
    [4, 0.0, 0.0, 0.0, 0.0, -352.325, -328.509, -82.5, 0.0],
    [5, 0.0, 0.0, 0.0, 0.0, -311.603, -289.231, -77.5, 0.0],
    [6, 0.0, 0.0, 0.0, 0.0, -273.381, -252.452, -72.5, 0.0],
    [7, 0.0, 0.0, 0.0, 0.0, -237.659, -218.174, -67.5, 0.0],
    [8, 0.0, 0.0, 0.0, 0.0, -204.438, -186.396, -62.5, 0.0],
    [9, 0.0, 0.0, 0.0, 0.0, -173.716, -157.117, -57.5, 0.0],
    [10, 0.0, 0.0, 0.0, 0.0, -145.494, -130.339, -52.5, 0.0],
    [11, 0.0, 0.0, 0.0, 0.0, -119.773, -106.061, -47.5, 0.0],
    [12, 0.0, 0.0, 0.0, 0.0, -96.551, -84.2823, -42.5, 0.0],
    [13, 0.0, 0.0, 0.0, 0.0, -75.8293, -65.004, -37.5, 0.0],
    [14, 0.0, 0.0, 0.0, 0.0, -57.6076, -48.2257, -32.5, 0.0],
    [15, 0.0, 0.0, 0.0, 0.0, -41.8859, -33.9474, -27.5, 0.0],
    [16, 0.0, 0.0, 0.0, 0.0, -28.6643, -22.1691, -22.5, 0.0],
    [17, 0.0, 0.0, 0.0, 0.0, -17.9426, -12.8908, -17.5, 0.0],
    [18, 0.0, 0.0, 0.0, 0.0, -9.72089, -6.11245, -12.5, 0.0],
    [19, 0.0, 0.0, 0.0, 0.0, -3.9992, -1.83413, -7.5, 0.0],
    [20, 0.0, 0.0, 0.0, 0.0, -0.777511, -0.0558228, -2.5, 0.0]
])


# Extracting data for all mesh densities
all_data = {
    'Mesh_density=1': extract_data(displacement_data_1, force_data_1),
    'Mesh_density=2': extract_data(displacement_data_2, force_data_2),
    'Mesh_density=5': extract_data(displacement_data_5, force_data_5),
    'Mesh_density=10': extract_data(displacement_data_10, force_data_10),
    'Mesh_density=20': extract_data(displacement_data_20, force_data_20)
}



# Assuming `all_data` is already populated with the extracted data
mesh_densities = {
    'Mesh_density=1': 1,
    'Mesh_density=2': 2,
    'Mesh_density=5': 5,
    'Mesh_density=10': 10,
    'Mesh_density=20': 20
}


'''
# Plot deflection (v) for different mesh densities  
plt.figure(figsize=(12, 8))

for key, density in mesh_densities.items():
    data = all_data[key]
    num_nodes = len(data['v'])
    x_coords = np.linspace(0, 10, num_nodes)  # x coordinates based on mesh density
    plt.plot(x_coords, data['v'], label=key)
# Plot exact solution
x_exact = np.linspace(0, 10, 100)
v_exact_values = exact_deflection(x_exact)
plt.plot(x_exact, v_exact_values, 'k--', label='Exact solution (v)')

plt.xlabel('X Coordinate')
plt.ylabel('Y-Deflection (v)')
plt.title('Y-Deflection for Different Mesh Densities')
plt.legend()
plt.grid(True)
plt.show()


# Plot Z-Angle for different mesh densities
plt.figure(figsize=(12, 8))

for key, density in mesh_densities.items():
    data = all_data[key]
    num_nodes = len(data['theta_z'])
    x_coords = np.linspace(0, 10, num_nodes)  # x coordinates based on mesh density
    plt.plot(x_coords, data['theta_z'], label=key)

# Plot exact solution
theta_z_exact_values = exact_angle(x_exact)
plt.plot(x_exact, theta_z_exact_values, 'k--', label='Exact solution (theta_z)')


plt.xlabel('X Coordinate')
plt.ylabel('Z-Angle (theta_z)')
plt.title('Z-Angle for Different Mesh Densities')
plt.legend()
plt.grid(True)
plt.show()



# Plotting Mz (Moment)
plt.figure(figsize=(12, 8))

for key, density in mesh_densities.items():
    data = all_data[key]
    num_elements = len(data['moment_z_1'])  # Number of elements
    x_coords = calculate_moment_x_coords(density, num_elements)
    moment_z = np.empty(2 * num_elements)
    moment_z[0::2] = data['moment_z_1']
    moment_z[1::2] = data['moment_z_2']

    plt.plot(x_coords, moment_z, label=f"{key} (numerical)")

# Plot exact solution
x_exact = np.linspace(0, 10, 100)
Mz_exact_values = exact_Mz(x_exact)
plt.plot(x_exact, Mz_exact_values, 'k--', label='Exact solution (Mz)')

plt.xlabel('X Coordinate')
plt.ylabel('Moment (Mz)')
plt.title('Moment (Mz) for Different Mesh Densities')
plt.legend()
plt.grid(True)
plt.show()



# Plot Qy
plt.figure(figsize=(12, 8))

for key, density in mesh_densities.items():
    data = all_data[key]
    num_elements = len(data['shear_y'])
    x_coords = np.linspace(0, 10, num_elements)  # x coordinates based on mesh density
    plt.plot(x_coords, data['shear_y'], label=key)

# Plot exact solution
Qy_exact_values = exact_Qy(x_exact)
plt.plot(x_exact, Qy_exact_values, 'k--', label='Exact solution (Qy)')

plt.xlabel('X Coordinate')
plt.ylabel('SHEAR_Y')
plt.title('SHEAR_Y for Different Mesh Densities')
plt.legend()
plt.grid(True)
plt.show()
'''

'''TO BE FIXED

# Calculate L2 error of deflection (v) for different mesh densities
def calculate_deflection_error(data, mesh_densities):
    # Calculate exact solution
    L2_error = np.zeros(len(mesh_densities))
    i = 0
    for key, density in mesh_densities.items():
        data = all_data[key]
        num_nodes = len(data['v'])

        # Calculate exact solution
        v_exact_values = exact_deflection(x_coords)
        L2_error[i] = np.sqrt(np.sum((data['v'] - v_exact_values)**2))
        i = i+1
    return L2_error

deflection_error = calculate_deflection_error(all_data, mesh_densities)
print(deflection_error)
'''