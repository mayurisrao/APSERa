import numpy as np
import sys
from pathlib import Path

print("Enter filename without extension")
filename = sys.argv[1]+"_hdr"

# Initialize variables
system_date = None
start_time = None
end_time = None

# Open and read the file
with open(filename, 'r') as file:
    for line in file:
        if "System date:" in line:
            system_date = line.split(":")[1].strip()
        elif "Start time:" in line:
            start_time = line.split(":")[1].strip()
        elif "End time:" in line:
            end_time = line.split(":")[1].strip()

# Convert strings to NumPy variables (if needed)
# For example, you might want to convert 'system_date' to a NumPy datetime object
# You can use np.datetime64 if the date format is compatible
# Example: system_date_np = np.datetime64(system_date)

# Print the extracted values
print("System date:", system_date)
print("Start time:", start_time)
print("End time:", end_time)
