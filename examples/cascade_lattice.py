import os
import subprocess


def main():
    # Example usage:
    quad_order = 19
    cl_fine_vals = [0.02, 0.01, 0.005]
    cl_coarse = 0.075
    cl_boundary = 0.2
    length = 0.075

    for i in range(3):
        cl_fine = cl_fine_vals[i]
        modify_and_run_geo_file(
            cl_fine=cl_fine, cl_coarse=cl_coarse, cl_boundary=cl_boundary, length=length
        )
        modify_quad_order(quad_order=quad_order)
        new_log_file, log_dir = modify_and_save_log_file(
            quad_order, cl_fine, cl_coarse, cl_boundary, length
        )
        print(f"New log file name: {new_log_file}")
        print(f"Log directory: {log_dir}")
        run_cpp_program(config_file_path="lattice.cfg")

    return 0


def modify_and_run_geo_file(cl_fine, cl_coarse, cl_boundary, length):
    # Define the file path
    geo_file_path = "meshes/geo_files/lattice_refined.geo"

    # Read the contents of the .geo file
    with open(geo_file_path, "r") as file:
        lines = file.readlines()

    # Modify the variables in the first 4 lines
    lines[0] = f"cl_fine = {cl_fine};\n"
    lines[1] = f"cl_coarse = {cl_coarse};\n"
    lines[2] = f"cl_boundary = {cl_boundary};\n"
    lines[3] = f"length = {length};\n"

    # Save the modified .geo file
    with open(geo_file_path, "w") as file:
        file.writelines(lines)

    # Execute the command
    command = f"gmsh {geo_file_path} -2 -format su2 -save_all"
    os.system(command)


def modify_quad_order(quad_order):
    # Define the file path
    cfg_file_path = "lattice.cfg"

    # Read the contents of the lattice.cfg file
    with open(cfg_file_path, "r") as file:
        lines = file.readlines()

    # Find the line containing "QUAD_ORDER" and modify it
    for i, line in enumerate(lines):
        if "QUAD_ORDER" in line:
            lines[i] = f"QUAD_ORDER = {quad_order}\n"
            break

    # Save the modified lattice.cfg file
    with open(cfg_file_path, "w") as file:
        file.writelines(lines)


def modify_and_save_log_file(quad_order, cl_fine, cl_coarse, cl_boundary, length):
    # Define the file path
    cfg_file_path = "lattice.cfg"

    # Read the contents of the lattice.cfg file
    with open(cfg_file_path, "r") as file:
        lines = file.readlines()

    # Identify the lines containing "LOG_DIR" and "LOG_FILE"
    log_dir_line = next(line for line in lines if "LOG_DIR" in line)
    log_file_line = next(line for line in lines if "LOG_FILE" in line)

    print(log_file_line)
    # Extract the current log directory and file name
    current_log_dir = log_dir_line.split("=")[1].strip()
    current_log_file = log_file_line.split("=")[1].strip()

    # Construct the new log file name
    new_log_file = f"lattice_quad_order_{quad_order}_cl_fine_{cl_fine}_cl_coarse_{cl_coarse}_cl_boundary_{cl_boundary}_length_{length}"

    # Update the log file name in the lines list
    log_file_line_new = f"LOG_FILE = {new_log_file}\n"
    lines[lines.index(log_file_line)] = log_file_line_new

    # Save the modified lattice.cfg file
    with open(cfg_file_path, "w") as file:
        file.writelines(lines)

    # Return the new log file name and log directory
    return new_log_file, current_log_dir


def run_cpp_program(config_file_path):
    # Define the path to the C++ program
    cpp_program_path = "../build/KiT-RT"

    # Construct the command to run the C++ program with the specified configuration file
    command = [cpp_program_path, config_file_path]

    try:
        # Run the C++ program
        subprocess.run(command, check=True)
        print("C++ program executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running C++ program: {e}")


if __name__ == "__main__":
    main()
