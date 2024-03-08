import subprocess
import os
import pandas as pd
import umbridge
import datetime


class KiTRT(umbridge.Model):
    def __init__(self):
        super().__init__("KiT-RT")
        # self._model = GenzFunction()

    def get_input_sizes(self, config):
        return [config.get("nvars", 5)]

    def get_output_sizes(self, config):
        return [1]

    def __call__(self, parameters, config):
        # CFG needs to be created here
        quad_order_value = 20  # Replace with the desired value
        self.write_lattice_script("lattice.cfg", quad_order_value)

        # KiT-RT needs to be called here

        # Step 1: Run the C++ application
        cpp_application_path = "../../build/KiT-RT"
        cfg_file_path = "lattice.cfg"
        log_file_folder = "result/logs/"

        # Run the C++ application
        subprocess.run([cpp_application_path, cfg_file_path])

        # Step 2: Read in the logfile.csv as a DataFrame
        log_file_path = os.path.join(log_file_folder, "lattice_log.csv")

        # Check if the logfile exists
        if os.path.exists(log_file_path):
            # Read the CSV file into a DataFrame
            logfile_df = pd.read_csv(log_file_path)

            # Step 3: Identify columns with names in the 'qois' list
            qois = [
                "TOTAL_OUTFLOW",
                "MAX_OUTFLOW",
                "CUR_PARTICLE_ABSORPTION",
            ]  # Add your desired column names here

            # Filter columns based on 'qois'
            selected_columns_df = logfile_df[qois]

            # Display the selected columns DataFrame
            print("Selected Columns DataFrame:")
            print(selected_columns_df)
        else:
            print(f"Error: Log file not found at {log_file_path}")

        """ ==============
        sample = np.asarray(parameters).T
        assert sample.shape[0] == config.get("nvars", 5)
        self._model.set_coefficients(
            sample.shape[0],
            config.get("c_factor", 1),
            config.get("coef_type", "sqexp"),
            config.get("w_factor", 0.5),
        )
        name = config.get("name", "oscillatory")
        val = self._model(name, sample)[0, 0]
        """
        return [[0]]

    def supports_evaluate(self):
        return True

    def write_lattice_script(filename, quad_order_value):
        current_date = datetime.datetime.now().strftime("%d. %b. %Y")

        script_content = f"""\
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Lattice Benchmarking File SN        %
        %  Author UM-Bridge bechmarking tool   %
        %  Date   {current_date}               %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %
        % ----IO specification ----
        %
        OUTPUT_DIR = result
        OUTPUT_FILE = lattice
        LOG_DIR = result/logs
        LOG_FILE = lattice_log
        MESH_FILE = meshes/lattice_unstructured.su2
        %
        % --- Problem definition ---
        %
        PROBLEM = LATTICE
        TIME_FINAL = 10
        SPATIAL_DIM = 3
        SOURCE_MAGNITUDE = 10.0
        %
        % ---- Solver specifications ----
        %
        % Solver type
        SOLVER = SN_SOLVER
        % CFL number
        CFL_NUMBER = 0.45
        % Reconstruction order
        RECONS_ORDER = 1
        %
        % ---- Boundary Conditions ----
        %
        BC_NEUMANN = ( void )
        %
        % ----- Quadrature Specification ---
        %
        QUAD_TYPE = GAUSS_LEGENDRE_TENSORIZED
        QUAD_ORDER = {quad_order_value}
        %
        % ----- Output ---- 
        %
        VOLUME_OUTPUT = (MINIMAL)
        VOLUME_OUTPUT_FREQUENCY = 1
        SCREEN_OUTPUT = (ITER, CUR_OUTFLOW)
        % (ITER, MASS,RMS_FLUX, VTK_OUTPUT, CSV_OUTPUT, CUR_OUTFLOW, TOTAL_OUTFLOW, MAX_OUTFLOW, CUR_PARTICLE_ABSORPTION, TOTAL_PARTICLE_ABSORPTION, MAX_PARTICLE_ABSORPTION )
        SCREEN_OUTPUT_FREQUENCY = 1
        HISTORY_OUTPUT = (ITER, MASS, RMS_FLUX, VTK_OUTPUT, CSV_OUTPUT, TOTAL_OUTFLOW, MAX_OUTFLOW, CUR_PARTICLE_ABSORPTION, TOTAL_PARTICLE_ABSORPTION, MAX_PARTICLE_ABSORPTION)
        HISTORY_OUTPUT_FREQUENCY = 1
        """

        with open(filename, "w") as f:
            f.write(script_content)


print("here")
print(umbridge.supported_models("http://0.0.0.0:80"))

models = [KiTRT()]
print("here")
umbridge.serve_models(models, 80)
