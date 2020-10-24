import CyRSoXS as cy
import numpy as np
import pandas as pd


class CyRsoxsRunner:
    """Framework object to support running CyRSoXS

    Usage example:

    input_values = {
        "energies": [280.0, 285.0, 288.0],
        "physical_size": 5.0,
        "dimensions": (64, 64, 8),
        "angle_start": 0.0,
        "angle_end": 360.0,
        "angle_increment": 1.0,
        "interpolation_type": cy.InterpolationType.Linear,  # optional
        "windowing_type": cy.FFTWindowing.NoPadding,  # optional
    }

    runner = CyRsoxsRunner(input_values)

    material_filenames = ["material_1.txt", "material_2.txt"]
    runner.set_optical_constant_from_files(material_filenames)

    # create the objects to use
    voxels = runner.make_voxel_object()
    scattering = runner.make_scattering_object()

    # create the morphology and fill the voxel data
    morphology = get_morphology()  # user supplied function
    set_voxels(morphology, voxels)  # user supplied function

    # run the simulation on the voxels, sending result to scatering
    runner.run(voxels, scattering)

    """

    def __init__(self, values, require_optional=False,
                 print_values=True):
        """Initialize the class with input data items from given dict.

        The value keys supported are:
            "energies": list of energies to simulate (required)
            "physical_size": physical size of the system (required)
            "dimensions": number of voxels in x,y,z directions (required)
            "angle_start": starting angle for rotation (required)
            "angle_end": ending angle for rotation (required)
            "angle_increment": angle step for rotation (required)
            "omp_threads": number of openMP threads (optional)
            "rotation_mask": whether to use masking for rotation (optional)
            "interpolation_type": type of Ewalds interpolation (optional)
            "windowing_type": windowing type (optional)

        All other keys are silently ignored.

        Args:
            values: dict with values to set.
            require_optional: If True, code will raise a KeyError if any
                optional values are not provided.
            print_values: If True, prints input data values

        Raises:
            ValueError: The input data validation fails, despite
                having all required items
        """
        self.refractive_index = None
        self.input_data = cy.InputData()
        self.n_materials = cy.get_n_materials()

        # required values. If missing, an exception is raised
        try:
            self.input_data.setEnergies(values['energies'])
            # Not sure we can get the energies out of the input data
            # object, so just store them locally
            self.energies = values['energies'][:]
            self.input_data.physSize(values['physical_size'])
            self.input_data.dimensions(X=values['dimensions'][0],
                                       Y=values['dimensions'][1],
                                       Z=values['dimensions'][2])
            self.input_data.setRotationAngle(
                StartAngle=values['angle_start'],
                EndAngle=values['angle_end'],
                IncrementAngle=values['angle_increment'])
        except KeyError as e:
            print(f"Missing required input item when setting input: {str(e)}")
            raise e

        # optional items, if missing they are ignored unless
        # require_optional is True
        for function, key in [
            (self.input_data.interpolationType, 'interpolation_type'),
            (self.input_data.windowingType,'windowing_type'),
            (self.input_data.rotMask, 'rotation_mask'),
            (self.input_data.openMP, 'omp_threads'),
        ]:
            try:
                function = values[key]
            except KeyError as e:
                print(f"warning: missing {str(e)}")
                if require_optional:
                    raise e

        if not self.input_data.validate():
            raise ValueError('Input data validation failed')

        if print_values:
            self.print_input_data()

    def print_input_data(self):
        """Prints the current input data values"""
        self.input_data.print()

    def generate_data_frame(self, filename, labelEnergy=None, sep='\s+'):
        """Returns a data frame to use for energies.

        Args:
            filename: The path where the filename is located. (string)
            labelEnergy: A dict with label Energy for each file.
                This maps the names for each column to the column index.
                (uses default values if None)
            sep: Seperator used for the file. (string)
        """
        if labelEnergy is None:
            labelEnergy = {
                'BetaPara': 0,
                'BetaPerp': 1,
                'DeltaPara': 2,
                'DeltaPerp': 3,
                'Energy': 6,
            }

        column_names = [
            'DeltaPara', 'BetaPara', 'DeltaPerp', 'BetaPerp', 'Energy'
        ]
        indices = [labelEnergy[column_name] for column_name in column_names]

        try:
            EnergyFile = pd.read_csv(filename, sep)
        except FileNotFoundError as e:
            print(f"{str(e)}")
            return

        data_frame = EnergyFile.iloc[:, indices].copy()
        data_frame.columns = column_names
        data_frame.sort_values(by=['Energy'], inplace=True)
        data_frame.drop_duplicates(subset=['Energy'],
                                   keep=False,
                                   ignore_index=True,
                                   inplace=True)
        return data_frame

    def get_interpolated_energies(self, value, data_frame):
        """Returns linearly interpolated value after removing any duplicates.

        Args:
            value: Energy value at which to interpolate. (double/float)
            data_frame: Data frame of energy file
        
        Returns:
            A list of interpolated optical properties for the given energy. 
        """
        energy_id = 4
        nearest_id = data_frame['Energy'].sub(value).abs().idxmin()
        numColumns = len(data_frame.columns)

        valArray = np.zeros(numColumns);

        if data_frame.iloc[nearest_id][energy_id] > value:
            xp = [
                data_frame.iloc[nearest_id - 1][energy_id],
                data_frame.iloc[nearest_id][energy_id]
            ]
            for i in range(numColumns):
                yp = [
                    data_frame.iloc[nearest_id - 1][i],
                    data_frame.iloc[nearest_id][i]
                ]
                valArray[i] = np.interp(value, xp, yp)
        elif data_frame.iloc[nearest_id][energy_id] < value:
            xp = [
                data_frame.iloc[nearest_id][energy_id],
                data_frame.iloc[nearest_id + 1][energy_id]
            ]
            for i in range(numColumns):
                yp = [
                    data_frame.iloc[nearest_id][i],
                    data_frame.iloc[nearest_id + 1][i]
                ]
                valArray[i] = np.interp(value, xp, yp)
        else:
            for i in range(numColumns):
                valArray[i] = data_frame.iloc[nearest_id][i]

        return valArray[0:4].tolist()

    def set_optical_constant_from_files(self, material_files, print_values=False):
        """Sets the optical constants from the given files

        This creates and fills refractive_index with the energies from
        the input data and the data from the files given. The data is added
        in the order the files are given.

        Args:
            material_files: Files with optical data (list of filenames)
                The number of files must be the same as the number of
                materials.
            print_values: If True, prints the optical values after they
                have been set.
        """
        if not self.input_data.validate():
            print('[ERROR] input data must be valid to set optical constants')
            return

        if len(material_files) != self.n_materials:
            print('[ERROR] number of material files must match number of materials.')
            return

        if self.refractive_index is None:
            self.refractive_index = cy.RefractiveIndex(self.input_data)

        # create data_frames for each file
        frames = [self.generate_data_frame(fname) for fname in material_files]

        for energy in self.energies:
            vals = [self.get_interpolated_energies(energy, frame) 
                for frame in frames]
            self.refractive_index.addData(OpticalConstants=vals, Energy=energy)

        if not self.refractive_index.validate():
            raise ValueError('Refractive index validation failed')

        if print_values:
            self.print_optical_constants()
        
    def print_optical_constants(self):
        """Prints the optical constants."""
        self.refractive_index.print()

    def make_scattering_object(self):
        """Makes and returns an empty scattering pattern object."""
        if not self.input_data.validate():
            print('[ERROR] input data must be valid to make scattering object')
            return None
        scattering_pattern = cy.ScatteringPattern(self.input_data)
        return scattering_pattern

    def make_voxel_object(self):
        """Makes and returns an empty voxel object."""
        if not self.input_data.validate():
            print('[ERROR] input data must be valid to make voxel object')
            return None
        return cy.VoxelData(self.input_data)
    
    def run(self, voxel_data, scattering_pattern, stdout=True, stderr=True,
            is_batch=False, write_meta_data=True):
        """Run the simulation.

        Args:
            voxel_data: Files with optical data (list of filenames)
                The number of files must be the same as the number of
                materials.
            print_values: If True, prints the optical values after they
                have been set.
            stdout
            stderr:
        """
        with cy.ostream_redirect(stdout=stdout, stderr=stderr):
            cy.launch(VoxelData=voxel_data,
                      RefractiveIndexData=self.refractive_index,
                      InputData=self.input_data,
                      ScatteringPattern=scattering_pattern,
                      IsBatch=is_batch,
                      WriteMetaData=write_meta_data)
