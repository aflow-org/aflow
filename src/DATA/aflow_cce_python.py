import json
import os
import subprocess
import warnings

VERBOSE = False


class CCE:
    def __init__(self, aflow_executable='aflow'):
        '''Class implementing CCE python functionality

        Methods:
        get_corrections -- get CCE corrections
        Arguments:
        struct_file_path -- path to the structure file
        functionals -- functionals to get corrections for (optional)
        enthalpies_formation_dft -- input DFT formation enthalpies (optional)
        oxidation_numbers -- input oxidation numbers for each atom (optional)

        get_oxidation_numbers -- get CCE oxidation numbers
        Arguments:
        struct_file_path -- path to the structure file

        get_cation_coordination_numbers -- get CCE cation coordination numbers
        Arguments:
        struct_file_path -- path to the structure file
        '''

        self.aflow_executable = aflow_executable
        # checks whether aflow executable exists
        try:
            subprocess.check_output(self.aflow_executable, shell=True)
        except subprocess.CalledProcessError:
            raise OSError('aflow executable not found: ' +
                          self.aflow_executable)

    def run_aflow_command(self, cmd):
        '''checks output of aflow command'''
        try:
            if VERBOSE:
                print('aflow_command(): cmd=' + self.aflow_executable + cmd)
            output = subprocess.check_output(self.aflow_executable + cmd,
                                             shell=True)
            return output
        except subprocess.CalledProcessError:
            raise OSError('An error occurred while executing: ' +
                          self.aflow_executable + cmd)

    def add_structure_file_path(self, cmd, struct_file):
        '''adds path to the structure file to the command'''
        if (os.path.exists(struct_file)):
            cmd += ' < ' + struct_file
            return cmd
        else:
            raise OSError(struct_file + ' not found')

    def get_json(self, cmd):
        '''sets json as output format, executes aflow command,
        and checks and returns json'''
        # return json output by default
        cmd += ' --print=json'
        # execute command
        output = self.run_aflow_command(cmd)
        # validate json
        try:
            res_json = json.loads(output)
        except json.decoder.JSONDecodeError:
            print(output)
            raise RuntimeError('Output is not valid json.')
        return res_json

    def validate_functionals(self, input_functionals):
        '''checks whether for all given functionals corrections are
        available'''
        available_functionals = ['PBE', 'LDA', 'SCAN']
        invalid_functionals = [n for n in input_functionals if n not in
                               available_functionals]
        valid_functionals = [n for n in input_functionals if n in
                             available_functionals]
        if invalid_functionals:
            warnings.warn('No corrections available for: ' +
                          ','.join(invalid_functionals))
        if not valid_functionals:
            raise ValueError('No valid functionals provided. '
                             'Please choose from: ' +
                             ','.join(available_functionals))
        return valid_functionals

    def add_functionals(self, functionals, cmd):
        '''adds functionals to command'''
        functionals_str = ','.join(functionals)
        cmd += ' --functionals=' + functionals_str
        return cmd

    def get_corrections(self, struct_file, enthalpies_formation_dft=[],
                        functionals=[], oxidation_numbers=[]):
        '''determines CCE corrections and corrected formation enthalpies if
        precalculated DFT values are given'''
        # initial cce command
        command = ' --get_cce_corrections'

        # adding path to structure file
        command = self.add_structure_file_path(command, struct_file)

        # handling functionals
        if type(functionals) == str:
            if functionals:
                # convert to list
                functionals = functionals.replace(' ', '')
                functionals = functionals.split(',')
        elif type(functionals) == tuple:
            if len(functionals) > 0:
                # convert to list
                functionals = list(functionals)
        elif type(functionals) == list:
            pass
        else:
            raise TypeError('The functionals argument must be either a list, '
                            'a tuple, or a string with functionals '
                            'separated by commas.')
        if len(functionals) > 0:
            functionals = [x.upper() for x in functionals]
            valid_functionals = self.validate_functionals(functionals)
            command = self.add_functionals(valid_functionals, command)

        # handling enthalpies_formation_dft
        if type(enthalpies_formation_dft) == str:
            if enthalpies_formation_dft:
                # convert to list
                enthalpies_formation_dft = enthalpies_formation_dft.replace(
                                           ' ', '')
                enthalpies_formation_dft = enthalpies_formation_dft.split(',')
        elif type(enthalpies_formation_dft) == tuple:
            if len(enthalpies_formation_dft) > 0:
                # convert to list
                enthalpies_formation_dft = list(enthalpies_formation_dft)
        elif type(enthalpies_formation_dft) == list:
            pass
        else:
            raise TypeError('The enthalpies argument must be either a list, '
                            'a tuple, or a string with enthalpies '
                            'separated by commas.')
        if len(enthalpies_formation_dft) > 0:
            if len(enthalpies_formation_dft) != len(functionals):
                raise ValueError('number of input formation enthalpies '
                                 'does not equal number of functionals')
            # consider only enthalpies corresponding to valid functionals
            if len(enthalpies_formation_dft) != len(valid_functionals):
                enthalpies_new = []
                for i in range(len(functionals)):
                    for j in range(len(valid_functionals)):
                        if functionals[i] == valid_functionals[j]:
                            enthalpies_new.append(enthalpies_formation_dft[i])
                enthalpies_formation_dft = enthalpies_new
            enthalpies_formation = ','.join([str(n) for n in
                                            enthalpies_formation_dft])
            command += ' --enthalpies_formation_dft=' + enthalpies_formation

        # handling oxidation_numbers
        if type(oxidation_numbers) == str:
            if oxidation_numbers:
                # convert to list
                oxidation_numbers = oxidation_numbers.replace(' ', '')
                oxidation_numbers = oxidation_numbers.split(',')
        elif type(oxidation_numbers) == tuple:
            if len(oxidation_numbers) > 0:
                # convert to list
                oxidation_numbers = list(oxidation_numbers)
        elif type(oxidation_numbers) == list:
            pass
        else:
            raise TypeError('The oxidation numbers must be either a list, '
                            'a tuple, or a string with numbers '
                            'separated by commas.')
        if len(oxidation_numbers) > 0:
            input_oxidation_numbers = ','.join([str(n) for n in
                                               oxidation_numbers])
            command += ' --oxidation_numbers=' + input_oxidation_numbers

        # get json
        return self.get_json(command)

    def get_oxidation_numbers(self, struct_file):
        '''determines oxidation numbers'''
        # initial cce command
        command = ' --get_oxidation_numbers'

        # adding path to structure file
        command = self.add_structure_file_path(command, struct_file)

        # get json
        return self.get_json(command)

    def get_cation_coordination_numbers(self, struct_file):
        '''determines cation coordination numbers'''
        # initial cce command
        command = ' --get_cation_coordination_numbers'

        # adding path to structure file
        command = self.add_structure_file_path(command, struct_file)

        # get json
        return self.get_json(command)
